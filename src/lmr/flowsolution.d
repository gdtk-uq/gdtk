/** flowsolution.d
 * Eilmer compressible-flow simulation code, postprocessing functions.
 *
 * The role of the post-processing functions is just to pick up data
 * from a previously-run simulation and either write plotting files
 * or extract interesting pieces of data.  To do this, we really don't
 * need or want all of the code machinery in the classes built for the
 * simulation code so we rebuild a bit of the core data handling here.
 * This allows us to be a bit flexible in what variables are held.
 * We might want to add flow variables, such as Mach number or Pitot
 * pressure, whish are not normally held in the main simulation data
 * structures.  It also frees us of the internal numbering of cells in
 * the simulation code that must allow for ghost-cells.
 *
 * Author: Peter J. and Rowan G.
 * First code: 2015-06-09
 */

module flowsolution;

import std.stdio;
import std.conv;
import std.format;
import std.string;
import std.regex;
import std.algorithm;
import std.array;
import std.math;
import std.json;
import std.file;
import gzip;
import fileutil;

import ntypes.complex;
import nm.number;
import geom;
import gas;
import gas.vib_specific_nitrogen;
import globalconfig;
import json_helper;
import flowstate;
import vtk_writer;
import lmrconfig;
import blockio;

import util.lua;
import geom.luawrap;
import gas.luagas_model;
import luaflowstate;
import luaflowsolution;
import luaidealgasflow;
import luagasflow;
import geom.misc.kdtree;


class FlowSolution {
    // The collection of flow blocks and grid blocks that define the flow
    // and the domain at one particular instant in time.
public:
    double simTime;
    size_t nBlocks;
    GridMotion grid_motion;
    FluidBlockLite[] flowBlocks;
    Grid[] gridBlocks;

    this(int snapshot, size_t nBlocks, double simTime=-1.0, string dir=".", bool make_kdtree=true)
    {
        string cfgFile = (dir == ".") ? lmrCfg.cfgFile : dir ~ "/" ~ lmrCfg.cfgFile;
        string content;
        try {
            content = readText(cfgFile);
        } catch (Exception e) {
            writeln("Message is: ", e.msg);
            throw new Error(text("Failed to read config file: ", cfgFile));
        }
        JSONValue jsonData;
        try {
            jsonData = parseJSON!string(content);
        } catch (Exception e) {
            writeln("Message is: ", e.msg);
            throw new Error(text("Failed to parse JSON from config file: ", cfgFile));
        }
        string fieldFmt = jsonData["field_format"].str;
        GlobalConfig.field_format = fieldFmt;
        string gridFmt = jsonData["grid_format"].str;
        GlobalConfig.grid_format = gridFmt;
        grid_motion = grid_motion_from_name(jsonData["grid_motion"].str);
        // -- end initialising from JSONData

	    // Find out variables from metadata file
        string fluidMetadataFile = (dir == ".") ? lmrCfg.fluidMetadataFile : dir ~ "/" ~ lmrCfg.fluidMetadataFile;
        auto variables = readFluidVariablesFromFluidMetadata(fluidMetadataFile);
        //
        // Use <block-list-filename> to get a hint of the type of each block.
        string blkListFile = (dir == ".") ? lmrCfg.blkListFile : dir ~ "/" ~ lmrCfg.blkListFile;
        auto listFile = File(blkListFile);
        auto listFileLine = listFile.readln().chomp(); // only comments on the first line
        //
        foreach (ib; 0 .. nBlocks) {
            listFileLine = listFile.readln().chomp();
            int ib_chk;
            string gridTypeName;
            string label;
            int ncells;
            formattedRead(listFileLine, " %d %s %s %d", &ib_chk, &gridTypeName, &label, &ncells);
            if (ib != ib_chk) {
                string msg = format("Reading %s file ib=%d ib_chk=%d", blkListFile, ib, ib_chk);
                throw new FlowSolverException(msg);
            }
            auto gridType = gridTypeFromName(gridTypeName);
            string gName;
            if (grid_motion != GridMotion.none) {
                gName = gridFilename(snapshot, to!int(ib));
            }
            else {
                gName = gridFilename(lmrCfg.initialFieldDir, to!int(ib));
            }
            gName = (dir == ".") ? gName : dir ~ "/" ~ gName;
            final switch (gridType) {
            case Grid_t.structured_grid:
                gridBlocks ~= new StructuredGrid(gName, gridFmt);
                break;
            case Grid_t.unstructured_grid:
                gridBlocks ~= new UnstructuredGrid(gName, gridFmt, false);
            }
            gridBlocks[$-1].sort_cells_into_bins();
            string fName = fluidFilename(snapshot, to!int(ib));
            fName = (dir == ".") ? fName : dir ~ "/" ~ fName;
            auto fb = new FluidBlockLite(fName, ib, simTime, jsonData, gridType, fieldFmt, variables, ncells);
            final switch (gridType) {
            case Grid_t.structured_grid:
                auto g = gridBlocks[$-1];
                fb.nic = g.niv - 1;
                fb.njc = g.njv - 1;
                fb.nkc = max(g.nkv-1, 1);
                break;
            case Grid_t.unstructured_grid:
                fb.nic = ncells;
                fb.njc = 1;
                fb.nkc = 1;
            }
            flowBlocks ~= fb;
        } // end foreach ib
        this.nBlocks = nBlocks;
        this.simTime = simTime;

        if (make_kdtree) {
            construct_kdtree();
        }
    } // end constructor


    // This method is used to free up memory. The object shell still remains,
    // but all of its internals have been deallocated. You should not attempt
    // to use the object after calling this method. This method is designed
    // to be called via the Lua wrapper. It is for callers who know what they
    // are doing with the objects underneath.
    void releaseMemory()
    {
        destroy(flowBlocks);
        destroy(gridBlocks);
    }

    override string toString()
    {
        string str = "FlowSolution(";
        str ~= "nBlocks=" ~ to!string(nBlocks);
        str ~= ", variables=" ~ to!string(flowBlocks[0].variableNames);
        str ~= ")";
        return str;
    }

    void add_aux_variables(string[] addVarsList, string tag="field")
    {
        foreach (blk; flowBlocks) blk.add_aux_variables(addVarsList, tag);
    }

    size_t[] find_enclosing_cell(ref const(Vector3) p)
    {
        size_t[] cell_identity = [0, 0, 0]; // blk_id, i, found_flag

        // If the user has opted into the fast, but memory intensive kdtree approach
        if (kdtree_made) {
            Node cellnode = {[p.x.re, p.y.re, p.z.re]};
            const(Node)* found = null;
            double bestDist = 0.0;
            size_t nVisited = 0;
            root.fast_nearest(cellnode, 0, found, bestDist, nVisited);
            cell_identity = [found.blkid, found.cellid, 1];
        }
        else { // Otherwise we do the standard check for an enclosing cell in each block
            foreach (ib; 0 .. nBlocks) {
                bool found = false;
                size_t indx = 0;
                gridBlocks[ib].find_enclosing_cell(p, indx, found);
                if (found) {
                    cell_identity = [ib, indx, (found)?1:0];
                    break;
                }
            } // foreach ib
        }
        return cell_identity;
    } // end find_enclosing_cell()

    size_t[] find_enclosing_cell(double x, double y, double z=0.0)
    {
        Vector3 p = Vector3(x, y, z);
        return find_enclosing_cell(p);
    }

    size_t find_enclosing_cells_along_line(ref const(Vector3) p0, ref const(Vector3) p1,
                                           size_t n, ref size_t[2][] cells_found)
    // Locate cells along the line p0 --> p1, accumulating the block index and the cell index
    // into the cells_found array.
    {
        size_t count = 0;
        foreach (ip; 0 .. n) {
            double frac = double(ip) / double(n-1);
            Vector3 p = p0*(1.0-frac) + p1*frac;
            auto identity = find_enclosing_cell(p);
            size_t ib = identity[0]; size_t idx = identity[1];
            size_t found = identity[2];
            if (found == 0) { // out of domain bounds
                writeln("# Info: Cell not found for point ", p);
                continue;
            } else { // maybe store cell data
                // It is convenient to omit repeated cells so that we can specify
                // tiny steps and be sure of sampling all cells and also be able
                // to specify multiple lines that get concatenated, again without
                // repeated cells in the output stream.
                if ((cells_found.length == 0) ||
                    (cells_found[$-1][0] != ib) || (cells_found[$-1][1] != idx)) {
                    cells_found ~= [ib, idx]; // Add a "new" cell.
                    count += 1;
                }
            }
        } // end foreach ip
        return count;
    } // end find_enclosing_cells_along_line()

    size_t[] find_nearest_cell_centre(double x, double y, double z=0.0)
    {
        // [TODO] Think about delegating this work to the underlying grid objects.
        // It seems that we should not be repeating code here, but it also seems that
        // we should not be repeating the cell positions in the flow solution file
        // when they can be obtained (via computation) from the underlying grid object.
        // The counter argument is that having the cell centre position available in
        // the solution file has been handy, over time.
        // PJ, 2016-11-13
        bool withZpos = canFind(flowBlocks[0].variableNames, "pos.z");
        size_t[] nearest = [0, 0]; // blk_id, i
        double dx = x - flowBlocks[0]["pos.x", 0];
        double dy = y - flowBlocks[0]["pos.y", 0];
        double dz = (withZpos) ? (z - flowBlocks[0]["pos.z", 0]) : 0.0;
        double minDist = sqrt(dx*dx + dy*dy + dz*dz);
        foreach (ib; 0 .. nBlocks) {
            auto flow = flowBlocks[ib];
            foreach (i; 0 .. flow.ncells) {
                dx = x - flowBlocks[ib]["pos.x", i];
                dy = y - flowBlocks[ib]["pos.y", i];
                dz = (withZpos) ? (z - flowBlocks[ib]["pos.z", i]) : 0.0;
                double dist = sqrt(dx*dx + dy*dy + dz*dz);
                if (dist < minDist) {
                    minDist = dist; nearest[0] = ib; nearest[1] = i;
                }
            } // foreach i
        } // foreach ib
        return nearest;
    } // end find_nearest_cell_centre()

    void subtract_ref_soln(string fileName)
    {
        lua_State* L = luaL_newstate();
        luaL_openlibs(L);
        lua_pushglobaltable(L);
        registerVector3(L);
        registerGlobalConfig(L);
        registerFlowSolution(L);
        registerFlowState(L);
        registerPaths(L);
        registerGpathUtils(L);
        registerSurfaces(L);
        registerVolumes(L);
        registerUnivariateFunctions(L);
        registerStructuredGrid(L);
        registerUnstructuredGrid(L);
        registerGasModel(L);
        registeridealgasflowFunctions(L);
        registergasflowFunctions(L);
        if ( luaL_dofile(L, toStringz(fileName)) != 0 ) {
            writeln("Problem in the user-supplied Lua script: ", fileName);
            string errMsg = to!string(lua_tostring(L, -1));
            throw new FlowSolverException(errMsg);
        }
        foreach (ib; 0 .. nBlocks) {
            flowBlocks[ib].subtract_ref_soln(L);
        }
    } // end subtract_ref_soln()

    void subtract(FlowSolution other)
    {
        // We assume that the flow solutions align geometrically.
        foreach (ib; 0 .. nBlocks) {
            flowBlocks[ib].subtract(other.flowBlocks[ib]);
        }
    }

    double[] compute_volume_weighted_norms(string varName, string regionStr)
    {
        double x0, y0, z0, x1, y1, z1;
        bool limitRegion = false;
        regionStr = regionStr.strip();
        regionStr = regionStr.replaceAll(regex("\""), "");
        if (regionStr.length > 0) {
            auto items = regionStr.split(",");
            x0 = to!double(items[0]); y0 = to!double(items[1]); z0 = to!double(items[2]);
            x1 = to!double(items[3]); y1 = to!double(items[4]); z1 = to!double(items[5]);
            limitRegion = true;
        }
        double L1 = 0.0;
        double L2 = 0.0;
        double Linf = 0.0;
        double[] peak_pos = [0.0, 0.0, 0.0];
        double volume_sum = 0.0;
        foreach (ib; 0 .. nBlocks) {
            auto flow = flowBlocks[ib];
            foreach (i; 0 .. flow.ncells) {
                double x = flowBlocks[ib]["pos.x", i];
                double y = flowBlocks[ib]["pos.y", i];
                double z = (GlobalConfig.dimensions == 3) ? flowBlocks[ib]["pos.z", i] : 0.0;
                if (GlobalConfig.dimensions == 3) {
                    if (limitRegion &&
                        (x < x0 || y < y0 || z < z0 ||
                        x > x1 || y > y1 || z > z1)) continue;
                }
                else {
                    if (limitRegion &&
                        (x < x0 || y < y0 ||
                        x > x1 || y > y1)) continue;
                }
                double volume = flowBlocks[ib]["vol", i];
                double value = flowBlocks[ib][varName, i];
                volume_sum += volume;
                L1 += volume * abs(value);
                L2 += volume * value * value;
                if (abs(value) > Linf) {
                    Linf = abs(value);
                    peak_pos[0] = x; peak_pos[1] = y; peak_pos[2] = z;
                }
            } // foreach i
        } // foreach ib
        L1 /= volume_sum;
        L2 = sqrt(L2/volume_sum);
        return [L1, L2, Linf, peak_pos[0], peak_pos[1], peak_pos[2]];
    } // end compute_volume_weighted_norms()

    string get_value_str(size_t ib, size_t i, string varName)
    {
        string value;
        if (canFind(flowBlocks[ib].variableNames, varName)) {
            value ~= format("%g", flowBlocks[ib][varName, i]);
        } else {
            value ~= "nil";
        }
        return value;
    } // end get_value_str()

    string get_massf_str(size_t ib, size_t i)
    {
        string txt;
        size_t count = 0;
        foreach (varName; flowBlocks[ib].variableNames) {
            if (startsWith(varName, "massf")) {
                if (count > 0) { txt ~= ", "; }
                string shortName = varName.replaceAll(regex("massf"), "");
                txt ~= format("%s=%g", shortName, flowBlocks[ib][varName, i]);
                count += 1;
            }
        }
        return txt;
    } // end get_massf_str()

    void write_vtk_files(string plotDir, string plotName, bool binary_format=false)
    {
        ensure_directory_is_present(plotDir);
        writefln("Writing VTK-XML files to directory \"%s\" with name \"%s\"", plotDir, plotName);
        File visitFile = File(plotDir~"/"~plotName~".visit", "w");
        // For each time index, the visit justs lists the names of the files for individual blocks.
        visitFile.writef("!NBLOCKS %d\n", nBlocks);
        string pvtuFileName = plotName~".pvtu";
        File pvtuFile = begin_PVTU_file(plotDir~"/"~pvtuFileName, flowBlocks[0].variableNames);
        foreach (jb; 0 .. nBlocks) {
            string vtuFileName = plotName~format("-b%04d.vtu", jb);
            add_piece_to_PVTU_file(pvtuFile, vtuFileName);
            visitFile.writef("%s\n", vtuFileName);
            write_VTU_file(flowBlocks[jb], gridBlocks[jb], plotDir~"/"~vtuFileName, binary_format);
        }
        finish_PVTU_file(pvtuFile);
        visitFile.close();
    } // end write_vtk_files()

private:
    bool kdtree_made = false;
    Node[] nodes;
    Node* root;

   /**
    * Assemble the cells in this flow solution into a kdtree for fast nearest neighbour matching.
    *
    *   @author: Nick Gibbons
    */
    void construct_kdtree() {
        // Avoid use of ~= by preallocating the nodes array, since it could get quite large.
        size_t totalcells = 0;
        foreach (ib; 0 .. nBlocks) totalcells += gridBlocks[ib].ncells;
        nodes.length = totalcells;

        size_t j = 0;
        foreach (ib; 0 .. nBlocks) {
            auto blk = gridBlocks[ib];
            foreach (i; 0 .. blk.ncells) {
                Vector3 p = blk.cell_barycentre(i);
                nodes[j].x[0] = p.x.re;
                nodes[j].x[1] = p.y.re;
                nodes[j].x[2] = p.z.re;
                nodes[j].blkid = ib;
                nodes[j].cellid = i;
                j += 1;
            }
        }
        root = makeTree(nodes);
        kdtree_made = true;
    }

} // end class FlowSolution

class FluidBlockLite {
    // Much like the Python library for postprocessing in Eilmer3,
    // we are going to handle the data as a big chunk of numbers,
    // with the label for each variable coming the top of the file.
    //
    // Note that this class is like the FluidBlock class but does not
    // have all of the data space needed for a simulation.
    // The intention is that it is "lighter weight" and so allow
    // postprocessing of workstations that may be considerably smaller
    // than the computer used for generating the simulation data.
public:
    size_t dimensions;
    Grid_t gridType;
    size_t ncells;
    size_t nic;
    size_t njc;
    size_t nkc;
    string[] bcGroups;
    string[] variableNames;
    size_t[string] variableIndex;
    double simTime;
    string flow_format;
    double omegaz; // Angular velocity (in rad/s) of the rotating frame.
                   // There is only one component, about the z-axis.

    size_t single_index(size_t i, size_t j, size_t k=0) const
    in {
        assert (gridType == Grid_t.structured_grid, "invalid index operation for grid");
        assert (i < nic, text("index i=", i, " is invalid, nic=", nic));
        assert (j < njc, text("index j=", j, " is invalid, njc=", njc));
        assert (k < nkc, text("index k=", k, " is invalid, nkc=", nkc));
    }
    do {
        return i + nic*(j + njc*k);
    }

    this(size_t blkID, JSONValue jsonData, Grid_t gridType, string flow_format) {
        this.gridType = gridType;
        this.flow_format = flow_format;
        JSONValue jsonDataForBlk = jsonData["block_" ~ to!string(blkID)];
        // Fill boundary group list if unstructured
        if (gridType == Grid_t.unstructured_grid) {
            size_t nboundaries = getJSONint(jsonDataForBlk, "nboundaries", 0);
            for (size_t i=0; i < nboundaries; i++) {
                auto myGroup = getJSONstring(jsonDataForBlk["boundary_" ~ to!string(i)], "group", "");
                bcGroups ~= myGroup;
            }
        }
        // rotating-frame angular velocity
        omegaz = getJSONdouble(jsonDataForBlk, "omegaz", 0.0);
    } // end "partial" constructor

    this(string filename, size_t blkID, double simTime, JSONValue jsonData, Grid_t gridType, string flow_format, string[] variables, int ncells)
    {
        this.simTime = simTime;
        variableNames = variables.dup;
        foreach (i, var; variables) variableIndex[var] = i;
        this.ncells = ncells;
        this(blkID, jsonData, gridType, flow_format);
        this.readFluidVariablesFromFile(filename, variables, ncells);
    } // end constructor from file

    this(ref const(FluidBlockLite) other, size_t[] cellList, size_t new_dimensions,
         size_t new_nic, size_t new_njc, size_t new_nkc)
    // Construct by extracting a subset of cells from another FluidBlockLite object.
    // Note that the values for the nic, njc and nkc counts should be consistent with
    // cellList.length == new_nic * new_njc * new_nkc
    // For an unstructured grid, new_njc=1 and new_nkc=1 should be provided.
    {
        flow_format = other.flow_format;
        gridType = other.gridType;
        dimensions = new_dimensions;
        omegaz = other.omegaz;
        ncells = cellList.length;
        nic = new_nic; njc = new_njc; nkc = new_nkc;
        simTime = other.simTime;
        variableNames = other.variableNames.dup();
        foreach(i, var; variableNames) { variableIndex[var] = i; }
        _data.length = ncells;
        foreach (i, iother; cellList) {
            _data[i].length = variableNames.length;
            foreach (ivar; 0 .. variableNames.length) {
                _data[i][ivar] = other._data[iother][ivar];
            }
        } // foreach i
        return;
    } // end constructor of subset

    ref double opIndex(string varName, size_t i)
    {
        return _data[i][variableIndex[varName]];
    }

    ref double opIndex(string varName, size_t i, size_t j, size_t k=0)
    {
        return _data[single_index(i,j,k)][variableIndex[varName]];
    }

    string variable_names_as_string(bool with_column_pos=false,
                                    bool with_quote_chars=false,
                                    bool as_comment=false)
    {
        auto writer = appender!string();
        if (as_comment) { formattedWrite(writer, "#"); }
        foreach (i, name; variableNames) {
            string txt = name;
            // Gnuplot column numbers start at 1.
            if (with_column_pos) { txt = format("%d:%s", i+1, txt); }
            if (with_quote_chars) { txt = format("\"%s\"", txt); }
            if ((i==0 && as_comment) || (i>0)) { formattedWrite(writer, " "); }
            formattedWrite(writer, "%s", txt);
        }
        return writer.data;
    }

    string values_as_string(size_t i)
    {
        auto writer = appender!string();
        formattedWrite(writer, "%.18e", _data[i][0]);
        foreach (ivar; 1 .. variableNames.length) {
            formattedWrite(writer, " %.18e", _data[i][ivar]);
        }
        return writer.data;
    }

    string values_as_string(size_t i, size_t j, size_t k)
    {
        return values_as_string(single_index(i,j,k));
    }

    void add_aux_variables(string[] addVarsList, string tag="field")
    // Adds variables to the data for each cell.
    {
        if (tag != "field") return;
        GasModel gmodel = GlobalConfig.gmodel_master;
        GasState Q = GasState(gmodel);
        // Gather massf species names in a list for later use as keys.
        string[] massf_names;
        foreach (isp; 0 .. gmodel.n_species) {
            auto name = cast(char[]) gmodel.species_name(to!int(isp));
            name = tr(name, " \t", "--", "s");
            massf_names ~= "massf-" ~ to!string(name);
        }
        number[] molef; // mole-fractions may be needed
        number[] conc; // concentrations or number-densities may be needed
        //
        // We assume a lot about the data that has been read in so,
        // we need to skip this function if all is not in place
        bool ok_to_proceed = true;
        foreach (name; ["pos.x", "pos.y", "a", "rho", "p", "vel.x", "vel.y", "e"]) {
            if (!canFind(variableNames, name)) { ok_to_proceed = false; }
        }
        if (!ok_to_proceed) {
            writeln("FluidBlockLite.add_aux_variables(): Some essential variables not found.");
            return;
        }
        bool threeD = canFind(variableNames, "pos.z");
        //
        bool add_nrf_velocities = canFind(addVarsList, "nrf"); // nonrotating-frame velocities
        bool add_cyl_coords = canFind(addVarsList, "cyl"); // cylindrical coordinates, r and theta
        bool add_mach = canFind(addVarsList, "mach");
        bool add_pitot_p = canFind(addVarsList, "pitot");
        bool add_total_p = canFind(addVarsList, "total-p");
        bool add_total_h = canFind(addVarsList, "total-h");
        bool add_total_T = canFind(addVarsList, "total-T");
        bool add_enthalpy = canFind(addVarsList, "enthalpy");
        bool add_entropy = canFind(addVarsList, "entropy");
        bool add_molef = canFind(addVarsList, "molef");
        bool add_conc = canFind(addVarsList, "conc"); // concentrations
        VibSpecificNitrogen gmodel2 = cast(VibSpecificNitrogen) gmodel;
        bool add_Tvib = gmodel2 && canFind(addVarsList, "Tvib");
        //
        if (add_cyl_coords) {
            variableNames ~= "r";
            variableIndex["r"] = variableNames.length - 1;
            variableNames ~= "theta";
            variableIndex["theta"] = variableNames.length - 1;
            variableNames ~= "vel.r";
            variableIndex["vel.r"] = variableNames.length - 1;
            variableNames ~= "vel.theta";
            variableIndex["vel.theta"] = variableNames.length - 1;
        }
        if (add_nrf_velocities) {
            variableNames ~= "nrfv.x";
            variableIndex["nrfv.x"] = variableNames.length - 1;
            variableNames ~= "nrfv.y";
            variableIndex["nrfv.y"] = variableNames.length - 1;
            // Note that the z-component is the same for both frames.
        }
        if (add_mach) {
            variableNames ~= "M_local";
            variableIndex["M_local"] = variableNames.length - 1;
        }
        if (add_pitot_p) {
            variableNames ~= "pitot_p";
            variableIndex["pitot_p"] = variableNames.length - 1;
        }
        if (add_total_p) {
            variableNames ~= "total_p";
            variableIndex["total_p"] = variableNames.length - 1;
        }
        if (add_total_h) {
            variableNames ~= "total_h";
            variableIndex["total_h"] = variableNames.length - 1;
        }
        if (add_total_T) {
            variableNames ~= "total_T";
            variableIndex["total_T"] = variableNames.length - 1;
        }
        if (add_enthalpy) {
            variableNames ~= "enthalpy";
            variableIndex["enthalpy"] = variableNames.length - 1;
        }
        if (add_entropy) {
            variableNames ~= "entropy";
            variableIndex["entropy"] = variableNames.length - 1;
        }
        if (add_molef) {
            molef.length = gmodel.n_species;
            foreach (name; massf_names) {
                string new_name = name.replaceAll(regex("massf"), "molef");
                variableNames ~= new_name;
                variableIndex[new_name] = variableNames.length - 1;
            }
        }
        if (add_conc) {
            conc.length = gmodel.n_species;
            foreach (name; massf_names) {
                string new_name = name.replaceAll(regex("massf"), "conc");
                variableNames ~= new_name;
                variableIndex[new_name] = variableNames.length - 1;
            }
        }
        if (add_Tvib) {
            variableNames ~= "Tvib";
            variableIndex["Tvib"] = variableNames.length - 1;
        }
        //
        // Be careful to add auxiliary variable values in the code below
        // in the same order as the list of variable names in the code above.
        //
        foreach (i; 0 .. ncells) {
            double x = _data[i][variableIndex["pos.x"]];
            double y = _data[i][variableIndex["pos.y"]];
            double z = (threeD) ? _data[i][variableIndex["pos.z"]] : 0.0;
            double a = _data[i][variableIndex["a"]];
            double p = _data[i][variableIndex["p"]];
            double rho = _data[i][variableIndex["rho"]];
            double T = _data[i][variableIndex["T"]];
            double g = a*a*rho/p; // approximation for gamma
            // Velocity in the block frame of reference that
            // may be rotating for turbomachinery calculations.
            double wx = _data[i][variableIndex["vel.x"]];
            double wy = _data[i][variableIndex["vel.y"]];
            double wz = (threeD) ? _data[i][variableIndex["vel.z"]] : 0.0;
            double w = sqrt(wx*wx + wy*wy + wz*wz);
            double M = w/a;
            if (add_cyl_coords) {
                // Position and velocity components in cylindrical coordinates.
                double r = hypot(x, y);
                double theta = atan2(y, x);
                double sn = y/r; // sin(theta)
                double cs = x/r; // cos(theta)
                double vel_r = cs*wx + sn*wy;
                double vel_theta = -sn*wx + cs*wy;
                _data[i] ~= r;
                _data[i] ~= theta;
                _data[i] ~= vel_r;
                _data[i] ~= vel_theta;
            }
            if (add_nrf_velocities) {
                // Nonrotating-frame velocities.
                Vector3 pos = Vector3(x, y, z);
                Vector3 vel = Vector3(wx, wy, wz);
                into_nonrotating_frame(vel, pos, omegaz);
                _data[i] ~= vel.x.re;
                _data[i] ~= vel.y.re;
                // Note that z-component is same for both frames.
            }
            if (add_mach) { _data[i] ~= M; }
            if (add_pitot_p) {
                // Rayleigh Pitot formula
                double pitot_p;
                if (M > 1.0) {
                    // Go through the shock and isentropic compression.
                    double t1 = (g+1)*M*M/2;
                    double t2 = (g+1)/(2*g*M*M - (g-1));
                    pitot_p = p * pow(t1,(g/(g-1))) * pow(t2,(1/(g-1)));
                } else {
                    // Isentropic compression only.
                    double t1 = 1 + 0.5*(g-1)*M*M;
                    pitot_p = p * pow(t1,(g/(g-1)));
                }
                _data[i] ~= pitot_p;
            }
            if (add_total_p) {
                // Isentropic process only.
                double t1 = 1 + 0.5*(g-1)*M*M;
                double total_p = p * pow(t1,(g/(g-1)));
                _data[i] ~= total_p;
            }
            if (add_total_h) {
                double e0 = _data[i][variableIndex["e"]];
                double tke;
		double e_int = 0.0;
                if (canFind(variableNames, "tke")) {
                    tke = _data[i][variableIndex["tke"]];
                } else {
                    tke = 0.0;
                }
		// We need to be greedy: search for as many e_modes as we can find.
		// And tally into e_int.
		int imode = 0;
		string e_mode_str = "e_modes[0]";
		while (canFind(variableNames, e_mode_str)) {
		    e_int += _data[i][variableIndex[e_mode_str]];
		    imode++;
		    e_mode_str = format("e_modes[%d]", imode);
		}
		// Sum up the bits of energy.
                double total_h = p/rho + e0 + e_int + 0.5*w*w + tke;
                _data[i] ~= total_h;
            }
            if (add_total_T) {
                double total_T = T * (1.0 + 0.5 * (g - 1.0) * M*M);
                _data[i] ~= total_T;
            }
            if (add_enthalpy) {
                foreach (isp; 0 .. Q.massf.length) {
                    Q.massf[isp] = _data[i][variableIndex[massf_names[isp]]];
                }
                Q.p = p; Q.T = T; Q.rho = rho; Q.u = _data[i][variableIndex["u"]];
                double enthalpy = gmodel.enthalpy(Q).re;
                _data[i] ~= enthalpy;
            }
            if (add_entropy) {
                foreach (isp; 0 .. Q.massf.length) {
                    Q.massf[isp] = _data[i][variableIndex[massf_names[isp]]];
                }
                Q.p = p; Q.T = T;
                double entropy = gmodel.entropy(Q).re;
                _data[i] ~= entropy;
            }
            if (add_molef) {
                foreach (isp; 0 .. Q.massf.length) {
                    Q.massf[isp] = _data[i][variableIndex[massf_names[isp]]];
                }
                Q.p = p; Q.rho = rho; Q.T = T;
                gmodel.massf2molef(Q, molef);
                foreach (mf; molef) { _data[i] ~= mf.re; }
            }
            if (add_conc) {
                foreach (isp; 0 .. Q.massf.length) {
                    Q.massf[isp] = _data[i][variableIndex[massf_names[isp]]];
                }
                Q.p = p; Q.rho = rho; Q.T = T;
                gmodel.massf2conc(Q, conc);
                foreach (c; conc) { _data[i] ~= c.re; }
            }
            if (add_Tvib) {
                foreach (isp; 0 .. Q.massf.length) {
                    Q.massf[isp] = _data[i][variableIndex[massf_names[isp]]];
                }
                Q.p = p; Q.rho = rho; Q.T = T;
                gmodel2.update_thermo_from_pT(Q);
                number Tvib = gmodel2.compute_Tvib(Q, to!number(T));
                _data[i] ~= Tvib.re;
            }
        } // foreach i
    } // end add_aux_variables()


    void subtract_ref_soln(lua_State* L)
    {
        string luaFnName = "refSoln";
        // Test if the user has supplied a reference solution for the flow domain.
        lua_getglobal(L, luaFnName.toStringz);
        if ( lua_isnil(L, -1) ) {
            // Do nothing. Just return.
            lua_pop(L, 1);
            return;
        }
        lua_pop(L, 1);
        foreach (i; 0 .. ncells) {
            // Call back to the Lua function to get a table of values.
            // function refSoln(x, y, z)
            lua_getglobal(L, luaFnName.toStringz);
            lua_pushnumber(L, simTime);
            lua_pushnumber(L, _data[i][variableIndex["pos.x"]]);
            lua_pushnumber(L, _data[i][variableIndex["pos.y"]]);
            if (GlobalConfig.dimensions == 3) {
                lua_pushnumber(L, _data[i][variableIndex["pos.z"]]);
            }
            else {
                lua_pushnumber(L, 0.0);
            }

            if ( lua_pcall(L, 4, 1, 0) != 0 ) {
                string errMsg = "Error in call to " ~ luaFnName ~
                    " from FlowSolution.subtract_ref_soln(): " ~
                    to!string(lua_tostring(L, -1));
                luaL_error(L, errMsg.toStringz);
            }
            // We are expecting a table, containing labelled values.
            if ( !lua_istable(L, -1) ) {
                string errMsg = "Error in FlowSolution.subtract_ref_soln().;\n" ~
                    "A table containing values is expected, but no table was found.";
                luaL_error(L, errMsg.toStringz);
                return;
            }
            // Subtract the ones that are common to the table and the cell.
            foreach (ivar; 0 .. variableNames.length) {
                lua_getfield(L, -1, variableNames[ivar].toStringz);
                double value = 0.0;
                if ( lua_isnumber(L, -1) ) value = lua_tonumber(L, -1);
                lua_pop(L, 1);
                _data[i][ivar] -= value;
            }
            lua_settop(L, 0); // clear the stack
        } // foreach i
    } // end subtract_ref_soln()

    void subtract(FluidBlockLite other)
    {
        foreach (ivar; 0 .. variableNames.length) {
            if (other.variableNames.length > ivar && (other.variableNames[ivar] == variableNames[ivar])) {
                // We assume that the flow solutions align geometrically.
                foreach (i; 0 .. ncells) { _data[i][ivar] -= other._data[i][ivar]; }
            } else {
                writeln("Warning: variable-name mismatch at ivar=%d, name=%s", ivar, variableNames[ivar]);
            }
        }
    } // end subtract()

    // Now that we have moved the function that reads the detailed data to another
    // module, fluidblockio_old, we need to declare the following storage as public.
public:
    double[][] _data;
} // end class FluidBlockLite
