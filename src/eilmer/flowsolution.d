/** flowsolution.d
 * Eilmer4 compressible-flow simulation code, postprocessing functions.
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
import util.lua;
import nm.complex;
import nm.number;
import geom;
import gas;
import gas.vib_specific_nitrogen;
import fvcore;
import globalconfig;
import json_helper;
import core.stdc.stdlib : exit;


class FlowSolution {
    // The collection of flow blocks and grid blocks that define the flow
    // and the domain at one particular instant in time.
public:
    string jobName;
    double sim_time;
    size_t nBlocks;
    BlockFlow[] flowBlocks;
    Grid[] gridBlocks;
    
    this(string jobName, string dir, int tindx, size_t nBlocks)
    {
        // -- initialising JSONData
        // We need to set a few things here if we are constructing this object
        // in a custom-postprocessing context.
        string configFileName = dir ~ "/config/" ~ jobName ~ ".config";
        string content;
        try {
            content = readText(configFileName);
        } catch (Exception e) {
            writeln("Failed to read config file: ", configFileName);
            writeln("Message is: ", e.msg);
            exit(1);
        }
        JSONValue jsonData;
        try {
            jsonData = parseJSON!string(content);
        } catch (Exception e) {
            writeln("Failed to parse JSON from config file: ", configFileName);
            writeln("Message is: ", e.msg);
            exit(1);
        }
        GlobalConfig.grid_format = jsonData["grid_format"].str;
        GlobalConfig.flow_format = jsonData["flow_format"].str;
        string gridFileExt = "gz"; if (GlobalConfig.grid_format == "rawbinary") { gridFileExt = "bin"; }
        string flowFileExt = "gz"; if (GlobalConfig.flow_format == "rawbinary") { flowFileExt = "bin"; }
        GlobalConfig.grid_motion = grid_motion_from_name(jsonData["grid_motion"].str);
        // -- end initialising JSONData
        //
        // Use job.list to get a hint of the type of each block.
        auto listFile = File(dir ~ "/config/" ~ jobName ~ ".list");
        auto listFileLine = listFile.readln().chomp(); // only comments on the first line
        //
        foreach (ib; 0 .. nBlocks) {
            listFileLine = listFile.readln().chomp();
            int ib_chk;
            string gridTypeName;
            string label;
            formattedRead(listFileLine, " %d %s %s", &ib_chk, &gridTypeName, &label);
            if (ib != ib_chk) {
                string msg = format("Reading .list file ib=%d ib_chk=%d", ib, ib_chk);
                throw new FlowSolverException(msg);
            }
            auto gridType = gridTypeFromName(gridTypeName);
            string fileName;
            if (GlobalConfig.grid_motion != GridMotion.none) {
                fileName = make_file_name!"grid"(jobName, to!int(ib), tindx, gridFileExt);
            } else {
                fileName = make_file_name!"grid"(jobName, to!int(ib), 0, gridFileExt);
            }
            fileName = dir ~ "/" ~ fileName;
            final switch (gridType) {
            case Grid_t.structured_grid:
                gridBlocks ~= new StructuredGrid(fileName, GlobalConfig.grid_format);
                break;
            case Grid_t.unstructured_grid:
                gridBlocks ~= new UnstructuredGrid(fileName, GlobalConfig.grid_format);
            }
            gridBlocks[$-1].sort_cells_into_bins();
            fileName = make_file_name!"flow"(jobName, to!int(ib), tindx, flowFileExt);
            fileName = dir ~ "/" ~ fileName;
            flowBlocks ~= new BlockFlow(fileName, ib, jsonData, gridType);
        } // end foreach ib
        this.jobName = jobName;
        this.nBlocks = nBlocks;
        sim_time = flowBlocks[0].sim_time;
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
        str ~= "jobName=" ~ jobName ~ ", sim_time=" ~ to!string(sim_time);
        str ~= ", nBlocks=" ~ to!string(nBlocks);
        str ~= ")";
        return str;
    }

    void add_aux_variables(string[] addVarsList)
    {
        foreach (blk; flowBlocks) blk.add_aux_variables(addVarsList);
    }

    size_t[] find_enclosing_cell(ref const(Vector3) p)
    {
        size_t[] cell_identity = [0, 0, 0]; // blk_id, i, found_flag
        foreach (ib; 0 .. nBlocks) {
            bool found = false;
            size_t indx = 0;
            gridBlocks[ib].find_enclosing_cell(p, indx, found);
            if (found) {
                cell_identity = [ib, indx, (found)?1:0];
                break;
            }
        } // foreach ib
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
        size_t[] nearest = [0, 0]; // blk_id, i
        double dx = x - flowBlocks[0]["pos.x", 0];
        double dy = y - flowBlocks[0]["pos.y", 0];
        double dz = z - flowBlocks[0]["pos.z", 0];
        double minDist = sqrt(dx*dx + dy*dy + dz*dz);
        foreach (ib; 0 .. nBlocks) {
            auto flow = flowBlocks[ib];
            foreach (i; 0 .. flow.ncells) {
                dx = x - flowBlocks[ib]["pos.x", i];
                dy = y - flowBlocks[ib]["pos.y", i];
                dz = z - flowBlocks[ib]["pos.z", i];
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
        if ( luaL_dofile(L, toStringz(fileName)) != 0 ) {
            writeln("Problem in the user-supplied Lua script: ", fileName);
            string errMsg = to!string(lua_tostring(L, -1));
            throw new FlowSolverException(errMsg);
        }
        foreach (ib; 0 .. nBlocks) {
            flowBlocks[ib].subtract_ref_soln(L);
        }
    } // end subtract_ref_soln()

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
                double z = flowBlocks[ib]["pos.z", i];
                if (limitRegion && 
                    (x < x0 || y < y0 || z < z0 ||
                     x > x1 || y > y1 || z > z1)) continue;
                double volume = flowBlocks[ib]["volume", i];
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
} // end class FlowSolution

class BlockFlow {
    // Much like the Python library for postprocessing in Eilmer3,
    // we are going to handle the data as a big chunk of numbers,
    // with the label for each variable coming the top of the file.
    //
    // Note that this class is like the SFluidBlock class but does not 
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
    double sim_time;

    size_t single_index(size_t i, size_t j, size_t k=0) const
    in {
        assert (gridType == Grid_t.structured_grid, "invalid index operation for grid");
        assert (i < nic, text("index i=", i, " is invalid, nic=", nic));
        assert (j < njc, text("index j=", j, " is invalid, njc=", njc));
        assert (k < nkc, text("index k=", k, " is invalid, nkc=", nkc));
    }
    body {
        return i + nic*(j + njc*k);
    }

    this(string filename, size_t blkID, JSONValue jsonData, Grid_t gridType)
    {
        this.gridType = gridType;
        string myLabel;
        size_t nvariables;
        // Read in the flow data for a single block.
        //
        // Keep in sync with: 
        // 1. SFluidBlock.write_solution(),
        // 2. UFluidBlock.write_solution()
        // 3. write_initial_sg_flow_file_from_lua() in luaflowstate.d
        // 4. write_initial_usg_flow_file_from_lua() in luaflowstate.d.
        //
        switch (GlobalConfig.flow_format) {
        case "gziptext": goto default;
        case "rawbinary":
            File fin = File(filename, "rb");
            string expected_header;
            final switch (gridType) {
            case Grid_t.structured_grid: expected_header = "structured_grid_flow 1.0"; break;
            case Grid_t.unstructured_grid: expected_header = "unstructured_grid_flow 1.0";
            }
            char[] found_header = new char[expected_header.length];
            fin.rawRead(found_header);
            if (found_header != expected_header) {
                throw new FlowSolverException("BlockFlow constructor, read_solution from raw_binary_file: "
                                              ~ "unexpected header: " ~ to!string(found_header)); 
            }
            int[1] int1; fin.rawRead(int1);
            int label_length = int1[0];
            if (label_length > 0) {
                char[] found_label = new char[label_length];
                fin.rawRead(found_label);
                myLabel = to!string(found_label);
            }
            double[1] dbl1; fin.rawRead(dbl1); sim_time = dbl1[0];
            fin.rawRead(int1); nvariables = int1[0];
            foreach(i; 0 .. nvariables) {
                char[] varname; fin.rawRead(int1); varname.length = int1[0]; 
                fin.rawRead(varname); variableNames ~= to!string(varname);
            }
            fin.rawRead(int1); int my_dimensions = int1[0];
            final switch (gridType) {
            case Grid_t.structured_grid:
                int[3] int3; fin.rawRead(int3);
                nic = int3[0]; njc = int3[1]; nkc = int3[2];
                ncells = nic*njc*nkc;
                break;
            case Grid_t.unstructured_grid:
                fin.rawRead(int1); ncells = int1[0];
                nic = ncells; njc = 1; nkc = 1;
            }
            // Assume the remaining data is all double type and is in standard cell order.
            _data.length = ncells;
            foreach (i; 0 .. ncells) {
                _data[i].length = nvariables;
                fin.rawRead(_data[i]);
            }
            fin.close();
            break;
        default:
            string[] tokens;
            auto byLine = new GzipByLine(filename);
            auto line = byLine.front; byLine.popFront();
            string format_version;
            final switch (gridType) {
            case Grid_t.structured_grid:
                formattedRead(line, "structured_grid_flow %s", &format_version);
                break;
            case Grid_t.unstructured_grid:
                formattedRead(line, "unstructured_grid_flow %s", &format_version);
            }
            if (format_version != "1.0") {
                string msg = text("BlockFlow.read_solution(): " ~
                                  "format version found: " ~ format_version);
                throw new FlowSolverException(msg); 
            }
            line = byLine.front; byLine.popFront();
            formattedRead(line, "label: %s", &myLabel);
            line = byLine.front; byLine.popFront();
            formattedRead(line, "sim_time: %g", &sim_time);
            line = byLine.front; byLine.popFront();
            formattedRead(line, "variables: %d", &nvariables);
            line = byLine.front; byLine.popFront();
            variableNames = line.strip().split();
            foreach (ref var; variableNames) { var = var.replaceAll(regex("\""), ""); }
            line = byLine.front; byLine.popFront();
            formattedRead(line, "dimensions: %d", &dimensions);
            final switch (gridType) {
            case Grid_t.structured_grid:
                line = byLine.front; byLine.popFront();
                formattedRead(line, "nicell: %d", &nic);
                line = byLine.front; byLine.popFront();
                formattedRead(line, "njcell: %d", &njc);
                line = byLine.front; byLine.popFront();
                formattedRead(line, "nkcell: %d", &nkc);
                ncells = nic*njc*nkc;
                break;
            case Grid_t.unstructured_grid:
                line = byLine.front; byLine.popFront();
                formattedRead(line, "ncells: %d", &ncells);
                nic = ncells; njc = 1; nkc = 1;
            }
            // Scan the remainder of the file, extracting our data.
            // Assume it is in standard cell order.
            _data.length = ncells;
            foreach (i; 0 .. ncells) {
                line = byLine.front; byLine.popFront();
                tokens = line.strip().split();
                assert(tokens.length == variableNames.length,
                       "wrong number of items for variable data");
                _data[i].length = variableNames.length;
                foreach (ivar; 0 .. variableNames.length) {
                    _data[i][ivar] = to!double(tokens[ivar]);
                }
            } // foreach i
            byLine.range.f.close();
        } // end switch flow_format
        foreach (i; 0 .. variableNames.length) { variableIndex[variableNames[i]] = i; }
        // Fill boundary group list
        size_t nboundaries = getJSONint(jsonData["block_" ~ to!string(blkID)], "nboundaries", 0);
        for (size_t i=0; i < nboundaries; i++) {
            auto myGroup = getJSONstring(jsonData["block_" ~ to!string(blkID)]["boundary_" ~ to!string(i)], "group", "");
            bcGroups ~= myGroup;
        }
    } // end constructor from file

    this(ref const(BlockFlow) other, size_t[] cellList, size_t new_dimensions,
         size_t new_nic, size_t new_njc, size_t new_nkc)
    // Construct by extracting a subset of cells from another BlockFlow object.
    // Note that the values for the nic, njc and nkc counts should be consistent with
    // cellList.length == new_nic * new_njc * new_nkc
    // For an unstructured grid, new_njc=1 and new_nkc=1 should be provided.
    {
        gridType = other.gridType;
        dimensions = new_dimensions;
        ncells = cellList.length;
        nic = new_nic; njc = new_njc; nkc = new_nkc;
        sim_time = other.sim_time;
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

    string variable_names_as_string(bool with_column_pos=false)
    {
        auto writer = appender!string();
        formattedWrite(writer, "#");
        foreach (i, name; variableNames) {
            if (with_column_pos) {
                formattedWrite(writer, " %d:%s", i+1, name);
            }
            else {
                formattedWrite(writer, " \"%s\"", name);
            }
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

    void add_aux_variables(string[] addVarsList)
    // Adds variables to the data for each cell.
    {
        GasModel gmodel = GlobalConfig.gmodel_master;
        GasState Q = new GasState(gmodel);
        // Gather massf species names in a list for later use as keys.
        string[] massf_names;
        foreach (isp; 0 .. gmodel.n_species) {
            auto name = cast(char[]) gmodel.species_name(to!int(isp));
            name = tr(name, " \t", "--", "s");
            massf_names ~= "massf[" ~ to!string(isp) ~ "]-" ~ to!string(name);
        }
        number[] molef; // mole-fractions may be needed
        number[] conc; // concentrations or number-densities may be needed
        //
        // We assume a lot about the data that has been read in so,
        // we need to skip this function if all is not in place
        bool ok_to_proceed = true;
        foreach (name; ["a", "rho", "p", "vel.x", "vel.y", "vel.z", "u", "tke"]) {
            if (!canFind(variableNames, name)) { ok_to_proceed = false; }
        }
        if (!ok_to_proceed) {
            writeln("BlockFlow.add_aux_variables(): Some essential variables not found.");
            return;
        }
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
            double a = _data[i][variableIndex["a"]];
            double p = _data[i][variableIndex["p"]];
            double rho = _data[i][variableIndex["rho"]];
            double T = _data[i][variableIndex["T"]];
            double g = a*a*rho/p; // approximation for gamma
            // Velocity in the block frame of reference that
            // may be rotating for turbomachinery calculations.
            double wx = _data[i][variableIndex["vel.x"]];
            double wy = _data[i][variableIndex["vel.y"]];
            double wz = _data[i][variableIndex["vel.z"]];
            double w = sqrt(wx*wx + wy*wy + wz*wz);
            double M = w/a;
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
                double e0 = _data[i][variableIndex["u"]];
                double tke = _data[i][variableIndex["tke"]];
                // Sum up the bits of energy,
                // forgetting the multiple energy modes, for the moment.
                double total_h = p/rho + e0 + 0.5*w*w + tke;
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
                // [TODO] 2018-09-15: Rowan, should we set Tvib to zero if massf[1] is essentially zero?
                number Tvib = gmodel2.compute_Tvib(Q, to!number(T), to!number(T+100.0), 0.1);
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
            lua_pushnumber(L, sim_time);
            lua_pushnumber(L, _data[i][variableIndex["pos.x"]]);
            lua_pushnumber(L, _data[i][variableIndex["pos.y"]]);
            lua_pushnumber(L, _data[i][variableIndex["pos.z"]]);
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

    void read_extra_vars_from_file(string filename, string[] extraVars)
    {
        string[] tokens;
        auto byLine = new GzipByLine(filename);
        auto line = byLine.front; byLine.popFront();
        string format_version;
        string dummyName;
        final switch (gridType) {
        case Grid_t.structured_grid:
            formattedRead(line, "%s  %s", &dummyName, &format_version);
            break;
        case Grid_t.unstructured_grid:
            formattedRead(line, "%s %s", &dummyName, &format_version);
        }
        if (format_version != "1.0") {
            string msg = text("BlockFlow.read_extra_vars_from_file(): " ~
                              "format version found: " ~ format_version);
            throw new FlowSolverException(msg); 
        }
        // Throw away label line.
        byLine.popFront();
        // Throw away variables line.
        byLine.popFront();
        // Read actual list of variables.
        line = byLine.front; byLine.popFront();
        string[] varNames = line.strip().split();
        if (varNames.length != extraVars.length) {
            string msg = text("BlockFlow.read_extra_vars_from_file(): " ~
                              "number of variables in file does not match user-supplied variables list.");
            throw new FlowSolverException(msg); 
        }
        foreach (var; extraVars) {
            variableIndex[var] = variableNames.length;
            variableNames ~= var;
        }
        line = byLine.front; byLine.popFront();
        formattedRead(line, "dimensions: %d", &dimensions);
        final switch (gridType) {
        case Grid_t.structured_grid:
            line = byLine.front; byLine.popFront();
            formattedRead(line, "nicell: %d", &nic);
            line = byLine.front; byLine.popFront();
            formattedRead(line, "njcell: %d", &njc);
            line = byLine.front; byLine.popFront();
            formattedRead(line, "nkcell: %d", &nkc);
            ncells = nic*njc*nkc;
            break;
        case Grid_t.unstructured_grid:
            line = byLine.front; byLine.popFront();
            formattedRead(line, "ncells: %d", &ncells);
            nic = ncells; njc = 1; nkc = 1;
        }
        if (_data.length != ncells) {
            string msg = text("BlockFlow.read_extra_vars_from_file(): " ~
                              "number of cells in file does not match number of existing cells in block.");
            throw new FlowSolverException(msg); 
        }
        // Scan the remainder of the file, extracting our data.
        // Assume it is in standard cell order.
        foreach (i; 0 .. ncells) {
            line = byLine.front; byLine.popFront();
            tokens = line.strip().split();
            assert(tokens.length == extraVars.length,
                   "wrong number of items for variable data");
            foreach (ivar; 0 .. extraVars.length) {
                _data[i] ~= to!double(tokens[ivar]);
            }
        } // foreach i
        byLine.range.f.close();
    } // end read_extra_vars_from_file

private:
    double[][] _data;
} // end class BlockFlow
