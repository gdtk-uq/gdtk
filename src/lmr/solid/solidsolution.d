/**
 * solidsolution.d
 * Postprocessing functions related to solid domain.
 *
 * This borrows the main idea from the flowsolution module,
 * namely a class to handle some of the solid data without
 * using all of the core code machinery.
 *
 * Author: Rowan G. and Peter J.
 * Date: 2015-08-04
 */

module solidsolution;

import std.stdio;
import std.conv;
import std.format;
import std.string;
import std.regex;
import std.algorithm;
import std.array;
import std.math;
import std.json;
import std.file : readText;

import util.lua;
import geom;
import gzip;
import fileutil;
import globalconfig;
import lmrconfig;
import blockio : readVariablesFromMetadata, readSolidVariablesFromFile;

class SolidSolution {
    // This holds the collection of solid blocks and solid grid blocks
    // that define the solid domain at a particular instant in time.
public:
    double simTime;
    size_t nBlocks;
    SolidBlockLite[] solidBlocks;
    StructuredGrid[] gridBlocks;

    this(int snapshot, size_t nBlocks, double simTime=-1.0, string dir=".")
    {
        size_t nFluidBlocks = GlobalConfig.nFluidBlocks;
        // We'll open the config file ourselves to read
        // rather than invoke the entire readConfig()
        // because we just want two things:
        //   1. field_format; and
        //   2. grid_format
        string cfgFile =  dir ~ "/" ~ lmrCfg.cfgFile;
        string content;
        try {
            content = readText(cfgFile);
        }
        catch (Exception e) {
            writeln("Message is: ", e.msg);
            throw new Error(text("Failed to read config file: ", cfgFile));
        }
        JSONValue jsonData;
        try {
            jsonData = parseJSON!string(content);
        }
        catch (Exception e) {
            writeln("Message is: ", e.msg);
            throw new Error(text("Failed to parse JSON from config file: ", cfgFile));
        }
        string fieldFmt = jsonData["field_format"].str;
        GlobalConfig.field_format = fieldFmt;
        string gridFmt = jsonData["grid_format"].str;
        GlobalConfig.grid_format = gridFmt;

        // Find variables from metadata file
        string solidMetadataFile = dir ~ "/" ~ lmrCfg.solidMetadataFile;
        auto variables = readVariablesFromMetadata(solidMetadataFile);

        foreach (ib; 0 .. nBlocks) {
            auto bid = ib + nFluidBlocks;
            string gName = dir ~ "/" ~ gridFilename(lmrCfg.initialFieldDir, to!int(bid));
            gridBlocks ~= new StructuredGrid(gName, gridFmt);
            auto g = gridBlocks[$-1];
            g.sort_cells_into_bins();
            int ncells = to!int(g.ncells);
            string sName = dir ~ "/" ~ solidFilename(snapshot, to!int(bid));
            auto sb = new SolidBlockLite(sName, simTime, fieldFmt, variables, ncells);
            sb.nic = g.niv - 1;
            sb.njc = g.njv - 1;
            sb.nkc = max(g.nkv - 1, 1);
            solidBlocks ~= sb;
        }
        this.nBlocks = nBlocks;
        this.simTime = simTime;
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
        size_t[] nearest = [0, 0, 0, 0]; // blk_id, i, j, k
        double dx = x - solidBlocks[0]["pos.x", 0, 0, 0];
        double dy = y - solidBlocks[0]["pos.y", 0, 0, 0];
        double dz = z - solidBlocks[0]["pos.z", 0, 0, 0];
        double minDist = sqrt(dx*dx + dy*dy + dz*dz);
        foreach (ib; 0 .. nBlocks) {
            auto solid = solidBlocks[ib];
            foreach (k; 0 .. solid.nkc) {
                foreach (j; 0 .. solid.njc) {
                    foreach (i; 0 .. solid.nic) {
                        dx = x - solidBlocks[ib]["pos.x", i, j, k];
                        dy = y - solidBlocks[ib]["pos.y", i, j, k];
                        dz = z - solidBlocks[ib]["pos.z", i, j, k];
                        double dist = sqrt(dx*dx + dy*dy + dz*dz);
                        if (dist < minDist) {
                            minDist = dist;
                            nearest = [ib, i, j, k];
                        }
                    } // foreach i
                } // foreach j
            } // foreach k
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
            throw new Error(errMsg);
        }
        foreach (ib; 0 .. nBlocks) {
            solidBlocks[ib].subtract_ref_soln(L);
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
            auto solid = solidBlocks[ib];
            foreach (k; 0 .. solid.nkc) {
                foreach (j; 0 .. solid.njc) {
                    foreach (i; 0 .. solid.nic) {
                        double x = solidBlocks[ib]["pos.x", i, j, k];
                        double y = solidBlocks[ib]["pos.y", i, j, k];
                        double z = (GlobalConfig.dimensions == 3) ? solidBlocks[ib]["pos.z", i, j, k] : 0.0;
                        if (limitRegion && 
                            (x < x0 || y < y0 || z < z0 ||
                             x > x1 || y > y1 || z > z1)) continue;
                        double volume = solidBlocks[ib]["vol", i, j, k];
                        double value = solidBlocks[ib][varName, i, j, k];
                        volume_sum += volume;
                        L1 += volume * abs(value);
                        L2 += volume * value * value;
                        if (abs(value) > Linf) {
                            Linf = abs(value);
                            peak_pos[0] = x; peak_pos[1] = y; peak_pos[2] = z;
                        }
                    } // foreach i
                } // foreach j
            } // foreach k
        } // foreach ib
        L1 /= volume_sum;
        L2 = sqrt(L2/volume_sum);
        return [L1, L2, Linf, peak_pos[0], peak_pos[1], peak_pos[2]];
    } // end compute_volume_weighted_norms()

    string get_value_str(size_t ib, size_t i, size_t j, size_t k, string varName)
    {
        string value;
        if (canFind(solidBlocks[ib].variableNames, varName)) {
            value ~= format("%g", solidBlocks[ib][varName, i, j, k]);
        } else {
            value ~= "nil";
        }
        return value;
    } // end get_value_str()


} // end class SolidSolution

class SolidBlockLite {
public:
    double[][] data; // public because we read into it in blockio module.
    double simTime;
    size_t ncells;
    size_t nic;
    size_t njc;
    size_t nkc;
    string[] variableNames;
    size_t[string] variableIndex;
    string fieldFmt;
    double sim_time;

    size_t single_index(size_t i, size_t j, size_t k=0) const
    in {
        assert (i < nic, text("index i=", i, " is invalid, nic=", nic));
        assert (j < njc, text("index j=", j, " is invalid, njc=", njc));
        assert (k < nkc, text("index k=", k, " is invalid, nkc=", nkc));
    }
    do {
        return i + nic*(j + njc*k);
    }

    this(string filename, double simTime, string fieldFmt, string[] variables, int ncells)
    {
        this.simTime = simTime;
        variableNames = variables.dup;
        foreach (i, var; variables) variableIndex[var] = i;
        this.ncells = ncells;
        this.fieldFmt = fieldFmt;
        this.readSolidVariablesFromFile(filename, variables, ncells);
    }

    size_t[] to_ijk_indices(size_t gid) const
    {
        size_t k = gid / (njc * nic);
        size_t j = (gid - k * (njc * nic)) / nic;
        size_t i = gid - k * (njc * nic) - j * nic;
        return [i, j, k];
    }

    ref double opIndex(string varName, size_t i)
    {
        return data[i][variableIndex[varName]];
    }

    ref double opIndex(string varName, size_t i, size_t j, size_t k=0)
    {
        return data[single_index(i,j,k)][variableIndex[varName]];
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
        formattedWrite(writer, "%.18e", data[i][0]);
        foreach (ivar; 1 .. variableNames.length) {
            formattedWrite(writer, " %.18e", data[i][ivar]);
        }
        return writer.data;
    }

    string values_as_string(size_t i, size_t j, size_t k)
    {
        return values_as_string(single_index(i,j,k));
    }

    void subtract_ref_soln(lua_State* L)
    {
        string luaFnName = "refSolidSoln";
        // Test if the user has supplied a reference solution for the flow domain.
        lua_getglobal(L, luaFnName.toStringz);
        if ( lua_isnil(L, -1) ) {
            // Do nothing. Just return.
            lua_pop(L, 1);
            return;
        }
        lua_pop(L, 1);
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. nic) {
                    // Call back to the Lua function to get a table of values.
                    // function refSolidSoln(x, y, z)
                    lua_getglobal(L, luaFnName.toStringz);
                    lua_pushnumber(L, simTime);
                    lua_pushnumber(L, data[single_index(i,j,k)][variableIndex["pos.x"]]);
                    lua_pushnumber(L, data[single_index(i,j,k)][variableIndex["pos.y"]]);
                    if (GlobalConfig.dimensions == 3) {
                        lua_pushnumber(L, data[single_index(i,j,k)][variableIndex["pos.z"]]);
                    }
                    else {
                        lua_pushnumber(L, 0.0);
                    }
                    if ( lua_pcall(L, 4, 1, 0) != 0 ) {
                        string errMsg = "Error in call to " ~ luaFnName ~ 
                            " from user-supplied reference solution file: " ~ 
                            to!string(lua_tostring(L, -1));
                        luaL_error(L, errMsg.toStringz);
                    }
                    // We are expecting a table, containing labelled values.
                    if ( !lua_istable(L, -1) ) {
                        string errMsg = `Error in SolidSolution.subtract_ref_soln().;
A table containing arguments is expected, but no table was found.`;
                        luaL_error(L, errMsg.toStringz);
                        return;
                    }
                    // Subtract the ones that are common to the table and the cell.
                    foreach (ivar; 0 .. variableNames.length) {
                        lua_getfield(L, -1, variableNames[ivar].toStringz);
                        double value = 0.0;
                        if ( lua_isnumber(L, -1) ) value = lua_tonumber(L, -1);
                        lua_pop(L, 1);
                        data[single_index(i,j,k)][ivar] -= value;
                    }
                    lua_settop(L, 0); // clear the stack
                } // foreach i
            } // foreach j
        } // foreach k
    } // end subtract_ref_soln()

} // end class SolidBlockLite
