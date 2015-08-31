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
import std.algorithm;
import std.array;
import std.math;

import util.lua;
import gzip;
import fileutil;
import sgrid;
import globalconfig;

class SolidSolution {
    // This holds the collection of solid blocks and solid grid blocks
    // that define the solid domain at a particular instant in time.
public:
    double sim_time;
    size_t nBlocks;
    SBlockSolid[] solidBlocks;
    StructuredGrid[] gridBlocks;

    this(string jobName, string dir, int tindx, size_t nBlocks)
    {
	foreach (ib; 0 .. nBlocks) {
	    string fileName;
	    // Presently, don't allow moving grid for solid domain,
	    // so tindx is always 0 for solid-grid
	    fileName = make_file_name!"solid-grid"(jobName, to!int(ib), 0);
	    fileName = dir ~ "/" ~ fileName;
	    gridBlocks ~= new StructuredGrid(fileName, GridFileFormat.gziptext);
	    fileName = make_file_name!"solid"(jobName, to!int(ib), tindx);
	    fileName = dir ~ "/" ~ fileName;
	    solidBlocks ~= new SBlockSolid(fileName);
	} // end foreach ib
	this.nBlocks = nBlocks;
	sim_time = solidBlocks[0].sim_time;
    } // end constructor

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
	regionStr = removechars(regionStr, "\"");
	if (regionStr.length > 0) {
	    auto items = regionStr.split(",");
	    x0 = to!double(items[0]); y0 = to!double(items[1]); z0 = to!double(items[2]);
	    x1 = to!double(items[3]); y1 = to!double(items[4]); z1 = to!double(items[5]);
	    limitRegion = true;
	    // writeln(format("    limit region to x0=%g y0=%g z0=%g x1=%g y1=%g z1=%g",
	    //		   x0, y0, z0, x1, y1, z1));
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
			double z = solidBlocks[ib]["pos.z", i, j, k];
			if (limitRegion && 
			    (x < x0 || y < y0 || z < z0 ||
			     x > x1 || y > y1 || z > z1)) continue;
			double volume = solidBlocks[ib]["volume", i, j, k];
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

} // end class SolidSolution

class SBlockSolid {
public:
    size_t nic;
        size_t njc;
    size_t nkc;
    string[] variableNames;
    size_t[string] variableIndex;
    double sim_time;

    this(string filename)
    {
	// Read in the flow data for a single structured block.
	string[] tokens;
	auto byLine = new GzipByLine(filename);
	auto line = byLine.front; byLine.popFront();
	formattedRead(line, " %g", &sim_time);
	line = byLine.front; byLine.popFront();
	variableNames = line.strip().split();
	foreach (ref var; variableNames) { var = removechars(var, "\""); }
	foreach (i; 0 .. variableNames.length) { variableIndex[variableNames[i]] = i; }
	line = byLine.front; byLine.popFront();
	formattedRead(line, "%d %d %d", &nic, &njc, &nkc);
	_data.length = nic;
	// Resize the storage for our block of data.
	foreach (i; 0 .. nic) {
	    _data[i].length = njc;
	    foreach (j; 0 .. njc) {
		_data[i][j].length = nkc;
		foreach (k; 0 .. nkc) {
		    _data[i][j][k].length = variableNames.length;
		} // foreach k
	    } // foreach j
	} // foreach i
	// Scan the remainder of the file, extracting our data.
	foreach (k; 0 .. nkc) {
	    foreach (j; 0 .. njc) {
		foreach (i; 0 .. nic) {
		    line = byLine.front; byLine.popFront();
		    tokens = line.strip().split();
		    assert(tokens.length == variableNames.length, "wrong number of items");
		    foreach (ivar; 0 .. variableNames.length) {
			_data[i][j][k][ivar] = to!double(tokens[ivar]);
		    }
		} // foreach i
	    } // foreach j
	} // foreach k
    } // end constructor from file

    ref double opIndex(string varName, size_t i, size_t j, size_t k=0)
    {
	return _data[i][j][k][variableIndex[varName]];
    }

    string variable_names_as_string()
    {
	auto writer = appender!string();
	formattedWrite(writer, "#");
	foreach(name; variableNames) {
	    formattedWrite(writer, " \"%s\"", name);
	}
	return writer.data;
    }

    string values_as_string(size_t i, size_t j, size_t k)
    {
	auto writer = appender!string();
	formattedWrite(writer, "%e", _data[i][j][k][0]);
	foreach (ivar; 1 .. variableNames.length) {
	    formattedWrite(writer, " %e", _data[i][j][k][ivar]);
	}
	return writer.data;
    }

    void subtract_ref_soln(lua_State* L)
    {
	string luaFnName = "refSolidSoln";
	foreach (k; 0 .. nkc) {
	    foreach (j; 0 .. njc) {
		foreach (i; 0 .. nic) {
		    // Call back to the Lua function to get a table of values.
		    // function refSolidSoln(x, y, z)
		    lua_getglobal(L, luaFnName.toStringz);
		    lua_pushnumber(L, _data[i][j][k][variableIndex["pos.x"]]);
		    lua_pushnumber(L, _data[i][j][k][variableIndex["pos.y"]]);
		    lua_pushnumber(L, _data[i][j][k][variableIndex["pos.z"]]);
		    if ( lua_pcall(L, 3, 1, 0) != 0 ) {
			string errMsg = "Error in call to " ~ luaFnName ~ 
			    " from LuaFnPath:opCall(): " ~ 
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
			_data[i][j][k][ivar] -= value;
		    }
		    lua_settop(L, 0); // clear the stack
		} // foreach i
	    } // foreach j
	} // foreach k
    } // end subtract_ref_soln()

private:
    double[][][][] _data;
} // end class SBlockSolid
