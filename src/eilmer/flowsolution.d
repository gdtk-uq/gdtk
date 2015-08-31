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
import std.algorithm;
import std.array;
import std.math;
import gzip;
import fileutil;
import util.lua;
import geom;
import sgrid;
import gas;
import globalconfig;


class FlowSolution {
    // The collection of flow blocks and grid blocks that define the flow
    // and the domain at one particular instant in time.
public:
    double sim_time;
    size_t nBlocks;
    SBlockFlow[] flowBlocks;
    StructuredGrid[] gridBlocks;
    
    this(string jobName, string dir, int tindx, size_t nBlocks)
    {
	foreach (ib; 0 .. nBlocks) {
	    string fileName;
	    if (GlobalConfig.moving_grid) {
		fileName = make_file_name!"grid"(jobName, to!int(ib), tindx);
	    } else {
		fileName = make_file_name!"grid"(jobName, to!int(ib), 0);
	    }
	    fileName = dir ~ "/" ~ fileName;
	    gridBlocks ~= new StructuredGrid(fileName, GridFileFormat.gziptext);
	    fileName = make_file_name!"flow"(jobName, to!int(ib), tindx);
	    fileName = dir ~ "/" ~ fileName;
	    flowBlocks ~= new SBlockFlow(fileName);
	} // end foreach ib
	this.nBlocks = nBlocks;
	sim_time = flowBlocks[0].sim_time;
    } // end constructor

    void add_aux_variables(string[] addVarsList)
    {
	foreach (blk; flowBlocks) blk.add_aux_variables(addVarsList);
    }

    size_t[] find_nearest_cell_centre(double x, double y, double z=0.0)
    {
	size_t[] nearest = [0, 0, 0, 0]; // blk_id, i, j, k
	double dx = x - flowBlocks[0]["pos.x", 0, 0, 0];
	double dy = y - flowBlocks[0]["pos.y", 0, 0, 0];
	double dz = z - flowBlocks[0]["pos.z", 0, 0, 0];
	double minDist = sqrt(dx*dx + dy*dy + dz*dz);
	foreach (ib; 0 .. nBlocks) {
	    auto flow = flowBlocks[ib];
	    foreach (k; 0 .. flow.nkc) {
		foreach (j; 0 .. flow.njc) {
		    foreach (i; 0 .. flow.nic) {
			dx = x - flowBlocks[ib]["pos.x", i, j, k];
			dy = y - flowBlocks[ib]["pos.y", i, j, k];
			dz = z - flowBlocks[ib]["pos.z", i, j, k];
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
	    flowBlocks[ib].subtract_ref_soln(L);
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
	    auto flow = flowBlocks[ib];
	    foreach (k; 0 .. flow.nkc) {
		foreach (j; 0 .. flow.njc) {
		    foreach (i; 0 .. flow.nic) {
			double x = flowBlocks[ib]["pos.x", i, j, k];
			double y = flowBlocks[ib]["pos.y", i, j, k];
			double z = flowBlocks[ib]["pos.z", i, j, k];
			if (limitRegion && 
			    (x < x0 || y < y0 || z < z0 ||
			     x > x1 || y > y1 || z > z1)) continue;
			double volume = flowBlocks[ib]["volume", i, j, k];
			double value = flowBlocks[ib][varName, i, j, k];
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

} // end class FlowSolution

class SBlockFlow {
    // Much like the Python library for postprocessing in Eilmer3,
    // we are going to handle the data as a big chunk of numbers,
    // with the label for each variable coming the top of the file.
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

    void add_aux_variables(string[] addVarsList)
    // Adds variables to the data for each cell.
    {
	// We assume a lot about the data that has been read in so,
	// we need to skip this function if all is not in place
	bool ok_to_proceed = true;
	foreach (name; ["a", "rho", "p", "vel.x", "vel.y", "vel.z", "e[0]", "tke"]) {
	    if (!canFind(variableNames, name)) { ok_to_proceed = false; }
	}
	if (!ok_to_proceed) {
	    writeln("SBlockFlow.add_aux_variables(): Some essential variables not found.");
	    return;
	}
	bool add_mach = canFind(addVarsList, "mach");
	bool add_pitot_p = canFind(addVarsList, "pitot");
	bool add_total_p = canFind(addVarsList, "total-p");
	bool add_total_h = canFind(addVarsList, "total-h");
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
	//
	// Be careful to add auxiliary variable values in the code below 
	// in the same order as the list of variable names in the code above.
	//
	foreach (k; 0 .. nkc) {
	    foreach (j; 0 .. njc) {
		foreach (i; 0 .. nic) {
		    double a = _data[i][j][k][variableIndex["a"]];
		    double p = _data[i][j][k][variableIndex["p"]];
		    double rho = _data[i][j][k][variableIndex["rho"]];
                    double g = a*a*rho/p; // approximation for gamma
                    // Velocity in the block frame of reference that
		    // may be rotating for turbomachinery calculations.
		    double wx = _data[i][j][k][variableIndex["vel.x"]];
		    double wy = _data[i][j][k][variableIndex["vel.y"]];
		    double wz = _data[i][j][k][variableIndex["vel.z"]];
                    double w = sqrt(wx*wx + wy*wy + wz*wz);
                    double M = w/a;
 		    if (add_mach) { _data[i][j][k] ~= M; }
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
                        _data[i][j][k] ~= pitot_p;
		    }
                    if (add_total_p) {
                        // Isentropic process only.
                        double t1 = 1 + 0.5*(g-1)*M*M;
                        double total_p = p * pow(t1,(g/(g-1)));
                        _data[i][j][k] ~= total_p;
		    }
                    if (add_total_h) {
                        double e0 = _data[i][j][k][variableIndex["e[0]"]];
                        double tke = _data[i][j][k][variableIndex["tke"]];
                        // Sum up the bits of energy,
                        // forgetting the multiple energy modes, for the moment.
                        double total_h = p/rho + e0 + 0.5*w*w + tke;
                        _data[i][j][k] ~= total_h;
		    }
 		} // foreach i
	    } // foreach j
	} // foreach k
	
    } // end add_aux_variables()


    void subtract_ref_soln(lua_State* L)
    {
	string luaFnName = "refSoln";
	foreach (k; 0 .. nkc) {
	    foreach (j; 0 .. njc) {
		foreach (i; 0 .. nic) {
		    // Call back to the Lua function to get a table of values.
		    // function refSoln(x, y, z)
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
			string errMsg = `Error in FlowSolution.subtract_ref_soln().;
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
} // end class SBlockFlow
