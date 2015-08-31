/**
 * luasolidprops.d
 * Lua interface to access writeInitialSolidFile
 *
 * Authors: Rowan G. and Peter J.
 * Date: 2015-04-29
 */

module luasolidprops;

import std.array;
import std.format;
import std.stdio;
import std.conv;
import std.string;

import util.lua;
import util.lua_service;
import globalconfig;
import luasgrid;
import solidfvcell;
import solidprops;

extern(C) int writeInitialSolidFileFromLua(lua_State* L)
{
    auto fname = to!string(luaL_checkstring(L, 1));
    auto grid = checkStructuredGrid(L, 2);
    double rho = getDouble(L, 4, "rho");
    double kS = getDouble(L, 4, "k");
    double Cp = getDouble(L, 4, "Cp");
    auto sp = new SolidProps(rho, kS, Cp);
    double t0 = luaL_checknumber(L, 5);

    // First check we have a simple constant T value to use
    if ( lua_isnumber(L, 3) ) {
	double T_init = luaL_checknumber(L, 3);
	writeInitialSolidFile(fname, grid, T_init, sp, t0);
	return 0;
    }
    // Else we might have a function to use.
    if ( lua_isfunction(L, 3) ) {
	// We'll assume the function is good to use.
	// Numbers of cells derived from numbers of vertices.
	int nic = grid.niv - 1;
	int njc = grid.njv - 1;
	int nkc = grid.nkv - 1;
	if (GlobalConfig.dimensions == 2) nkc = 1;
	
	string cellDataToString(int i, int j, int k)
	{
	    auto p000 = grid[i,j,k];
	    auto p100 = grid[i+1,j,k];
	    auto p110 = grid[i+1,j+1,k];
	    auto p010 = grid[i,j+1,k];
	    // [TODO] provide better calculation using geom module.
	    // For the moment, it doesn't matter greatly because the solver 
	    // will compute it's own approximations
	    auto pos = 0.25*(p000 + p100 + p110 + p010);
	    auto volume = 0.0;
	    if (GlobalConfig.dimensions == 3) {
		auto p001 = grid[i,j,k+1];
		auto p101 = grid[i+1,j,k+1];
		auto p111 = grid[i+1,j+1,k+1];
		auto p011 = grid[i,j+1,k+1];
	    pos = 0.5*pos + 0.125*(p001 + p101 + p111 + p011);
	    }
	    // Now grab temperature via Lua function call
	    lua_pushvalue(L, 3);
	    lua_pushnumber(L, pos.x);
	    lua_pushnumber(L, pos.y);
	    lua_pushnumber(L, pos.z);
	    if ( lua_pcall(L, 3, 1, 0) != 0 ) {
		string errMsg = "Error in Lua function call for setting temperature\n";
		errMsg ~= "as a function of postion (x, y, z).\n";
		luaL_error(L, errMsg.toStringz);
	    }
	    double T = luaL_checknumber(L, -1);
	    double e = updateEnergy(sp, T);
	    lua_pop(L, 1);

	    // Should match SolidFVCell.writeValuesToString()
	    auto writer = appender!string();
	    formattedWrite(writer, "%.16e %.16e %.16e %.16e %.16e %.16e",
			   pos.x, pos.y, pos.z, volume, e, T);
	    
	    return writer.data;
	}

	// Write the data for the whole structured block.
	auto f = File(fname, "w");
	f.writefln("%20.16e", t0);
	// Variable list for cell on one line.
	auto writer = appender!string();
	foreach(varname; varListForSolidCell()) {
	    formattedWrite(writer, " \"%s\"", varname);
	}
	f.writeln(writer.data);
	// Numbers of cells.
	f.writefln("%d %d %d", nic, njc, nkc);
	// The actual cell data.
	foreach (k; 0 .. nkc) {
	    foreach (j; 0 .. njc) {
		foreach (i; 0 .. nic) {
		    f.writefln(" %s", cellDataToString(i,j,k));
		}
	    }
	}
	f.close();
	return 0;
    }
    // Something went wrong if we get here.
    return -1;
}

void registerSolidProps(lua_State* L)
{
    lua_pushcfunction(L, &writeInitialSolidFileFromLua);
    lua_setglobal(L, "writeInitialSolidFile");
}
