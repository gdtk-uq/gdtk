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
import geom;
import globalconfig;
import geom.luawrap;
import solidfvcell;
import solidprops;

extern(C) int writeInitialSolidFileFromLua(lua_State* L)
{
    auto fname = to!string(luaL_checkstring(L, 1));
    auto grid = checkStructuredGrid(L, 2);
    double t0 = luaL_checknumber(L, 5);
    SolidProps sp;
    // First check if we have simple constant values for the
    // initial temperature, and the properties.
    if ( lua_isnumber(L, 3) && lua_istable(L, 4) ) {
	if (!GlobalConfig.solid_has_homogeneous_properties) {
	    string errMsg = "\nERROR: solid has been declared as having inhomogeneous properties,\n" ~
		"and yet a table of constant material properties has been supplied.\n" ~
		"Instead, a function of (x,y,z) is expected.\n";
	    luaL_error(L, errMsg.toStringz);
	}
	double T_init = luaL_checknumber(L, 3);
	double rho = getDouble(L, 4, "rho");
	double Cp = getDouble(L, 4, "Cp");
	double kS = 0.0;
	if (GlobalConfig.solid_has_isotropic_properties) {
	    kS = getDouble(L, 4, "k");
	    sp = SolidProps(rho, kS, Cp);
	}
	else {
	    double k11 = 0.0;
            double k12 = 0.0;
	    double k13 = 0.0;
	    double k21 = 0.0; 
	    double k22 = 0.0;
	    double k23 = 0.0;
	    double k31 = 0.0;
	    double k32 = 0.0;
	    double k33 = 0.0;
	    k11 = getDouble(L, 4, "k11");
	    k12 = getDouble(L, 4, "k12");
	    k21 = getDouble(L, 4, "k21");
	    k22 = getDouble(L, 4, "k22");
	    if (GlobalConfig.dimensions == 3) {
		k13 = getDouble(L, 4, "k13");
		k23 = getDouble(L, 4, "k23");
		k31 = getDouble(L, 4, "k31");
		k32 = getDouble(L, 4, "k32");
		k33 = getDouble(L, 4, "k33");
	    }
	    sp = SolidProps(rho, kS, Cp,
			    k11, k12, k13,
			    k21, k22, k23,
			    k31, k32, k33);
	}
	writeInitialSolidFile(fname, grid, T_init, sp, t0);
	return 0;
    }
    // Else, we have a function to use for either the initial
    // temperature, or the material properties, or both.
    // In either case, we need to write out the file from this
    // function, that is, we can't delegate.
    if (GlobalConfig.solid_has_homogeneous_properties) {
	if (!lua_istable(L, 4)) {
	    string errMsg = "\nERROR: The solid domain is set with homogeneous properties,\n" ~
		"so a table of material 'properties' was expected, but not found.\n";
	    luaL_error(L, errMsg.toStringz);
	}
	// Else, we can setup the SolidProps object
	double rho = getDouble(L, 4, "rho");
	double Cp = getDouble(L, 4, "Cp");
	double kS = 0.0;
	if (GlobalConfig.solid_has_isotropic_properties) {
	    kS = getDouble(L, 4, "k");
	    sp = SolidProps(rho, kS, Cp);
	}
	else {
	    double k11 = 0.0;
            double k12 = 0.0;
	    double k13 = 0.0;
	    double k21 = 0.0; 
	    double k22 = 0.0;
	    double k23 = 0.0;
	    double k31 = 0.0;
	    double k32 = 0.0;
	    double k33 = 0.0;
	    k11 = getNumberFromTable(L, 4, "k11", true, 0.0, true, "\nproblem with k11 value");
	    k12 = getNumberFromTable(L, 4, "k12", true, 0.0, true, "\nproblem with k12 value");
	    k21 = getNumberFromTable(L, 4, "k21", true, 0.0, true, "\nproblem with k21 value");
	    k22 = getNumberFromTable(L, 4, "k22", true, 0.0, true, "\nproblem with k21 value");
	    if (GlobalConfig.dimensions == 3) {
		k13 = getNumberFromTable(L, 4, "k13", true, 0.0, true, "\nproblem with k13 value");
		k23 = getNumberFromTable(L, 4, "k23", true, 0.0, true, "\nproblem with k23 value");
		k31 = getNumberFromTable(L, 4, "k31", true, 0.0, true, "\nproblem with k31 value");
		k32 = getNumberFromTable(L, 4, "k32", true, 0.0, true, "\nproblem with k31 value");
		k33 = getNumberFromTable(L, 4, "k33", true, 0.0, true, "\nproblem with k33 value");
	    }
	    sp = SolidProps(rho, kS, Cp,
			    k11, k12, k13,
			    k21, k22, k23,
			    k31, k32, k33);
	}
    }
    else {
	// Check we have a function for later use.
	if (!lua_isfunction(L, 4) ) {
	    string errMsg = "ERROR: solid has been declared as having inhomogeneous properties,\n" ~
		"but no function of (x,y,z) as been supplied.\n";
	    luaL_error(L, errMsg.toStringz);
	}
    }
    // Let's check on the temperature value.
    double T_init = -1.0;
    if ( lua_isnumber(L, 3) ) {
	T_init = luaL_checknumber(L, 3);
    }
    
    auto nic = grid.niv - 1;
    auto njc = grid.njv - 1;
    auto nkc = grid.nkv - 1;
    if (GlobalConfig.dimensions == 2) nkc = 1;
    
    string cellDataToString(size_t i, size_t j, size_t k)
    {
	Vector3 p000 = *grid[i,j,k];
	Vector3 p100 = *grid[i+1,j,k];
	Vector3 p110 = *grid[i+1,j+1,k];
	Vector3 p010 = *grid[i,j+1,k];
	// [TODO] provide better calculation using geom module.
	// For the moment, it doesn't matter greatly because the solver 
	// will compute its own approximations
	auto pos = 0.25*(p000 + p100 + p110 + p010);
	auto volume = 0.0;
	if (GlobalConfig.dimensions == 3) {
	    Vector3 p001 = *grid[i,j,k+1];
	    Vector3 p101 = *grid[i+1,j,k+1];
	    Vector3 p111 = *grid[i+1,j+1,k+1];
	    Vector3 p011 = *grid[i,j+1,k+1];
	    pos = 0.5*pos + 0.125*(p001 + p101 + p111 + p011);
	}
	if ( lua_isfunction(L, 4) ) {
	    // Now grab material properties via Lua function call
	    lua_pushvalue(L, 4);
	    lua_pushnumber(L, pos.x);
	    lua_pushnumber(L, pos.y);
	    lua_pushnumber(L, pos.z);
	    if ( lua_pcall(L, 3, 1, 0) != 0 ) {
		string errMsg = "\nError in Lua function call for setting materical properties\n";
		errMsg ~= "as a function of postion (x, y, z).\n";
		luaL_error(L, errMsg.toStringz);
	    }
	    if (!lua_istable(L, -1)) {
		string errMsg = "ERROR: The material properties function should return a table.\n";
		luaL_error(L, errMsg.toStringz);
	    }
	    
	    double rho = getDouble(L, -1, "rho");
	    double Cp = getDouble(L, -1, "Cp");
	    double kS = 0.0;
	    if (GlobalConfig.solid_has_isotropic_properties) {
		kS = getDouble(L, -1, "k");
		sp = SolidProps(rho, kS, Cp);
	    }
	    else {
		double k11 = 0.0;
		double k12 = 0.0;
		double k13 = 0.0;
		double k21 = 0.0; 
		double k22 = 0.0;
		double k23 = 0.0;
		double k31 = 0.0;
		double k32 = 0.0;
		double k33 = 0.0;
		k11 = getNumberFromTable(L, -1, "k11", true, 0.0, true, "\nproblem with k11 value");
		k12 = getNumberFromTable(L, -1, "k12", true, 0.0, true, "\nproblem with k12 value");
		k21 = getNumberFromTable(L, -1, "k21", true, 0.0, true, "\nproblem with k21 value");
		k22 = getNumberFromTable(L, -1, "k22", true, 0.0, true, "\nproblem with k21 value");
		if (GlobalConfig.dimensions == 3) {
		    k13 = getNumberFromTable(L, -1, "k13", true, 0.0, true, "\nproblem with k13 value");
		    k23 = getNumberFromTable(L, -1, "k23", true, 0.0, true, "\nproblem with k23 value");
		    k31 = getNumberFromTable(L, -1, "k31", true, 0.0, true, "\nproblem with k31 value");
		    k32 = getNumberFromTable(L, -1, "k32", true, 0.0, true, "\nproblem with k31 value");
		    k33 = getNumberFromTable(L, -1, "k33", true, 0.0, true, "\nproblem with k33 value");
		}
		sp = SolidProps(rho, kS, Cp,
				k11, k12, k13,
				k21, k22, k23,
				k31, k32, k33);
	    }
	    lua_pop(L, 1);
	}
	double T, e;
	if ( lua_isfunction(L, 3) ) {
	    // Now grab temperature via Lua function call
	    lua_pushvalue(L, 3);
	    lua_pushnumber(L, pos.x);
	    lua_pushnumber(L, pos.y);
	    lua_pushnumber(L, pos.z);
	    if ( lua_pcall(L, 3, 1, 0) != 0 ) {
		string errMsg = "\nError in Lua function call for setting temperature\n";
		errMsg ~= "as a function of postion (x, y, z).\n";
		luaL_error(L, errMsg.toStringz);
	    }
	    T = luaL_checknumber(L, -1);
	    e = updateEnergy(sp, T);
	    lua_pop(L, 1);
	}
	else {
	    T = T_init;
	    e = updateEnergy(sp, T);
	}
	// Should match SolidFVCell.writeValuesToString()
	auto writer = appender!string();
	formattedWrite(writer, "%.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e",
		       pos.x, pos.y, pos.z, volume, e, T,
		       sp.rho, sp.Cp, sp.k,
		       sp.k11, sp.k12, sp.k13,
		       sp.k21, sp.k22, sp.k23,
		       sp.k31, sp.k32, sp.k33);
	return writer.data;
    }

    // Write the data for the whole structured block.
    auto f = File(fname, "w");
    f.writefln("%.18e", t0);
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

void registerSolidProps(lua_State* L)
{
    lua_pushcfunction(L, &writeInitialSolidFileFromLua);
    lua_setglobal(L, "writeInitialSolidFile");
}
