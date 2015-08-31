/**
 * solid_udf_source_terms.d
 *
 * This module handles user-defined source terms
 * that the user might specify with a Lua script.
 *
 * Authors: RG & PJ
 * Date: 2015-05-06
 **/

module solid_udf_source_terms;

import std.stdio;
import std.string;
import util.lua;
import util.lua_service;

import solidfvcell;

void initUDFSolidSourceTerms(lua_State* L, string fname)
{
    luaL_dofile(L, fname.toStringz);
}

void addUDFSourceTermsToSolidCell(lua_State* L, SolidFVCell cell, double t)
{
    // Push user function onto TOS
    lua_getglobal(L, "solidSourceTerms");
    // Push sim_time onto TOS
    lua_pushnumber(L, t);
    // Push useful data into an arguments table
    lua_newtable(L);
    lua_pushnumber(L, cell.pos.x); lua_setfield(L, -2, "x");
    lua_pushnumber(L, cell.pos.y); lua_setfield(L, -2, "y");
    lua_pushnumber(L, cell.pos.z); lua_setfield(L, -2, "z");
    lua_pushnumber(L, cell.volume); lua_setfield(L, -2, "vol");
    // Call solidSourceTerms function with (args)
    int number_args = 2;
    int number_results = 1;

    if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
	luaL_error(L, "error running solid source terms function: %s\n",
		   lua_tostring(L, -1));
    }
    
    // Grab the energy source term
    cell.Q = luaL_checknumber(L, -1);
    lua_pop(L, 1);
}
