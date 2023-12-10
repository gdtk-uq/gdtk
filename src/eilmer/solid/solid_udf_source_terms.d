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
import std.conv;
import util.lua;
import util.lua_service;

import globalconfig;
import solidfvcell;
import ssolidblock;

void initUDFSolidSourceTerms(lua_State* L, string fname)
{
    luaL_dofile(L, fname.toStringz);
}

void addUDFSourceTermsToSolidCell(lua_State* L, SolidFVCell cell, double t, SSolidBlock blk)
{
    // Push user function onto TOS
    lua_getglobal(L, "solidSourceTerms");
    // Push sim_time onto TOS
    lua_pushnumber(L, t);
    // Push useful data into an arguments table
    lua_newtable(L);
    lua_pushnumber(L, cell.T); lua_setfield(L, -2, "T");
    lua_pushnumber(L, cell.pos.x); lua_setfield(L, -2, "x");
    lua_pushnumber(L, cell.pos.y); lua_setfield(L, -2, "y");
    lua_pushnumber(L, cell.pos.z); lua_setfield(L, -2, "z");
    lua_pushnumber(L, cell.volume); lua_setfield(L, -2, "vol");
    // We want cell indices in the Lua domain to look like cell indices for a SFluidBlock.
    auto ijk = blk.toIJKIndices(cell.id);
    int i = to!int(ijk[0]);
    int j = to!int(ijk[1]);
    int k = 0;
    if (GlobalConfig.dimensions == 3) {
        k = to!int(ijk[2]);
    }
    lua_pushinteger(L, i); lua_setfield(L, -2, "i");
    lua_pushinteger(L, j); lua_setfield(L, -2, "j");
    lua_pushinteger(L, k); lua_setfield(L, -2, "k");
    lua_pushinteger(L, blk.id); lua_setfield(L, -2, "blkId");
    // Call solidSourceTerms function with (args)
    int number_args = 2;
    int number_results = 1;
    //
    if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
        luaL_error(L, "error running solid source terms function: %s\n",
                   lua_tostring(L, -1));
    }
    //
    // Grab the energy source term
    cell.Q = luaL_checknumber(L, -1);
    lua_pop(L, 1);
} // end addUDFSourceTermsToSolidCell()
