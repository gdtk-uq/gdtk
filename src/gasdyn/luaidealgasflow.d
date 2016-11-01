/**
 * A Lua interface for the D idealgasflow module.
 *
 * Authors: Peter J. and Rowan G.
 * Date: 2016-10-16, just enough for the Billig shock shape correlation
 */

module luaidealgasflow;

import std.stdio;
import std.string;
import std.conv;
import util.lua;
import util.lua_service;
import idealgasflow;

// Name of metatable
immutable string idealgasflowMT = "idealgasflow";

extern(C) int idealgasflow_A_Astar(lua_State* L)
{
    if (!lua_isnumber(L, 1)) {
	string errMsg = "Expected a number for M";
	luaL_error(L, errMsg.toStringz);
    }
    double mach = to!double(luaL_checknumber(L, 1));
    double g = 1.4; // default value
    if (lua_isnumber(L, 2)) {
	g = to!double(luaL_checknumber(L, 2));
    }
    lua_pushnumber(L, A_Astar(mach, g));
    return 1;
}

extern(C) int idealgasflow_beta_obl(lua_State* L)
{
    if (!lua_isnumber(L, 1)) {
	string errMsg = "Expected a number for M1";
	luaL_error(L, errMsg.toStringz);
    }
    double M1 = to!double(luaL_checknumber(L, 1));
    if (!lua_isnumber(L, 2)) {
	string errMsg = "Expected a number for theta";
	luaL_error(L, errMsg.toStringz);
    }
    double theta = to!double(luaL_checknumber(L, 2));
    double g = 1.4; // default value
    if (lua_isnumber(L, 3)) {
	g = to!double(luaL_checknumber(L, 3));
    }
    double tol = 1.0e-6; // default value
    if (lua_isnumber(L, 4)) {
	tol = to!double(luaL_checknumber(L, 4));
    }
    lua_pushnumber(L, beta_obl(M1, theta, g, tol));
    return 1;
}

extern(C) int idealgasflow_M2_obl(lua_State* L)
{
    if (!lua_isnumber(L, 1)) {
	string errMsg = "Expected a number for M1";
	luaL_error(L, errMsg.toStringz);
    }
    double M1 = to!double(luaL_checknumber(L, 1));
    if (!lua_isnumber(L, 2)) {
	string errMsg = "Expected a number for beta";
	luaL_error(L, errMsg.toStringz);
    }
    double beta = to!double(luaL_checknumber(L, 2));
    if (!lua_isnumber(L, 3)) {
	string errMsg = "Expected a number for theta";
	luaL_error(L, errMsg.toStringz);
    }
    double theta = to!double(luaL_checknumber(L, 3));
    double g = 1.4; // default value
    if (lua_isnumber(L, 4)) {
	g = to!double(luaL_checknumber(L, 4));
    }
    lua_pushnumber(L, M2_obl(M1, beta, theta, g));
    return 1;
}

extern(C) int idealgasflow_p2_p1_obl(lua_State* L)
{
    if (!lua_isnumber(L, 1)) {
	string errMsg = "Expected a number for M1";
	luaL_error(L, errMsg.toStringz);
    }
    double M1 = to!double(luaL_checknumber(L, 1));
    if (!lua_isnumber(L, 2)) {
	string errMsg = "Expected a number for beta";
	luaL_error(L, errMsg.toStringz);
    }
    double beta = to!double(luaL_checknumber(L, 2));
    double g = 1.4; // default value
    if (lua_isnumber(L, 3)) {
	g = to!double(luaL_checknumber(L, 3));
    }
    lua_pushnumber(L, p2_p1_obl(M1, beta, g));
    return 1;
}


void registeridealgasflowFunctions(lua_State* L)
{
    // Register the idealgasflow table of functions.
    luaL_newmetatable(L, idealgasflowMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &idealgasflow_A_Astar);
    lua_setfield(L, -2, "A_Astar");
    lua_pushcfunction(L, &idealgasflow_beta_obl);
    lua_setfield(L, -2, "beta_obl");
    lua_pushcfunction(L, &idealgasflow_M2_obl);
    lua_setfield(L, -2, "M2_obl");
    lua_pushcfunction(L, &idealgasflow_p2_p1_obl);
    lua_setfield(L, -2, "p2_p1_obl");

    lua_setglobal(L, idealgasflowMT.toStringz);
}
    






