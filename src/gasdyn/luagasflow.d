/**
 * A Lua interface for the D gasflow module.
 *
 * Authors: Peter J. and Rowan G.
 * Date: 2017-04-02, just enough to begin the nenzfr2 demo
 */

module luagasflow;

import std.stdio;
import std.string;
import std.conv;
import std.algorithm;
import util.lua;
import util.lua_service;
import gas.gas_model;
import gas.luagas_model;
import idealgasflow;
import gasflow;

// Name of metatable
immutable string gasflowMT = "gasflow";

extern(C) int gasflow_normal_shock(lua_State* L)
{
    // Expected arguments:
    //   1 state1 table
    //   2 Vs
    //   3 rho_tol (optional)
    //   4 T_tol (optional)
    // We also expect the gasmodel filed in the state1 table.
    //
    // [TODO] consider accepting everything in a table
    // and returning all results in a table.
    //
    lua_getfield(L, 1, "gasmodel");
    GasModel gm = checkGasModel(L, -1);
    lua_pop(L, 1);
    GasState state1 = new GasState(gm);
    getGasStateFromTable(L, gm, 1, state1);
    writeln("state1: ", state1.toString);
    if (!lua_isnumber(L, 2)) {
	string errMsg = "Expected a number for Vs";
	luaL_error(L, errMsg.toStringz);
    }
    double Vs = to!double(luaL_checknumber(L, 2));
    double rho_tol=1.0e-6; // default value
    if (lua_isnumber(L, 3)) {
	rho_tol = to!double(luaL_checknumber(L, 3));
    }
    double T_tol = 0.1; // default value
    if (lua_isnumber(L, 4)) {
	T_tol = to!double(luaL_checknumber(L, 4));
    }
    GasState state2 =  new GasState(state1);
    double[] results = normal_shock(state1, Vs, state2, gm, rho_tol, T_tol);
    // [TODO] PJ, consider clearing stack at this point.
    lua_pushnumber(L, results[0]); // V2
    lua_pushnumber(L, results[1]); // Vg
    // [TODO] return a table as the third item
    // setGasStateInTable(L, gm, 3, state2);
    return 2;
} // end gasflow_normal_shock()

string registerfn(string fname)
{
    return "    lua_pushcfunction(L, &gasflow_"~fname~");\n" ~
	"    lua_setfield(L, -2, \""~fname~"\");";
}

void registergasflowFunctions(lua_State* L)
{
    // Register the gasflow table of functions.
    luaL_newmetatable(L, gasflowMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    mixin(registerfn("normal_shock"));

    lua_setglobal(L, gasflowMT.toStringz);
} // end registergasflowFunctions()
