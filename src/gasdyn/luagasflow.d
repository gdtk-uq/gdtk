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
    // Function signature in Lua domain:
    // V2, Vg, state2 = gasflow.normal_shock(state1, Vs, rho_tol, T_tol)
    // Input:
    //   state1: in a GasState table (with gasmodel field)
    //   Vs: velocity of shock into quiescent gas 
    //   rho_tol: optional toleralnce on density
    //   T_tol: optional tolerance on temperature
    // Returns:
    //   V2: gas velocity leaving shock (in shock frame)
    //   Vg: gas velocity in lab frame, for a moving shock
    //   state2: gas state following shock
    //
    lua_getfield(L, 1, "gasmodel");
    GasModel gm = checkGasModel(L, -1);
    lua_pop(L, 1);
    GasState state1 = new GasState(gm);
    getGasStateFromTable(L, gm, 1, state1);
    gm.update_thermo_from_pT(state1); // needed for cea_gas
    // Same values into state2, for now.
    GasState state2 =  new GasState(gm);
    getGasStateFromTable(L, gm, 1, state2);
    gm.update_thermo_from_pT(state2);
    //
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
    //
    double[] vel_results = normal_shock(state1, Vs, state2, gm, rho_tol, T_tol);
    //
    lua_settop(L, 0); // clear the stack, in preparation for pushing results
    lua_pushnumber(L, vel_results[0]); // V2
    lua_pushnumber(L, vel_results[1]); // Vg
    pushNewGasTable(L, state2, gm);
    return 3;
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
