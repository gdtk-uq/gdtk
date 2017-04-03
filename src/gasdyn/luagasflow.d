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
import gas.cea_gas;
import gas.luagas_model;
import idealgasflow;
import gasflow;

// Name of metatable
immutable string gasflowMT = "gasflow";

extern(C) int gasflow_normal_shock(lua_State* L)
{
    // Function signature in Lua domain:
    // state2, V2, Vg = gasflow.normal_shock(state1, Vs, rho_tol, T_tol)
    // Input:
    //   state1: a GasState table (with gasmodel field) for pre-shock state
    //   Vs: velocity of shock into quiescent gas 
    //   rho_tol: optional tolerance on density
    //   T_tol: optional tolerance on temperature
    // Returns:
    //   state2: a GasState table for the gas state following the shock
    //   V2: gas velocity leaving shock (in shock frame)
    //   Vg: gas velocity in lab frame, for a moving shock
    //
    lua_getfield(L, 1, "gasmodel");
    GasModel gm = checkGasModel(L, -1);
    lua_pop(L, 1);
    GasState state1 = new GasState(gm);
    getGasStateFromTable(L, gm, 1, state1);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state1); }
    // Same values into state2, for now.
    GasState state2 =  new GasState(gm);
    getGasStateFromTable(L, gm, 1, state2);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state2); }
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
    pushNewGasTable(L, state2, gm);
    lua_pushnumber(L, vel_results[0]); // V2
    lua_pushnumber(L, vel_results[1]); // Vg
    return 3;
} // end gasflow_normal_shock()

extern(C) int gasflow_normal_shock_p2p1(lua_State* L)
{
    // Function signature in Lua domain:
    // V1, V2, Vg = gasflow.normal_shock_p2p1(state1, p2p1)
    // Input:
    //   state1: a GasState table for pre-shock gas state
    //   p2p1: ratio of pressure across the shock 
    // Returns:
    //   V1: the incident shock speed (into quiescent gas)
    //   V2: gas velocity leaving shock (in shock frame)
    //   Vg: gas velocity in lab frame, for a moving shock
    //
    lua_getfield(L, 1, "gasmodel");
    GasModel gm = checkGasModel(L, -1);
    lua_pop(L, 1);
    GasState state1 = new GasState(gm);
    getGasStateFromTable(L, gm, 1, state1);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state1); }
    // Same values into state2, for now.
    GasState state2 =  new GasState(gm);
    getGasStateFromTable(L, gm, 1, state2);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state2); }
    //
    if (!lua_isnumber(L, 2)) {
	string errMsg = "Expected a number for p2p1";
	luaL_error(L, errMsg.toStringz);
    }
    double p2p1 = to!double(luaL_checknumber(L, 2));
    //
    double[] vel_results = normal_shock_p2p1(state1, p2p1, state2, gm);
    //
    lua_settop(L, 0); // clear the stack, in preparation for pushing results
    lua_pushnumber(L, vel_results[0]); // V1
    lua_pushnumber(L, vel_results[1]); // V2
    lua_pushnumber(L, vel_results[2]); // Vg
    return 3;
} // end gasflow_normal_shock_p2p1()

extern(C) int gasflow_reflected_shock(lua_State* L)
{
    // Function signature in Lua domain:
    // state5, Vr = gasflow.reflected_shock(state2, Vg)
    // Input:
    //   state2: a GasState table for gas following incident shock
    //   Vg: velocity (in lab frame) of gas following the incident shock 
    // Returns:
    //   state5: gas state between reflected shock and tube end
    //   Vr: velocity (in lab frame) of the shock moving upstream
    //
    lua_getfield(L, 1, "gasmodel");
    GasModel gm = checkGasModel(L, -1);
    lua_pop(L, 1);
    GasState state2 = new GasState(gm);
    getGasStateFromTable(L, gm, 1, state2);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state2); }
    // Same values into state5, for now.
    GasState state5 =  new GasState(gm);
    getGasStateFromTable(L, gm, 1, state5);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state5); }
    //
    if (!lua_isnumber(L, 2)) {
	string errMsg = "Expected a number for Vg";
	luaL_error(L, errMsg.toStringz);
    }
    double Vg = to!double(luaL_checknumber(L, 2));
    //
    double Vr = reflected_shock(state2, Vg, state5, gm);
    //
    lua_settop(L, 0); // clear the stack, in preparation for pushing results
    pushNewGasTable(L, state5, gm);
    lua_pushnumber(L, Vr);
    return 2;
} // end gasflow_reflected_shock()

extern(C) int gasflow_expand_from_stagnation(lua_State* L)
{
    // Function signature in Lua domain:
    // state1, V = gasflow.expand_from_stagnation(state0, p_over_p0)
    // Input:
    //   state0: a GasState table for stagnation gas
    //   p_over_p0: pressure ratio
    // Returns:
    //   state1: GasState table for expanded gas 
    //   V: velocity of expanded gas
    //
    lua_getfield(L, 1, "gasmodel");
    GasModel gm = checkGasModel(L, -1);
    lua_pop(L, 1);
    GasState state0 = new GasState(gm);
    getGasStateFromTable(L, gm, 1, state0);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state0); }
    // Same values into state1, for now.
    GasState state1 =  new GasState(gm);
    getGasStateFromTable(L, gm, 1, state1);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state1); }
    //
    if (!lua_isnumber(L, 2)) {
	string errMsg = "Expected a number for p_over_p0";
	luaL_error(L, errMsg.toStringz);
    }
    double p_over_p0 = to!double(luaL_checknumber(L, 2));
    //
    double V = expand_from_stagnation(p_over_p0, state0, state1, gm);
    //
    lua_settop(L, 0); // clear the stack, in preparation for pushing results
    pushNewGasTable(L, state1, gm);
    lua_pushnumber(L, V);
    return 2;
} // end gasflow_expand_from_stagnation()

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
    mixin(registerfn("normal_shock_p2p1"));
    mixin(registerfn("reflected_shock"));
    mixin(registerfn("expand_from_stagnation"));

    lua_setglobal(L, gasflowMT.toStringz);
} // end registergasflowFunctions()
