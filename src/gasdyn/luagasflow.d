/**
 * A Lua interface for the D gasflow module.
 *
 * Authors: Peter J. and Rowan G.
 * Date: 2017-04-02, just enough to begin the nenzfr2 demo
 */

module gasdyn.luagasflow;

import std.stdio;
import std.string;
import std.conv;
import std.algorithm;
import ntypes.complex;
import nm.number;
import util.lua;
import util.lua_service;
import gas;
import gas.cea_gas;
import gas.luagas_model;
import gasdyn.idealgasflow;
import gasdyn.gasflow;

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
    GasState state1 = GasState(gm);
    getGasStateFromTable(L, gm, 1, state1);
    // The CEAgas model is a bit special in that some of the state data
    // are stored in a table within the GasModel object.
    // We need to call the update_thermo function to get this table
    // up-to-date.
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state1); }
    // Same values into state2, for now.
    GasState state2 =  GasState(gm);
    getGasStateFromTable(L, gm, 1, state2);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state2); }
    //
    if (!lua_isnumber(L, 2)) {
        string errMsg = "Expected a number for Vs";
        luaL_error(L, errMsg.toStringz);
    }
    number Vs = to!number(luaL_checknumber(L, 2));
    double rho_tol=1.0e-6; // default value
    if (lua_isnumber(L, 3)) {
        rho_tol = to!double(luaL_checknumber(L, 3));
    }
    double T_tol = 0.1; // default value
    if (lua_isnumber(L, 4)) {
        T_tol = to!double(luaL_checknumber(L, 4));
    }
    try {
        number[] vel_results = normal_shock(state1, Vs, state2, gm, rho_tol, T_tol);
        lua_settop(L, 0); // clear the stack, in preparation for pushing results
        pushNewGasTable(L, state2, gm);
        lua_pushnumber(L, vel_results[0]); // V2
        lua_pushnumber(L, vel_results[1]); // Vg
        return 3;
    } catch(Exception e) {
        return luaL_error(L, e.msg.toStringz);
    }
} // end gasflow_normal_shock()

extern(C) int gasflow_normal_shock_1(lua_State* L)
{
    // Function signature in Lua domain:
    // state2, V2, Vg = gasflow.normal_shock_1(state1, Vs, rho_tol, T_tol)
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
    GasState state1 = GasState(gm);
    getGasStateFromTable(L, gm, 1, state1);
    // The CEAgas model is a bit special in that some of the state data
    // are stored in a table within the GasModel object.
    // We need to call the update_thermo function to get this table
    // up-to-date.
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state1); }
    // Same values into state2, for now.
    GasState state2 =  GasState(gm);
    getGasStateFromTable(L, gm, 1, state2);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state2); }
    //
    if (!lua_isnumber(L, 2)) {
        string errMsg = "Expected a number for Vs";
        luaL_error(L, errMsg.toStringz);
    }
    number Vs = to!number(luaL_checknumber(L, 2));
    double rho_tol=1.0e-6; // default value
    if (lua_isnumber(L, 3)) {
        rho_tol = to!double(luaL_checknumber(L, 3));
    }
    double T_tol = 0.1; // default value
    if (lua_isnumber(L, 4)) {
        T_tol = to!double(luaL_checknumber(L, 4));
    }
    try {
        number[] vel_results = normal_shock_1(state1, Vs, state2, gm, rho_tol, T_tol);
        lua_settop(L, 0); // clear the stack, in preparation for pushing results
        pushNewGasTable(L, state2, gm);
        lua_pushnumber(L, vel_results[0]); // V2
        lua_pushnumber(L, vel_results[1]); // Vg
        return 3;
    } catch(Exception e) {
        return luaL_error(L, e.msg.toStringz);
    }
} // end gasflow_normal_shock_1()

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
    GasState state1 = GasState(gm);
    getGasStateFromTable(L, gm, 1, state1);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state1); }
    // Same values into state2, for now.
    GasState state2 =  GasState(gm);
    getGasStateFromTable(L, gm, 1, state2);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state2); }
    //
    if (!lua_isnumber(L, 2)) {
        string errMsg = "Expected a number for p2p1";
        luaL_error(L, errMsg.toStringz);
    }
    number p2p1 = to!number(luaL_checknumber(L, 2));
    //
    try {
        number[] vel_results = normal_shock_p2p1(state1, p2p1, state2, gm);
        //
        lua_settop(L, 0); // clear the stack, in preparation for pushing results
        lua_pushnumber(L, vel_results[0]); // V1
        lua_pushnumber(L, vel_results[1]); // V2
        lua_pushnumber(L, vel_results[2]); // Vg
        return 3;
    } catch(Exception e) {
        return luaL_error(L, e.msg.toStringz);
    }
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
    GasState state2 = GasState(gm);
    getGasStateFromTable(L, gm, 1, state2);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state2); }
    // Same values into state5, for now.
    GasState state5 =  GasState(gm);
    getGasStateFromTable(L, gm, 1, state5);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state5); }
    //
    if (!lua_isnumber(L, 2)) {
        string errMsg = "Expected a number for Vg";
        luaL_error(L, errMsg.toStringz);
    }
    number Vg = to!number(luaL_checknumber(L, 2));
    //
    try {
        number Vr = reflected_shock(state2, Vg, state5, gm);
        //
        lua_settop(L, 0); // clear the stack, in preparation for pushing results
        pushNewGasTable(L, state5, gm);
        lua_pushnumber(L, Vr);
        return 2;
    } catch(Exception e) {
        return luaL_error(L, e.msg.toStringz);
    }
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
    GasState state0 = GasState(gm);
    getGasStateFromTable(L, gm, 1, state0);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state0); }
    // Same values into state1, for now.
    GasState state1 =  GasState(gm);
    getGasStateFromTable(L, gm, 1, state1);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state1); }
    //
    if (!lua_isnumber(L, 2)) {
        string errMsg = "Expected a number for p_over_p0";
        luaL_error(L, errMsg.toStringz);
    }
    number p_over_p0 = to!number(luaL_checknumber(L, 2));
    //
    try {
        number V = expand_from_stagnation(state0, p_over_p0, state1, gm);
        //
        lua_settop(L, 0); // clear the stack, in preparation for pushing results
        pushNewGasTable(L, state1, gm);
        lua_pushnumber(L, V);
        return 2;
    } catch(Exception e) {
        return luaL_error(L, e.msg.toStringz);
    }
} // end gasflow_expand_from_stagnation()

extern(C) int gasflow_expand_to_mach(lua_State* L)
{
    // Function signature in Lua domain:
    // state1, V = gasflow.expand_to_mach(state0, mach)
    // Input:
    //   state0: a GasState table for stagnation gas
    //   mach: mach number at expanded state
    // Returns:
    //   state1: GasState table for expanded gas
    //   V: velocity of expanded gas
    //
    lua_getfield(L, 1, "gasmodel");
    GasModel gm = checkGasModel(L, -1);
    lua_pop(L, 1);
    GasState state0 = GasState(gm);
    getGasStateFromTable(L, gm, 1, state0);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state0); }
    // Same values into state1, for now.
    GasState state1 =  GasState(gm);
    getGasStateFromTable(L, gm, 1, state1);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state1); }
    //
    if (!lua_isnumber(L, 2)) {
        string errMsg = "Expected a number for mach";
        luaL_error(L, errMsg.toStringz);
    }
    number mach = to!number(luaL_checknumber(L, 2));
    //
    try {
        number V = expand_to_mach(state0, mach, state1, gm);
        //
        lua_settop(L, 0); // clear the stack, in preparation for pushing results
        pushNewGasTable(L, state1, gm);
        lua_pushnumber(L, V);
        return 2;
    } catch(Exception e) {
        return luaL_error(L, e.msg.toStringz);
    }
} // end gasflow_expand_to_mach()

extern(C) int gasflow_total_condition(lua_State* L)
{
    // Function signature in Lua domain:
    // state0 = gasflow.total_condition(state1, V1)
    // Input:
    //   state1: a GasState table for free-stream gas
    //   V1: velocity of free-stream gas
    // Returns:
    //   state0: GasState table for stagnation condition
    //
    lua_getfield(L, 1, "gasmodel");
    GasModel gm = checkGasModel(L, -1);
    lua_pop(L, 1);
    GasState state1 = GasState(gm);
    getGasStateFromTable(L, gm, 1, state1);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state1); }
    // Same values into state0, for now.
    GasState state0 =  GasState(gm);
    getGasStateFromTable(L, gm, 1, state0);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state0); }
    //
    if (!lua_isnumber(L, 2)) {
        string errMsg = "Expected a number for V1";
        luaL_error(L, errMsg.toStringz);
    }
    number V1 = to!number(luaL_checknumber(L, 2));
    //
    try {
        total_condition(state1, V1, state0, gm);
        //
        lua_settop(L, 0); // clear the stack, in preparation for pushing results
        pushNewGasTable(L, state0, gm);
        return 1;
    } catch(Exception e) {
        return luaL_error(L, e.msg.toStringz);
    }
} // end gasflow_total_condition()

extern(C) int gasflow_pitot_condition(lua_State* L)
{
    // Function signature in Lua domain:
    // state2pitot = gasflow.pitot_condition(state1, V1)
    // Input:
    //   state1: a GasState table for free-stream gas
    //   V1: velocity of free-stream gas
    // Returns:
    //   state2pitot: GasState table for stagnation condition
    //
    lua_getfield(L, 1, "gasmodel");
    GasModel gm = checkGasModel(L, -1);
    lua_pop(L, 1);
    GasState state1 = GasState(gm);
    getGasStateFromTable(L, gm, 1, state1);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state1); }
    // Same values into state2pitot, for now.
    GasState state2pitot =  GasState(gm);
    getGasStateFromTable(L, gm, 1, state2pitot);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state2pitot); }
    //
    if (!lua_isnumber(L, 2)) {
        string errMsg = "Expected a number for V1";
        luaL_error(L, errMsg.toStringz);
    }
    number V1 = to!number(luaL_checknumber(L, 2));
    //
    try {
        pitot_condition(state1, V1, state2pitot, gm);
        //
        lua_settop(L, 0); // clear the stack, in preparation for pushing results
        pushNewGasTable(L, state2pitot, gm);
        return 1;
    } catch(Exception e) {
        return luaL_error(L, e.msg.toStringz);
    }
} // end gasflow_pitot_condition()

extern(C) int gasflow_steady_flow_with_area_change(lua_State* L)
{
    // Function signature in Lua domain:
    // state2, V2 = gasflow.steady_flow_with_area_change(state1, V1, A2_over_A1, tol)
    // Input:
    //   state1: a GasState table for condition at station 1
    //   V1: velocity of gas at station 1
    //   A2_over_A1: area ratio
    //   tol: (optional) tolerance for function solver
    // Returns:
    //   state2: GasState table for condition at station 2
    //   V2: velocity of gas at station 2
    //
    lua_getfield(L, 1, "gasmodel");
    GasModel gm = checkGasModel(L, -1);
    lua_pop(L, 1);
    GasState state1 = GasState(gm);
    getGasStateFromTable(L, gm, 1, state1);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state1); }
    // Same values into state2, for now.
    GasState state2 =  GasState(gm);
    getGasStateFromTable(L, gm, 1, state2);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state2); }
    //
    if (!lua_isnumber(L, 2)) {
        string errMsg = "Expected a number for V1";
        luaL_error(L, errMsg.toStringz);
    }
    number V1 = to!number(luaL_checknumber(L, 2));
    if (!lua_isnumber(L, 3)) {
        string errMsg = "Expected a number for A2_over_A1";
        luaL_error(L, errMsg.toStringz);
    }
    number A2_over_A1 = to!number(luaL_checknumber(L, 3));
    double tol=1.0e-4; // default value
    if (lua_isnumber(L, 4)) {
        tol = to!double(luaL_checknumber(L, 4));
    }
    double p2p1_min=1.0e-4; // default value
    if (lua_isnumber(L, 5)) {
        p2p1_min = to!double(luaL_checknumber(L, 5));
    }
    //
    try {
        number V2 = steady_flow_with_area_change(state1, V1, A2_over_A1, state2,
                                                 gm, tol, p2p1_min);
        //
        lua_settop(L, 0); // clear the stack, in preparation for pushing results
        pushNewGasTable(L, state2, gm);
        lua_pushnumber(L, V2);
        return 2;
    } catch(Exception e) {
        return luaL_error(L, e.msg.toStringz);
    }
} // end gasflow_steady_flow_with_area_change()

extern(C) int gasflow_finite_wave_dp(lua_State* L)
{
    // Function signature in Lua domain:
    // state2, V2 = gasflow.finite_wave_dp(state1, V1, characteristic, p2, steps)
    // Input:
    //   state1: a GasState table for condition at station 1
    //   V1: velocity of gas at station 1
    //   p2: pressure at target state 2, following processing
    //   characteristic: "cplus" or "cminus"
    //   steps: (optional) number of steps to take in pressure
    // Returns:
    //   state2: GasState table for condition at station 2
    //   V2: velocity of gas at station 2
    //
    lua_getfield(L, 1, "gasmodel");
    GasModel gm = checkGasModel(L, -1);
    lua_pop(L, 1);
    GasState state1 = GasState(gm);
    getGasStateFromTable(L, gm, 1, state1);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state1); }
    // Same values into state2, for now.
    GasState state2 =  GasState(gm);
    getGasStateFromTable(L, gm, 1, state2);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state2); }
    //
    if (!lua_isnumber(L, 2)) {
        string errMsg = "Expected a number for V1";
        luaL_error(L, errMsg.toStringz);
    }
    number V1 = to!number(luaL_checknumber(L, 2));
    if (!lua_isstring(L, 3)) {
        string errMsg = "Expected a string for characteristic";
        luaL_error(L, errMsg.toStringz);
    }
    string characteristic = to!string(luaL_checkstring(L, 3));
    if (!lua_isnumber(L, 4)) {
        string errMsg = "Expected a number for p2";
        luaL_error(L, errMsg.toStringz);
    }
    number p2 = to!number(luaL_checknumber(L, 4));
    int steps = 100; // default value
    if (lua_isnumber(L, 5)) {
        steps = to!int(luaL_checkint(L, 5));
    }
    //
    try {
        number V2 = finite_wave_dp(state1, V1, characteristic, p2, state2, gm, steps);
        //
        lua_settop(L, 0); // clear the stack, in preparation for pushing results
        pushNewGasTable(L, state2, gm);
        lua_pushnumber(L, V2);
        return 2;
    } catch(Exception e) {
        return luaL_error(L, e.msg.toStringz);
    }
} // end gasflow_finite_wave_dp()

extern(C) int gasflow_finite_wave_dv(lua_State* L)
{
    // Function signature in Lua domain:
    // state2, V2 = gasflow.finite_wave_dv(state1, V1, characteristic, V2_target,
    //                                     steps, Tmin)
    // Input:
    //   state1: a GasState table for condition at station 1
    //   V1: velocity of gas at station 1
    //   V2_target: velocity at target state 2, following processing
    //   characteristic: "cplus" or "cminus"
    //   steps: (optional) number of steps to take through the process
    //   Tmin: (optional) temperature (in Kelvin) below which we terminate the process
    // Returns:
    //   state2: GasState table for condition at station 2
    //   V2: velocity of gas at station 2
    //
    lua_getfield(L, 1, "gasmodel");
    GasModel gm = checkGasModel(L, -1);
    lua_pop(L, 1);
    GasState state1 = GasState(gm);
    getGasStateFromTable(L, gm, 1, state1);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state1); }
    // Same values into state2, for now.
    GasState state2 =  GasState(gm);
    getGasStateFromTable(L, gm, 1, state2);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state2); }
    //
    if (!lua_isnumber(L, 2)) {
        string errMsg = "Expected a number for V1";
        luaL_error(L, errMsg.toStringz);
    }
    number V1 = to!number(luaL_checknumber(L, 2));
    if (!lua_isstring(L, 3)) {
        string errMsg = "Expected a string for characteristic";
        luaL_error(L, errMsg.toStringz);
    }
    string characteristic = to!string(luaL_checkstring(L, 3));
    if (!lua_isnumber(L, 4)) {
        string errMsg = "Expected a number for V2_target";
        luaL_error(L, errMsg.toStringz);
    }
    number V2_target = to!number(luaL_checknumber(L, 4));
    int steps = 100; // default value
    if (lua_isnumber(L, 5)) {
        steps = to!int(luaL_checkint(L, 5));
    }
    double Tmin = 200.0; // default value
    if (lua_isnumber(L, 6)) {
        Tmin = to!double(luaL_checknumber(L, 6));
    }
    //
    try {
        number V2 = finite_wave_dv(state1, V1, characteristic, V2_target,
                                   state2, gm, steps, Tmin);
        //
        lua_settop(L, 0); // clear the stack, in preparation for pushing results
        pushNewGasTable(L, state2, gm);
        lua_pushnumber(L, V2);
        return 2;
    } catch(Exception e) {
        return luaL_error(L, e.msg.toStringz);
    }
} // end gasflow_finite_wave_dv()

extern(C) int gasflow_theta_oblique(lua_State* L)
{
    // Function signature in Lua domain:
    // state2, theta, V2 = gasflow.theta_oblique(state1, V1, beta)
    // Input:
    //   state1: a GasState table for condition before shock
    //   V1: velocity of gas before shock
    //   beta: shock wave angle (in radians) wrt stream direction
    // Returns:
    //   state2: GasState table for condition post-shock
    //   V2: post-shock speed of gas (in m/s)
    //   theta: stream deflection angle (in radians)
    //
    lua_getfield(L, 1, "gasmodel");
    GasModel gm = checkGasModel(L, -1);
    lua_pop(L, 1);
    GasState state1 = GasState(gm);
    getGasStateFromTable(L, gm, 1, state1);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state1); }
    // Same values into state2, for now.
    GasState state2 =  GasState(gm);
    getGasStateFromTable(L, gm, 1, state2);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state2); }
    //
    if (!lua_isnumber(L, 2)) {
        string errMsg = "Expected a number for V1";
        luaL_error(L, errMsg.toStringz);
    }
    number V1 = to!number(luaL_checknumber(L, 2));
    if (!lua_isnumber(L, 3)) {
        string errMsg = "Expected a number for beta";
        luaL_error(L, errMsg.toStringz);
    }
    number beta = to!number(luaL_checknumber(L, 3));
    //
    try {
        number[] results = theta_oblique(state1, V1, beta, state2, gm);
        //
        lua_settop(L, 0); // clear the stack, in preparation for pushing results
        pushNewGasTable(L, state2, gm);
        lua_pushnumber(L, results[0]); // theta
        lua_pushnumber(L, results[1]); // V2
        return 3;
    } catch(Exception e) {
        return luaL_error(L, e.msg.toStringz);
    }
} // end gasflow_theta_oblique()

extern(C) int gasflow_beta_oblique(lua_State* L)
{
    // Function signature in Lua domain:
    // beta = gasflow.theta_oblique(state1, V1, theta)
    // Input:
    //   state1: a GasState table for condition before shock
    //   V1: velocity of gas before shock
    //   theta: stream deflection angle (in radians)
    // Returns:
    //   beta: shock wave angle (in radians) wrt free-stream direction
    //
    lua_getfield(L, 1, "gasmodel");
    GasModel gm = checkGasModel(L, -1);
    lua_pop(L, 1);
    GasState state1 = GasState(gm);
    getGasStateFromTable(L, gm, 1, state1);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state1); }
    // Same values into state2, for now.
    GasState state2 =  GasState(gm);
    getGasStateFromTable(L, gm, 1, state2);
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(state2); }
    //
    if (!lua_isnumber(L, 2)) {
        string errMsg = "Expected a number for V1";
        luaL_error(L, errMsg.toStringz);
    }
    number V1 = to!number(luaL_checknumber(L, 2));
    if (!lua_isnumber(L, 3)) {
        string errMsg = "Expected a number for theta";
        luaL_error(L, errMsg.toStringz);
    }
    number theta = to!number(luaL_checknumber(L, 3));
    //
    try {
        number beta = beta_oblique(state1, V1, theta, gm);
        //
        lua_settop(L, 0); // clear the stack, in preparation for pushing results
        lua_pushnumber(L, beta);
        return 1;
    } catch(Exception e) {
        return luaL_error(L, e.msg.toStringz);
    }
} // end gasflow_beta_oblique()

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
    mixin(registerfn("normal_shock_1"));
    mixin(registerfn("normal_shock_p2p1"));
    mixin(registerfn("reflected_shock"));
    mixin(registerfn("expand_from_stagnation"));
    mixin(registerfn("expand_to_mach"));
    mixin(registerfn("total_condition"));
    mixin(registerfn("pitot_condition"));
    mixin(registerfn("steady_flow_with_area_change"));
    mixin(registerfn("finite_wave_dp"));
    mixin(registerfn("finite_wave_dv"));
    mixin(registerfn("theta_oblique"));
    mixin(registerfn("beta_oblique"));

    lua_setglobal(L, gasflowMT.toStringz);
} // end registergasflowFunctions()
