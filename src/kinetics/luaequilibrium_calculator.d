/**
 * Author: Nick N. Gibbons
 * Date: 2021-06-08
 */

module kinetics.luaequilibrium_calculator;

import std.stdio;
import std.conv;
import std.string;
import ntypes.complex;
import nm.number;
import util.lua;
import util.lua_service;
import gas;
import gas.luagas_model;

import kinetics.equilibrium_update;


// name for EquilibriumCalculator in Lua scripts
immutable string EquilibriumCalculatorMT = "EquilibriumCalculator";

// Since we have two garbage collectors at play
// in D and Lua, it simplifies things to hang
// onto a store of objects in D's memory space.
static const(EquilibriumCalculator)[] EquilibriumCalculatorStore;

EquilibriumCalculator checkEquilibriumCalculator(lua_State* L, int index)
{
    return checkObj!(EquilibriumCalculator, EquilibriumCalculatorMT)(L, index);
}

/**
 * This function implements the constructor for a EquilibriumCalculator
 * from the Lua interface.
 *
 * Construction of a EquilibriumCalculator requires a thermally perfect
 * gasmodel file.
 *
 * ----------------------------------------------------------------------
 * EquilibriumCalculator:new{filename='gmodelname'}
 * ----------------------------------------------------------------------
 */
extern(C) int newEquilibriumCalculator(lua_State* L)
{
    lua_remove(L, 1); // Remove first argument 'this'

    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
        string errMsg = "Error in call to EquilibriumCalculator:new{}. " ~
            "A table containing named arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect to find a 'filename' entry
    lua_getfield(L, 1, "filename");
    if ( !lua_isstring(L, -1) ) {
        string errMsg = "Error in call to EquilibriumCalculator:new{}. " ~
            "A string was expected as the filename argument. " ~
            "No valid string was found.";
        luaL_error(L, errMsg.toStringz);
    }
    auto fname = to!string(luaL_checkstring(L, -1));
    lua_pop(L, 1);

    auto myEquilibriumCalculator = new EquilibriumCalculator(fname);
    EquilibriumCalculatorStore ~= pushObj!(EquilibriumCalculator, EquilibriumCalculatorMT)(L, myEquilibriumCalculator);
    return 1;
} // end newEquilibriumCalculator

// ----------------------------------------------------
// Exposed methods of the EquilibriumCalculator class
// ----------------------------------------------------
extern(C) int set_massf_from_pT(lua_State* L)
{
    // Arg 1 is "self"
    auto eqCalc = checkEquilibriumCalculator(L, 1);
    // Arg 3 is a GasModel
    auto gm = checkGasModel(L, 3);
    // Arg 2 is GasState
    auto Q = GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);

    try {
        eqCalc.set_massf_from_pT(Q);
    }
    catch (GasModelException e) {
        string errMsg = "Error in call to equilibrium calculator. " ~
            "Caught exception: " ~ to!string(e);
        luaL_error(L, errMsg.toStringz);
    }
    // Update gas table
    setGasStateInTable(L, gm, 2, Q);
    return 0;
}
extern(C) int set_massf_and_T_from_rhou(lua_State* L)
{
    // Arg 1 is "self"
    auto eqCalc = checkEquilibriumCalculator(L, 1);
    // Arg 3 is a GasModel
    auto gm = checkGasModel(L, 3);
    // Arg 2 is GasState
    auto Q = GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);

    try {
        eqCalc.set_massf_and_T_from_rhou(Q);
    }
    catch (GasModelException e) {
        string errMsg = "Error in call to equilibrium calculator. " ~
            "Caught exception: " ~ to!string(e);
        luaL_error(L, errMsg.toStringz);
    }
    // Update gas table
    setGasStateInTable(L, gm, 2, Q);
    return 0;
}
extern(C) int set_massf_and_T_from_ps(lua_State* L)
{
    // Arg 1 is "self"
    auto eqCalc = checkEquilibriumCalculator(L, 1);
    // Arg 4 is a GasModel
    auto gm = checkGasModel(L, 4);
    // Arg 2 is GasState
    auto Q = GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    // Arg 3 is the entropy
    double s = luaL_checknumber(L, 3);

    try {
        eqCalc.set_massf_and_T_from_ps(Q,s);
    }
    catch (GasModelException e) {
        string errMsg = "Error in call to equilibrium calculator. " ~
            "Caught exception: " ~ to!string(e);
        luaL_error(L, errMsg.toStringz);
    }
    // Update gas table
    setGasStateInTable(L, gm, 2, Q);
    return 0;
}
extern(C) int set_massf_from_rhoT(lua_State* L)
{
    // Arg 1 is "self"
    auto eqCalc = checkEquilibriumCalculator(L, 1);
    // Arg 3 is a GasModel
    auto gm = checkGasModel(L, 3);
    // Arg 2 is GasState
    auto Q = GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);

    try {
        eqCalc.set_massf_from_rhoT(Q);
    }
    catch (GasModelException e) {
        string errMsg = "Error in call to equilibrium calculator. " ~
            "Caught exception: " ~ to!string(e);
        luaL_error(L, errMsg.toStringz);
    }
    // Update gas table
    setGasStateInTable(L, gm, 2, Q);
    return 0;
}
extern(C) int get_entropy(lua_State* L)
{
    // Arg 1 is "self"
    auto eqCalc = checkEquilibriumCalculator(L, 1);
    // Arg 3 is a GasModel
    auto gm = checkGasModel(L, 3);
    // Arg 2 is GasState
    auto Q = GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);

    double s = eqCalc.get_s(Q);
    // Return entropy
    lua_pushnumber(L, s);
    return 1;
}

// --------- end: exposed methods ----------------- //

void registerEquilibriumCalculator(lua_State* L)
{
    luaL_newmetatable(L, EquilibriumCalculatorMT.toStringz);

    // metatable.__index = metatable
    lua_pushvalue(L, -1); // duplicate current metatable
    lua_setfield(L, -2, "__index");
    // Register methods for use
    lua_pushcfunction(L, &newEquilibriumCalculator);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &set_massf_from_pT);
    lua_setfield(L, -2, "set_massf_from_pT");
    lua_pushcfunction(L, &set_massf_and_T_from_rhou);
    lua_setfield(L, -2, "set_massf_and_T_from_rhou");
    lua_pushcfunction(L, &set_massf_and_T_from_ps);
    lua_setfield(L, -2, "set_massf_and_T_from_ps");
    lua_pushcfunction(L, &set_massf_from_rhoT);
    lua_setfield(L, -2, "set_massf_from_rhoT");
    lua_pushcfunction(L, &get_entropy);
    lua_setfield(L, -2, "get_entropy");

    // Make class visisble
    lua_setglobal(L, EquilibriumCalculatorMT.toStringz);
}
