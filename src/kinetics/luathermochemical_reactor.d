/**
 * luathermochemical_reactor.d
 * Lua wrapper for the base class.
 *
 * Author: Rowan G. and Peter J.
 * Date: 2018-10-27 Adapted from luathermochemistry_update.d
 */

module kinetics.luathermochemical_reactor;

import std.stdio;
import std.conv;
import std.string;
import ntypes.complex;
import nm.number;
import util.lua;
import util.lua_service;
import gas;
import gas.luagas_model;

import kinetics.thermochemical_reactor;
import kinetics.init_thermochemical_reactor;
// import kinetics.chemistry_update;

// name for ThermochemicalReactor in Lua scripts
immutable string ThermochemicalReactorMT = "ThermochemicalReactor";

// Since we have two garbage collectors at play
// in D and Lua, it simplifies things to hang
// onto a store of objects in D's memory space.
static const(ThermochemicalReactor)[] ThermochemicalReactorStore;

ThermochemicalReactor checkThermochemicalReactor(lua_State* L, int index)
{
    return checkObj!(ThermochemicalReactor, ThermochemicalReactorMT)(L, index);
}

/**
 * This function implements the constructor for a ThermochemicalReactor
 * from the Lua interface.
 *
 * Construction of a ThermochemicalReactor is from a previously-constructed
 * GasModel and, maybe, one or two other Lua files.
 * The number of files needed and their content depends on the specific gas model.
 * ----------------------------------------------------------------------
 * rupdate = ThermochemicalReactor:new{gasmodel=gmodel,
 *                                     filename1='fname1', filename2='fname2'}
 * ----------------------------------------------------------------------
 */
extern(C) int newThermochemicalReactor(lua_State* L)
{
    lua_remove(L, 1); // Remove first argument 'this'

    int narg = lua_gettop(L);
    if (narg == 0 || !lua_istable(L, 1)) {
        string errMsg = "Error in call to ThermochemicalReactor:new{}. " ~
            "A table containing named arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect to find a 'gasmodel' entry
    lua_getfield(L, 1, "gasmodel");
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to ThermochemicalReactor:new{}. " ~
            "No gasmodel entry found in named arguments.";
        luaL_error(L, errMsg.toStringz());
    }
    auto gmodel = checkGasModel(L, -1);
    if ( gmodel is null ) {
        string errMsg = "Error in call to ThermochemicalReactor:new{}. " ~
            "A GasModel object was expected as the gasmodel argument. " ~
            "No valid GasModel was found.";
        luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1);
    // Look for an optional 'filename1' entry.
    string fileName1 = "";
    lua_getfield(L, 1, "filename1");
    if (lua_isstring(L, -1)) { fileName1 = to!string(luaL_checkstring(L, -1)); }
    lua_pop(L, 1);
    // Look for an optional 'filename2' entry.
    string fileName2 = "";
    lua_getfield(L, 1, "filename2");
    if (lua_isstring(L, -1)) { fileName2 = to!string(luaL_checkstring(L, -1)); }
    lua_pop(L, 1);

    auto myReactor = init_thermochemical_reactor(gmodel, fileName1, fileName2);
    ThermochemicalReactorStore ~= pushObj!(ThermochemicalReactor, ThermochemicalReactorMT)(L, myReactor);
    return 1;
} // end newThermochemicalReactor

// -------------------------------------------------------
// Advance the state of the thermochemical system in time.
// -------------------------------------------------------
extern(C) int updateThermochemicalState(lua_State* L)
{
    // Arg 1 is "self"
    auto myReactor = checkThermochemicalReactor(L, 1);
    GasModel gm = myReactor._gmodel;
    // Arg 2 is GasState
    auto Q = GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    // Arg 3 is tInterval
    double tInterval = luaL_checknumber(L, 3);
    // Arg 4 is dtSuggest
    double dtSuggest = luaL_checknumber(L, 4);
    // Extra parameters are not considered, presently. PJ 2017-04-22
    // Need to use these for JJ's fuel-air mix kinetics scheme.
    number[maxParams] params;

    try {
        myReactor(Q, tInterval, dtSuggest, params);
    }
    catch (ThermochemicalReactorUpdateException e) {
        string errMsg = "Error in call to thermochemical state advance. " ~
            "Caught exception: " ~ to!string(e);
        luaL_error(L, errMsg.toStringz);
    }
    // Update gas table
    setGasStateInTable(L, gm, 2, Q);

    // Return newly suggested dt.
    lua_pushnumber(L, dtSuggest);
    return 1;
}

// --------- end: exposed methods ----------------- //

void registerThermochemicalReactor(lua_State* L)
{
    luaL_newmetatable(L, ThermochemicalReactorMT.toStringz);

    // metatable.__index = metatable
    lua_pushvalue(L, -1); // duplicate current metatable
    lua_setfield(L, -2, "__index");
    // Register methods for use
    lua_pushcfunction(L, &newThermochemicalReactor);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &updateThermochemicalState);
    lua_setfield(L, -2, "updateState");
    lua_pushcfunction(L, &updateThermochemicalState);
    lua_setfield(L, -2, "__call");

    // Make class visisble
    lua_setglobal(L, ThermochemicalReactorMT.toStringz);
}
