/**
 * Author: Rowan G.
 * Date: 2016-02-22
 */

module kinetics.luachemistry_update;

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
import kinetics.chemistry_update;


// name for ChemistryUpdate in Lua scripts
immutable string ChemistryUpdateMT = "ChemistryUpdate";

// Since we have two garbage collectors at play
// in D and Lua, it simplifies things to hang
// onto a store of objects in D's memory space.
static const(ChemistryUpdate)[] ChemistryUpdateStore;

ChemistryUpdate checkChemistryUpdate(lua_State* L, int index)
{
    return checkObj!(ChemistryUpdate, ChemistryUpdateMT)(L, index);
}

/**
 * This function implements the constructor for a ChemistryUpdate
 * from the Lua interface.
 *
 * Construction of a ChemistryUpdate is from a filename and
 * a previously-constructed GasModel.
 * ----------------------------------------------------------------------
 * rupdate = ChemistryUpdate:new{filename='fname', gasmodel=gmodel}
 * ----------------------------------------------------------------------
 */
extern(C) int newChemistryUpdate(lua_State* L)
{
    lua_remove(L, 1); // Remove first argument 'this'

    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
        string errMsg = "Error in call to ChemistryUpdate:new{}. " ~
            "A table containing named arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect to find a 'filename' entry
    lua_getfield(L, 1, "filename");
    if ( !lua_isstring(L, -1) ) {
        string errMsg = "Error in call to ChemistryUpdate:new{}. " ~
            "A string was expected as the filename argument. " ~
            "No valid string was found.";
        luaL_error(L, errMsg.toStringz);
    }
    auto fname = to!string(luaL_checkstring(L, -1));
    lua_pop(L, 1);
    // Expect to find a 'gasmodel' entry
    lua_getfield(L, 1, "gasmodel");
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to ChemistryUpdate:new{}. " ~
            "No gasmodel entry found in named arguments.";
        luaL_error(L, errMsg.toStringz());
    }
    auto gmodel = checkGasModel(L, -1);
    if ( gmodel is null ) {
        string errMsg = "Error in call to ChemistryUpdate:new{}. " ~
            "A GasModel object was expected as the gasmodel argument. " ~
            "No valid GasModel was found.";
        luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1);

    auto myChemUpdate = new ChemistryUpdate(fname, gmodel);
    ChemistryUpdateStore ~= pushObj!(ChemistryUpdate, ChemistryUpdateMT)(L, myChemUpdate);
    return 1;
} // end newChemistryUpdate

// ----------------------------------------------------
// Exposed methods of the ReactionMechanism class
// ----------------------------------------------------
extern(C) int updateState(lua_State* L)
{
    // Arg 1 is "self"
    auto chemUpdate = checkChemistryUpdate(L, 1);
    // Arg 5 is gasmodel (grab this first for help with GasState)
    // [FIX-ME] Rowan, should be able to dispose of this argument from your scripts
    // now that a reference is packed inside the ChemistryUpdate object.
    // auto gm = checkGasModel(L, 5);
    GasModel gm = chemUpdate._gmodel;
    // Arg 2 is GasState
    auto Q = GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    // Arg 3 is tInterval
    double tInterval = luaL_checknumber(L, 3);
    // Arg 4 is dtSuggest
    double dtSuggest = luaL_checknumber(L, 4);
    // Extra parameters are not considered, presently. PJ 2017-04-22
    number[maxParams] params;

    try {
        chemUpdate(Q, tInterval, dtSuggest, params);
    }
    catch (ThermochemicalReactorUpdateException e) {
        string errMsg = "Error in call to chemistry update. " ~
            "Caught exception: " ~ to!string(e);
        luaL_error(L, errMsg.toStringz);
    }
    // Update gas table
    setGasStateInTable(L, gm, 2, Q);

    // Return new dtSuggest
    lua_pushnumber(L, dtSuggest);
    return 1;
}

// --------- end: exposed methods ----------------- //

void registerChemistryUpdate(lua_State* L)
{
    luaL_newmetatable(L, ChemistryUpdateMT.toStringz);

    // metatable.__index = metatable
    lua_pushvalue(L, -1); // duplicate current metatable
    lua_setfield(L, -2, "__index");
    // Register methods for use
    lua_pushcfunction(L, &newChemistryUpdate);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &updateState);
    lua_setfield(L, -2, "updateState");
    lua_pushcfunction(L, &updateState);
    lua_setfield(L, -2, "__call");

    // Make class visisble
    lua_setglobal(L, ChemistryUpdateMT.toStringz);
}
