/**
 * Author: Rowan G.
 * Date: 2016-02-22
 */

module kinetics.luachemistry_update;

import std.stdio;
import std.conv;
import std.string;
import util.lua;
import util.lua_service;
import gas;
import gas.luagas_model;

import kinetics.chemistry_update;

// name for ReactionUpdateScheme in Lua scripts
immutable string ReactionUpdateSchemeMT = "ReactionUpdateScheme";

// Since we have two garbage collectors at play
// in D and Lua, it simplifies things to hang
// onto a store of objects in D's memory space.
static const(ReactionUpdateScheme)[] ReactionUpdateSchemeStore;

ReactionUpdateScheme checkReactionUpdateScheme(lua_State* L, int index)
{
    return checkObj!(ReactionUpdateScheme, ReactionUpdateSchemeMT)(L, index);
}

/**
 * This function implements the constructor for a ReactionUpdateScheme
 * from the Lua interface.
 *
 * Construction of a ReactionUpdateScheme is from a filename and
 * a previously-constructed GasModel.
 * ----------------------------------------------------------------------
 * rupdate = ReactionUpdateScheme:new{filename='fname', gasmodel=gmodel}
 * ----------------------------------------------------------------------
 */
extern(C) int newReactionUpdateScheme(lua_State* L)
{
    lua_remove(L, 1); // Remove first argument 'this'

    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = "Error in call to ReactionMechanism:new{}. " ~
	    "A table containing named arguments is expected, but no table was found.";
	luaL_error(L, errMsg.toStringz);
    }
    // Expect to find a 'filename' entry
    lua_getfield(L, 1, "filename");
    if ( !lua_isstring(L, -1) ) {
	string errMsg = "Error in call to ReactionMechanism:new{}. " ~
	    "A string was expected as the filename argument. " ~
	    "No valid string was found.";
	luaL_error(L, errMsg.toStringz);
    }
    auto fname = to!string(luaL_checkstring(L, -1));
    lua_pop(L, 1);
    // Expect to find a 'gasmodel' entry
    lua_getfield(L, 1, "gasmodel");
    if ( lua_isnil(L, -1) ) {
	string errMsg = "Error in call to ReactionMechanism:new{}. " ~
	    "No gasmodel entry found in named arguments.";
	luaL_error(L, errMsg.toStringz());
    }
    auto gmodel = checkGasModel(L, -1);
    if ( gmodel is null ) {
	string errMsg = "Error in call to ReactionMechanisms:new{}. " ~
	    "A GasModel object was expected as the gasmodel argument. " ~
	    "No valid GasModel was found.";
	luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1);
    
    auto myReacUpdate = new ReactionUpdateScheme(fname, gmodel);
    ReactionUpdateSchemeStore ~= pushObj!(ReactionUpdateScheme, ReactionUpdateSchemeMT)(L, myReacUpdate);
    return 1;
} // end newReactionUpdateScheme

// ----------------------------------------------------
// Exposed methods of the ReactionMechanism class
// ----------------------------------------------------
extern(C) int updateState(lua_State* L)
{
    // Arg 1 is "self"
    auto rupdate = checkReactionUpdateScheme(L, 1);
    // Arg 5 is gasmodel (grab this first for help with GasState)
    auto gm = checkGasModel(L, 5);
    // Arg 2 is GasState
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    // Arg 3 is tInterval
    double tInterval = luaL_checknumber(L, 3);
    // Arg 4 is dtSuggest
    double dtSuggest = luaL_checknumber(L, 4);

    try {
	rupdate.update_state(Q, tInterval, dtSuggest, gm);
    }
    catch (Exception e) {
	string errMsg = "Error in call to updateState(). " ~
	    "Caught exception: " ~ to!string(e);
	luaL_error(L, errMsg.toStringz);
    }
    
    // Return new dtSuggest
    lua_pushnumber(L, dtSuggest);
    return 1;
}

// --------- end: exposed methods ----------------- //

void registerReactionUpdateScheme(lua_State* L, int tblIdx)
{
    luaL_newmetatable(L, ReactionUpdateSchemeMT.toStringz);

    // metatable.__index = metatable
    lua_pushvalue(L, -1); // duplicate current metatable
    lua_setfield(L, -2, "__index");
    // Register methods for use
    lua_pushcfunction(L, &newReactionUpdateScheme);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &updateState);
    lua_setfield(L, -2, "updateState");

    // Make class visisble
    lua_setfield(L, tblIdx, ReactionUpdateSchemeMT.toStringz);
}
