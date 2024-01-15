/**
 * Author: Rowan G.
 * Date: 2016-02-20
 */

module kinetics.luareaction_mechanism;

import std.stdio;
import std.conv;
import std.string;
import ntypes.complex;
import nm.number;
import util.lua;
import util.lua_service;
import gas;
import gas.luagas_model;

import kinetics.reaction_mechanism;

// name for ReactionMechanism in Lua scripts
immutable string ReactionMechanismMT = "ReactionMechanism";

// Since we have two garbage collectors at play
// in D and Lua, it simplifies things to hang
// onto a store of objects in D's memory space.
static const(ReactionMechanism)[] ReactionMechanismStore;

ReactionMechanism checkReactionMechanism(lua_State* L, int index)
{
    return checkObj!(ReactionMechanism, ReactionMechanismMT)(L, index);
}

/**
 * This function implements the constructor for a ReactionMechanism
 * from the Lua interface.
 *
 * Construction of a ReactionMechanism is from a filename,
 * a previously-constructed GasModel, and, optionally,
 * values for T_lower and T_upper for limits for evaluating the reaction
 * rate constants.
 * --------------------------------------------------------------------
 * rmech = ReactionMechanism:new{filename='fname', gasmodel=gmodel,
 *                               T_lower=300.0, T_upper=30000.0}
 * ---------------------------------------------------------------------
 */
extern(C) int newReactionMechanism(lua_State* L)
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
    // Optionally look for T_lower and T_upper
    double T_lower;
    lua_getfield(L, 1, "T_lower");
    if ( lua_isnumber(L, -1) ) {
        T_lower = luaL_checknumber(L, -1);
    }
    else {
        T_lower = 300.0;
    }
    lua_pop(L, 1);
    double T_upper;
    lua_getfield(L, 1, "T_upper");
    if ( lua_isnumber(L, -1) ) {
        T_upper = luaL_checknumber(L, -1);
    }
    else {
        T_upper = 30000.0;
    }
    lua_pop(L, 1);

    auto L2 = init_lua_State();
    doLuaFile(L2, fname);
    lua_getglobal(L2, "reaction");
    auto myReacMech = createReactionMechanism(L2, gmodel, T_lower, T_upper);
    lua_close(L2);
    ReactionMechanismStore ~= pushObj!(ReactionMechanism, ReactionMechanismMT)(L, myReacMech);
    return 1;
} // end newReactionMechanism()

// ----------------------------------------------------
// Exposed methods of the ReactionMechanism class
// ----------------------------------------------------
extern(C) int nReactions(lua_State* L)
{
    auto rmech = checkReactionMechanism(L, 1);
    lua_pushinteger(L, rmech.n_reactions);
    return 1;
}

extern(C) int evalRateConstants(lua_State* L)
{
    auto rmech = checkReactionMechanism(L, 1);
    auto gm = checkGasModel(L, 2);
    auto Q = GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 3, Q);
    rmech.eval_rate_constants(gm, Q);
    return 0;
}

extern(C) int evalRates(lua_State* L)
{
    auto rmech = checkReactionMechanism(L, 1);
    auto gm = checkGasModel(L, 2);
    auto Q = GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 3, Q);
    number[] conc;
    conc.length = gm.n_species;
    gm.massf2conc(Q, conc);
    number[] rates;
    rates.length = gm.n_species;
    rmech.eval_rates(conc, rates);
    lua_newtable(L);
    foreach (int isp; 1 .. gm.n_species+1) {
        lua_pushnumber(L, rates[isp-1]);
        lua_rawseti(L, -2, isp);
    }
    return 1;
}

extern(C) int k_f(lua_State* L)
{
    auto rmech = checkReactionMechanism(L, 1);
    int ir = luaL_checkint(L, 2);
    lua_pushnumber(L, rmech.k_f(ir));
    return 1;
}

extern(C) int k_b(lua_State* L)
{
    auto rmech = checkReactionMechanism(L, 1);
    int ir = luaL_checkint(L, 2);
    lua_pushnumber(L, rmech.k_b(ir));
    return 1;
}

extern(C) int rate(lua_State* L)
{
    auto rmech = checkReactionMechanism(L, 1);
    int ir = luaL_checkint(L, 2);
    int isp = luaL_checkint(L, 3);
    lua_pushnumber(L, rmech.rate(ir, isp));
    return 1;
}

// --------- end: exposed methods ----------------- //

void registerReactionMechanism(lua_State* L)
{
    luaL_newmetatable(L, ReactionMechanismMT.toStringz);

    // metatable.__index = metatable
    lua_pushvalue(L, -1); // duplicates current metatable
    lua_setfield(L, -2, "__index");
    // Register methods for use
    lua_pushcfunction(L, &newReactionMechanism);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &nReactions);
    lua_setfield(L, -2, "nReactions");
    lua_pushcfunction(L, &evalRateConstants);
    lua_setfield(L, -2, "evalRateConstants");
    lua_pushcfunction(L, &evalRates);
    lua_setfield(L, -2, "evalRates");
    lua_pushcfunction(L, &k_f);
    lua_setfield(L, -2, "k_f");
    lua_pushcfunction(L, &k_b);
    lua_setfield(L, -2, "k_b");
    lua_pushcfunction(L, &rate);
    lua_setfield(L, -2, "rate");

    // Make class visible
    lua_setglobal(L, ReactionMechanismMT.toStringz);
}
