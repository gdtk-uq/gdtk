/**
 * Author: Rowan G.
 * Date: 2018-08-12
 */

module kinetics.luapseudo_species_kinetics;

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
import kinetics.pseudo_species_kinetics;

immutable string PseudoSpeciesKineticsMT = "PseudoSpeciesKinetics";

static const(PseudoSpeciesKinetics)[] PseudoSpeciesKineticsStore;

PseudoSpeciesKinetics checkPseudoSpeciesKinetics(lua_State *L, int index)
{
    return checkObj!(PseudoSpeciesKinetics, PseudoSpeciesKineticsMT)(L, index);
}

extern(C) int newPseudoSpeciesKinetics(lua_State* L)
{
    lua_remove(L, 1); // Remove first argument 'this'

    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
        string errMsg = "Error in call to PseudoSpeciesKinetics:new{}. " ~
            "A table containing named arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }

    lua_getfield(L, 1, "gasModel");
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to PseudoSpeciesKinetics:new{}. " ~
            "No gasmodel entry found in named arguments.";
        luaL_error(L, errMsg.toStringz());
    }
    auto gmodel = checkGasModel(L, -1);
    if (gmodel is null) {
        string errMsg = "Error in call to TwoTemperatureAirKinetics:new{}. " ~
            "A GasModel object was expected as the gasmodel argument. " ~
            "No valid GasModel was found.";
        luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1);

    auto myPseudoSpeciesKinetics = new PseudoSpeciesKinetics(gmodel);
    PseudoSpeciesKineticsStore ~= pushObj!(PseudoSpeciesKinetics, PseudoSpeciesKineticsMT)(L, myPseudoSpeciesKinetics);
    return 1;
}

extern(C) int updatePseudoSpeciesState(lua_State *L)
{
    // arg 1 is "self"
    auto pseudoSpeciesKinetics = checkPseudoSpeciesKinetics(L, 1);
    GasModel gm = pseudoSpeciesKinetics._gmodel;
    // arg 2 is GasState
    auto Q = GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    // arg 3 is tInterval
    double tInterval = luaL_checknumber(L, 3);

    // Some dummy parameters to keep method signature happy.
    double dtSuggest;
    // We need a dummy array of empty extra params
    // for the function signature
    number[maxParams] params;

    try {
        pseudoSpeciesKinetics(Q, tInterval, dtSuggest, params);
    }
    catch (ThermochemicalReactorUpdateException e) {
        string errMsg = "Error in call to pseudo species kinetics update. " ~
            "Caught exception: " ~ to!string(e);
        luaL_error(L, errMsg.toStringz);
    }
    // Update gas table.
    setGasStateInTable(L, gm, 2, Q);

    // Return new suggested time step size.
    lua_pushnumber(L, dtSuggest);
    return 2;
}

void registerPseudoSpeciesKinetics(lua_State* L)
{
    luaL_newmetatable(L, PseudoSpeciesKineticsMT.toStringz);

    // metatable.__index = metatable
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    // Register methods for use.
    lua_pushcfunction(L, &newPseudoSpeciesKinetics);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &updatePseudoSpeciesState);
    lua_setfield(L, -2, "updateState");
    lua_pushcfunction(L, &updatePseudoSpeciesState);
    lua_setfield(L, -2, "__call");

    // Make class visible.
    lua_setglobal(L, PseudoSpeciesKineticsMT.toStringz);
}
