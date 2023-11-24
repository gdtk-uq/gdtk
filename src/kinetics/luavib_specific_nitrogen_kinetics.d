/**
 * Author: Rowan G.
 * Date: 2018-05-08
 */

module kinetics.luavib_specific_nitrogen_kinetics;

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
import kinetics.vib_specific_nitrogen_kinetics;

immutable string VibSpecNitrogenKineticsMT = "VibSpecNitrogenKinetics";

static const(VibSpecificNitrogenRelaxation)[] VibSpecNitrogenKineticsStore;

VibSpecificNitrogenRelaxation checkVibSpecNitrogenKinetics(lua_State *L, int index)
{
    return checkObj!(VibSpecificNitrogenRelaxation, VibSpecNitrogenKineticsMT)(L, index);
}

extern(C) int newVibSpecNitrogenKinetics(lua_State* L)
{
    lua_remove(L, 1); // Remove first argument 'this'

    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
        string errMsg = "Error in call to VibSpecNitrogenKinetics:new{}. " ~
            "A table containing named arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect to find a 'gasModel' entry
    lua_getfield(L, 1, "gasModel");
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to VibSpecNitrogenKinetics:new{}. " ~
            "No gasModel entry found in named arguments.";
        luaL_error(L, errMsg.toStringz());
    }
    auto gmodel = checkGasModel(L, -1);
    if ( gmodel is null ) {
        string errMsg = "Error in call to VibSpecNitrogenKinetics:new{}. " ~
            "A GasModel object was expected as the gasmodel argument. " ~
            "No valid GasModel was found.";
        luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1);

    auto myVibSpecNitrogenKinetics = new VibSpecificNitrogenRelaxation("", gmodel);
    VibSpecNitrogenKineticsStore ~= pushObj!(VibSpecificNitrogenRelaxation, VibSpecNitrogenKineticsMT)(L, myVibSpecNitrogenKinetics);
    return 1;
}

extern(C) int updateNitrogenStates(lua_State *L)
{
    // arg 1 is "self"
    auto vibSpecN2Kinetics = checkVibSpecNitrogenKinetics(L, 1);
    GasModel gm = vibSpecN2Kinetics._gmodel;
    // arg 2 is GasState
    auto Q = GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    // arg 3 is tInterval
    double tInterval = luaL_checknumber(L, 3);

    // We need a dummy array of empty extra params
    // for the function signature
    number[maxParams] params;
    // and dummy dtSuggest
    double dtSuggest = -1.0;

    try {
        vibSpecN2Kinetics(Q, tInterval, dtSuggest, params);
    }
    catch (ThermochemicalReactorUpdateException e) {
        string errMsg = "Error in call to two temperature air kinetics update. " ~
            "Caught exception: " ~ to!string(e);
        luaL_error(L, errMsg.toStringz);
    }
    // Update gas table.
    setGasStateInTable(L, gm, 2, Q);

    return 0;
}

void registerVibSpecNitrogenKinetics(lua_State* L)
{
    luaL_newmetatable(L, VibSpecNitrogenKineticsMT.toStringz);

    // metatable.__index = metatable
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    // Register methods for use.
    lua_pushcfunction(L, &newVibSpecNitrogenKinetics);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &updateNitrogenStates);
    lua_setfield(L, -2, "updateState");
    lua_pushcfunction(L, &updateNitrogenStates);
    lua_setfield(L, -2, "__call");

    // Make class visible.
    lua_setglobal(L, VibSpecNitrogenKineticsMT.toStringz);
}
