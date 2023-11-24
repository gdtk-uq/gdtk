/**
 * Author: Rowan G.
 * Date: 2018-01-02
 */

module kinetics.luatwo_temperature_air_kinetics;

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
import kinetics.two_temperature_air_kinetics;

immutable string TwoTemperatureAirKineticsMT = "TwoTemperatureAirKinetics";

static const(TwoTemperatureAirKinetics)[] TwoTemperatureAirKineticsStore;

TwoTemperatureAirKinetics checkTwoTemperatureAirKinetics(lua_State *L, int index)
{
    return checkObj!(TwoTemperatureAirKinetics, TwoTemperatureAirKineticsMT)(L, index);
}

extern(C) int newTwoTemperatureAirKinetics(lua_State* L)
{
    lua_remove(L, 1); // Remove first argument 'this'

    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
        string errMsg = "Error in call to TwoTemperatureAirKinetics:new{}. " ~
            "A table containing named arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect to find a 'chemFilen' entry
    lua_getfield(L, 1, "chemFile");
    if ( !lua_isstring(L, -1) ) {
        string errMsg = "Error in call to TwoTemperatureAirKinetics:new{}. " ~
            "A string was expected as the chemFile argument. " ~
            "No valid string was found.";
        luaL_error(L, errMsg.toStringz);
    }
    auto chemFile = to!string(luaL_checkstring(L, -1));
    lua_pop(L, 1);
    // Expect to find a 'energyExchFile' entry
    lua_getfield(L, 1, "energyExchFile");
    if ( !lua_isstring(L, -1) ) {
        string errMsg = "Error in call to TwoTemperatureAirKinetics:new{}. " ~
            "A string was expected as the energyExchFile argument. " ~
            "No valid string was found.";
        luaL_error(L, errMsg.toStringz);
    }
    auto energyExchFile = to!string(luaL_checkstring(L, -1));
    lua_pop(L, 1);
    // Expect to find a 'gasModel' entry
    lua_getfield(L, 1, "gasModel");
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to TwoTemperatureAirKinetics:new{}. " ~
            "No gasmodel entry found in named arguments.";
        luaL_error(L, errMsg.toStringz());
    }
    auto gmodel = checkGasModel(L, -1);
    if ( gmodel is null ) {
        string errMsg = "Error in call to TwoTemperatureAirKinetics:new{}. " ~
            "A GasModel object was expected as the gasmodel argument. " ~
            "No valid GasModel was found.";
        luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1);

    auto myTwoTemperatureAirKinetics = new TwoTemperatureAirKinetics(chemFile, energyExchFile, gmodel);
    TwoTemperatureAirKineticsStore ~= pushObj!(TwoTemperatureAirKinetics, TwoTemperatureAirKineticsMT)(L, myTwoTemperatureAirKinetics);
    return 1;
}

extern(C) int updateTwoTempAirState(lua_State *L)
{
    // arg 1 is "self"
    auto twoTempAirKinetics = checkTwoTemperatureAirKinetics(L, 1);
    GasModel gm = twoTempAirKinetics._gmodel;
    // arg 2 is GasState
    auto Q = GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    // arg 3 is tInterval
    double tInterval = luaL_checknumber(L, 3);
    // arg 4 is dtSuggest
    double dtSuggest = luaL_checknumber(L, 4);
    //
    // We need a dummy array of empty extra params
    // for the function signature
    number[maxParams] params;

    try {
        twoTempAirKinetics(Q, tInterval, dtSuggest, params);
    }
    catch (ThermochemicalReactorUpdateException e) {
        string errMsg = "Error in call to two temperature air kinetics update. " ~
            "Caught exception: " ~ to!string(e);
        luaL_error(L, errMsg.toStringz);
    }
    // Update gas table.
    setGasStateInTable(L, gm, 2, Q);

    // Return newly suggested dt.
    lua_pushnumber(L, dtSuggest);
    return 2;
}

void registerTwoTemperatureAirKinetics(lua_State* L)
{
    luaL_newmetatable(L, TwoTemperatureAirKineticsMT.toStringz);

    // metatable.__index = metatable
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    // Register methods for use.
    lua_pushcfunction(L, &newTwoTemperatureAirKinetics);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &updateTwoTempAirState);
    lua_setfield(L, -2, "updateState");
    lua_pushcfunction(L, &updateTwoTempAirState);
    lua_setfield(L, -2, "__call");

    // Make class visible.
    lua_setglobal(L, TwoTemperatureAirKineticsMT.toStringz);
}
