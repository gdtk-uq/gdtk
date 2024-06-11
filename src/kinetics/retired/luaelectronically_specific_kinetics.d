/**
 * Author: Rowan G.
 * Date: 2018-10-17
 */

module kinetics.luaelectronically_specific_kinetics;

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
import kinetics.electronically_specific_kinetics;

immutable string ElectronicallySpecificKineticsMT = "ElectronicallySpecificKinetics";

static const(ElectronicallySpecificKinetics)[] ESKStore;

ElectronicallySpecificKinetics checkElectronicallySpecificKinetics(lua_State *L, int index)
{
    return checkObj!(ElectronicallySpecificKinetics, ElectronicallySpecificKineticsMT)(L, index);
}

extern(C) int newElectronicallySpecificKinetics(lua_State* L)
{
    lua_remove(L, 1); // Remove first argument 'this'

    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
        string errMsg = "Error in call to ElectronicallySpecificKinetics:new{}. " ~
            "A table containing named arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect to find a 'listOfFiles' entry
    lua_getfield(L, 1, "listOfFiles");
    if ( !lua_isstring(L, -1) ) {
        string errMsg = "Error in call to ElectronicallySpecificKinetics:new{}. " ~
            "A string was expected as the filename1 argument. " ~
            "No valid string was found.";
        luaL_error(L, errMsg.toStringz);
    }
    auto listOfFiles = to!string(luaL_checkstring(L, -1));
    lua_pop(L, 1);

    // Expect to find a 'gasModel' entry
    lua_getfield(L, 1, "gasModel");
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to ElectronicallySpecificKinetics:new{}. " ~
            "No gasModel entry found in named arguments.";
        luaL_error(L, errMsg.toStringz());
    }
    auto gmodel = checkGasModel(L, -1);
    if ( gmodel is null ) {
        string errMsg = "Error in call to ElectronicallySpecificKinetics:new{}. " ~
            "A GasModel object was expected as the gasmodel argument. " ~
            "No valid GasModel was found.";
        luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1);

    auto myESK = new ElectronicallySpecificKinetics(listOfFiles,gmodel);
    ESKStore ~= pushObj!(ElectronicallySpecificKinetics, ElectronicallySpecificKineticsMT)(L, myESK);
    return 1;
}

extern(C) int updateElectronicStates(lua_State *L)
{
    // arg 1 is "self"
    auto myESK = checkElectronicallySpecificKinetics(L, 1);
    GasModel gm = myESK._gmodel;
    // arg 2 is GasState
    auto Q = GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, gm, 2, Q);
    // arg 3 is tInterval
    double tInterval = luaL_checknumber(L, 3);

    // We need a dummy array of empty extra params
    // for the function signature
    number[maxParams] params;
    // and dtSuggest
    double dtSuggest = luaL_checknumber(L, 4);
    try {
        myESK(Q, tInterval, dtSuggest, params);
    }
    catch (ThermochemicalReactorUpdateException e) {
        string errMsg = "Error in call to electronically-specific kinetics update. " ~
            "Caught exception: " ~ to!string(e);
        luaL_error(L, errMsg.toStringz);
    }
    // Update gas table.
    setGasStateInTable(L, gm, 2, Q);
    // Return new dtSuggest
    lua_pushnumber(L, dtSuggest);
    return 0;
}

void registerElectronicallySpecificKinetics(lua_State* L)
{
    luaL_newmetatable(L, ElectronicallySpecificKineticsMT.toStringz);

    // metatable.__index = metatable
    lua_pushvalue(L, -1);
    lua_setfield(L, -2, "__index");
    // Register methods for use.
    lua_pushcfunction(L, &newElectronicallySpecificKinetics);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &updateElectronicStates);
    lua_setfield(L, -2, "update");
    lua_pushcfunction(L, &updateElectronicStates);
    lua_setfield(L, -2, "__call");

    // Make class visible.
    lua_setglobal(L, ElectronicallySpecificKineticsMT.toStringz);
}
