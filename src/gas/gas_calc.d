// gas_calc.d
// Console program for running gas calculations via the Lua functions.
//
// Author: Rowan G.
// 2021-01-30: Extracted from luagas_model.d
//
module gas_calc;

import std.stdio;
import std.conv;
import std.string;
import util.lua;
import util.lua_service;

import gas.luagas_model;
import kinetics.luathermochemical_reactor;
import kinetics.luareaction_mechanism;
import kinetics.luachemistry_update;
import kinetics.luaequilibrium_calculator;
import kinetics.luatwo_temperature_air_kinetics;
import kinetics.luavib_specific_nitrogen_kinetics;
version (with_dvode)
{
    import kinetics.luapseudo_species_kinetics;
}
import luaidealgasflow;
import luagasflow;
import nm.luabbla;

int main(string[] args) {
    if ( args.length < 2 ) {
        writeln("ERROR: Wrong number of arguments.");
        writeln("");
        writeln("Usage: ");
        writeln("  gas-calc inputFile (+ script args)");
        return 1;
    }
    auto L = luaL_newstate();
    luaL_openlibs(L);
    lua_pushglobaltable(L);
    registerGasModel(L, -1);
    registerThermochemicalReactor(L, -1);
    registerReactionMechanism(L, -1);
    registerChemistryUpdate(L, -1);
    registerEquilibriumCalculator(L, -1);
    registerTwoTemperatureAirKinetics(L, -1);
    registerVibSpecNitrogenKinetics(L, -1);
    version(with_dvode) {
        registerPseudoSpeciesKinetics(L, -1);
    }
    registeridealgasflowFunctions(L);
    registergasflowFunctions(L);
    registerBBLA(L);
    // Pass on command line args to user's scripts.
    lua_newtable(L);
    int argIdx = 1;
    foreach (i; 2 .. args.length) {
        lua_pushstring(L, args[i].toStringz);
        lua_rawseti(L, -2, argIdx);
        argIdx++;
    }
    lua_setglobal(L, "arg");

    if ( luaL_dofile(L, toStringz(args[1])) != 0 ) {
        writeln(to!string(lua_tostring(L, -1)));
        return 1;
    }
    return 0;
}
