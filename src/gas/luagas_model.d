/**
 * luagas_model.d
 * Lua interface to the gas model.
 *
 * Author: Rowan G.
 * Version: Initial cut.
 */

module gas.luagas_model;

import std.math;
import std.algorithm;
import std.stdio;
import std.conv;
import std.string;
import util.lua;
import util.lua_service;
import ntypes.complex;
import nm.number;

import kinetics.luathermochemical_reactor;
import kinetics.luareaction_mechanism;
import kinetics.luachemistry_update;
import kinetics.luatwo_temperature_air_kinetics;
import kinetics.luavib_specific_nitrogen_kinetics;
version (with_dvode)
{
    import kinetics.luapseudo_species_kinetics;
}

import gas.gas_model;
import gas.gas_state;
import gas.init_gas_model: init_gas_model;
import gas.cea_gas;
import gas.physical_constants;

// name for GasModel class in Lua scripts
immutable string GasModelMT = "GasModel";

// Since we have two garbage collectors at play
// in D and Lua, it simplifies things to hang
// onto a store of objects in D's memory space.
static const(GasModel)[] GasModelStore;

GasModel checkGasModel(lua_State* L, int index)
{
    return checkObj!(GasModel, GasModelMT)(L, index);
}

/**
 * This function implements the constructor for the GasModel
 * from the Lua interface.
 *
 * Construction of GasModel object in Lua will accept a filename.
 * All of the configuration details needed to construct any
 * particular gas model are contained within that file.
 * -------------
 * gmodel = GasModel:new{fname}
 * -------------
 */
extern(C) int newGasModel(lua_State* L)
{
    lua_remove(L, 1); // Remove first argument this

    lua_rawgeti(L, 1, 1);
    string fname = to!string(luaL_checkstring(L, -1));
    lua_pop(L, 1);
    GasModel gm = init_gas_model(fname);
    GasModelStore ~= pushObj!(GasModel, GasModelMT)(L, gm);
    return 1;
}

extern(C) int nSpecies(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    lua_pushnumber(L, gm.n_species);
    return 1;
}

extern(C) int speciesIndex(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    string spName = to!string(luaL_checkstring(L, 2));
    lua_pushnumber(L, gm.species_index(spName));
    return 1;
}

extern(C) int nModes(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    lua_pushnumber(L, gm.n_modes);
    return 1;
}

extern(C) int molMasses(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    lua_newtable(L);
    foreach (isp, mMass; gm.mol_masses ) {
        lua_pushnumber(L, mMass);
        lua_setfield(L, -2, toStringz(gm.species_name(isp)));
    }
    return 1;
}

extern(C) int charge(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    lua_newtable(L);
    foreach (isp, mCharge; gm.charge ) {
        lua_pushnumber(L, mCharge);
        lua_setfield(L, -2, toStringz(gm.species_name(isp)));
    }
    return 1;
}

extern(C) int speciesName(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    int i = to!int(lua_tointeger(L, 2));
    lua_pushstring(L, toStringz(gm.species_name(i)));
    return 1;
}

extern(C) int thermoPU(lua_State* L)
{
    try {
        auto gm = checkGasModel(L, 1);
        auto Q = GasState(gm);
        getGasStateFromTable(L, gm, 2, Q);
        if ( Q.p <= 0.0 || isNaN(Q.p) ) {
            string errMsg = "ERROR: when calling 'updateThermoFromPU'\n";
            errMsg ~= "        The supplied pressure value is negative, 0 or has not been set.\n";
            errMsg ~= "        Check that the 'p' field is set with a valid value.\n";
            errMsg ~= "        The complete gas state is:\n";
            errMsg ~= Q.toString();
            errMsg ~= "\nBailing out\n";
            throw new Exception(errMsg);
        }
	if ( isNaN(Q.u) ) {
        string errMsg = "ERROR: when calling 'updateThermoFromPU'\n";
        errMsg ~= "        The energy value has not been set.\n";
        errMsg ~= "        Check that the 'u' field is set with a valid value.\n";
        errMsg ~= "        The complete gas state is:\n";
        errMsg ~= Q.toString();
        errMsg ~= "\nBailing out\n";
        throw new Error(errMsg);
	}
        gm.update_thermo_from_pu(Q);
        gm.update_sound_speed(Q);
        setGasStateInTable(L, gm, 2, Q);
    } catch(Exception e) {
        luaL_error(L, "%s", e.msg.toStringz);
    }
    return 0;
}


extern(C) int thermoPT(lua_State* L)
{
    try {
        auto gm = checkGasModel(L, 1);
        auto Q = GasState(gm);
        getGasStateFromTable(L, gm, 2, Q);
        if ( Q.p <= 0.0 || isNaN(Q.p) ) {
            string errMsg = "ERROR: when calling 'updateThermoFromPT'\n";
            errMsg ~= "        The supplied pressure value is negative, 0 or has not been set.\n";
            errMsg ~= "        Check that the 'p' field is set with a valid value.\n";
            errMsg ~= "        The complete gas state is:\n";
            errMsg ~= Q.toString();
            errMsg ~= "\nBailing out\n";
            throw new Exception(errMsg);
        }
        if ( Q.T <= 0.0 || isNaN(Q.T) ) {
            string errMsg = "ERROR: when calling 'updateThermoFromPT'\n";
            errMsg ~= "        Supplied temperature value is negative,  0 or has not been set.\n";
            errMsg ~= "        Check that the 'T' field is set with a valid value.\n";
            errMsg ~= "        The complete gas state is:\n";
            errMsg ~= Q.toString();
            errMsg ~= "\nBailing out\n";
            throw new Exception(errMsg);
        }
        gm.update_thermo_from_pT(Q);
        gm.update_sound_speed(Q);
        setGasStateInTable(L, gm, 2, Q);
    } catch(Exception e) {
        luaL_error(L, "%s", e.msg.toStringz);
    }
    return 0;
}

extern(C) int thermoRHOU(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 2, Q);
    if ( Q.rho <= 0.0 || isNaN(Q.rho) ) {
        string errMsg = "ERROR: when calling 'updateThermoFromRHOU'\n";
        errMsg ~= "        The supplied density value is negative, 0 or has not been set.\n";
        errMsg ~= "        Check that the 'rho' field is set with a valid value.\n";
        errMsg ~= "        The complete gas state is:\n";
        errMsg ~= Q.toString();
        errMsg ~= "\nBailing out\n";
        throw new Error(errMsg);
    }
    if ( isNaN(Q.u) ) {
        string errMsg = "ERROR: when calling 'updateThermoFromRHOU'\n";
        errMsg ~= "        The energy value has not been set.\n";
        errMsg ~= "        Check that the 'u' field is set with a valid value.\n";
        errMsg ~= "        The complete gas state is:\n";
        errMsg ~= Q.toString();
        errMsg ~= "\nBailing out\n";
        throw new Error(errMsg);
    }
    gm.update_thermo_from_rhou(Q);
    gm.update_sound_speed(Q);
    setGasStateInTable(L, gm, 2, Q);
    return 0;
}

extern(C) int thermoRHOT(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 2, Q);
    if ( Q.rho <= 0.0 || isNaN(Q.rho) ) {
        string errMsg = "ERROR: when calling 'updateThermoFromRHOT'\n";
        errMsg ~= "        The supplied density value is negative, 0 or has not been set.\n";
        errMsg ~= "        Check that the 'rho' field is set with a valid value.\n";
        errMsg ~= "        The complete gas state is:\n";
        errMsg ~= Q.toString();
        errMsg ~= "\nBailing out\n";
        throw new Error(errMsg);
    }
    if ( Q.T <= 0.0 || isNaN(Q.T) ) {
        string errMsg = "ERROR: when calling 'updateThermoFromRHOT'\n";
        errMsg ~= "        The supplied temperature value is negative, 0 or has not been set.\n";
        errMsg ~= "        Check that the 'T' field is set with a valid value.\n";
        errMsg ~= "        The complete gas state is:\n";
        errMsg ~= Q.toString();
        errMsg ~= "\nBailing out\n";
        throw new Error(errMsg);
    }
    gm.update_thermo_from_rhoT(Q);
    gm.update_sound_speed(Q);
    setGasStateInTable(L, gm, 2, Q);
    return 0;
}

extern(C) int thermoRHOP(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 2, Q);
    if ( Q.rho <= 0.0 || isNaN(Q.rho) ) {
        string errMsg = "ERROR: when calling 'updateThermoFromRHOP'\n";
        errMsg ~= "        The supplied density value is negative, 0 or has not been set.\n";
        errMsg ~= "        Check that the 'rho' field is set with a valid value.\n";
        errMsg ~= "        The complete gas state is:\n";
        errMsg ~= Q.toString();
        errMsg ~= "\nBailing out\n";
        throw new Error(errMsg);
    }
    if ( Q.p <= 0.0 || isNaN(Q.p) ) {
        string errMsg = "ERROR: when calling 'updateThermoFromRHOP'\n";
        errMsg ~= "        The supplied pressure value is negative, 0 or has not been set.\n";
        errMsg ~= "        Check that the 'p' field is set with a valid value.\n";
        errMsg ~= "        The complete gas state is:\n";
        errMsg ~= Q.toString();
        errMsg ~= "\nBailing out\n";
        throw new Error(errMsg);
    }
    gm.update_thermo_from_rhop(Q);
    gm.update_sound_speed(Q);
    setGasStateInTable(L, gm, 2, Q);
    return 0;
}

extern(C) int thermoPS(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 2, Q);
    if ( Q.p <= 0.0 || isNaN(Q.p) ) {
        string errMsg = "ERROR: when calling 'updateThermoFromPS'\n";
        errMsg ~= "        The supplied pressure value is negative, 0 or has not been set.\n";
        errMsg ~= "        Check that the 'p' field is set with a valid value.\n";
        errMsg ~= "        The complete gas state is:\n";
        errMsg ~= Q.toString();
        errMsg ~= "\nBailing out\n";
        throw new Error(errMsg);
    }
    number s = lua_tonumber(L, 3);
    gm.update_thermo_from_ps(Q, s);
    gm.update_sound_speed(Q);
    setGasStateInTable(L, gm, 2, Q);
    return 0;
}

extern(C) int thermoHS(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 2, Q);
    number h = lua_tonumber(L, 3);
    number s = lua_tonumber(L, 4);
    gm.update_thermo_from_hs(Q, h, s);
    gm.update_sound_speed(Q);
    setGasStateInTable(L, gm, 2, Q);
    return 0;
}

extern(C) int soundSpeed(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 2, Q);
    if ( Q.T <= 0.0 || isNaN(Q.T) ) {
        string errMsg = "ERROR: when calling 'updateSoundSpeed'\n";
        errMsg ~= "        The supplied temperature value is negative, 0 or has not been set.\n";
        errMsg ~= "        Check that the 'T' field is set with a valid value.\n";
        errMsg ~= "        The complete gas state is:\n";
        errMsg ~= Q.toString();
        errMsg ~= "\nBailing out\n";
        throw new Error(errMsg);
    }
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(Q); }
    gm.update_sound_speed(Q);
    setGasStateInTable(L, gm, 2, Q);
    return 0;
}

extern(C) int transCoeffs(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 2, Q);
    if ( Q.T <= 0.0 || isNaN(Q.T) ) {
        string errMsg = "ERROR: when calling 'updateTransCoeffs'\n";
        errMsg ~= "        The supplied temperature value is negative, 0 or has not been set.\n";
        errMsg ~= "        Check that the 'T' field is set with a valid value.\n";
        errMsg ~= "        The complete gas state is:\n";
        errMsg ~= Q.toString();
        errMsg ~= "\nBailing out\n";
        throw new Error(errMsg);
    }
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(Q); }
    gm.update_trans_coeffs(Q);
    setGasStateInTable(L, gm, 2, Q);
    return 0;
}

// We don't wrap dudT_const_v as we do Cv in its place.
// We don't wrap dhdT_const_p as we do Cp in its place.

extern(C) int dpdrhoConstT(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 2, Q);
    if ( Q.T <= 0.0 || isNaN(Q.T) ) {
        string errMsg = "ERROR: when calling 'dpdrhoConstT'\n";
        errMsg ~= "        The supplied temperature value is negative, 0 or has not been set.\n";
        errMsg ~= "        Check that the 'T' field is set with a valid value.\n";
        errMsg ~= "        The complete gas state is:\n";
        errMsg ~= Q.toString();
        errMsg ~= "\nBailing out\n";
        throw new Error(errMsg);
    }
    auto dpdrho = gm.dpdrho_const_T(Q);
    lua_pushnumber(L, dpdrho);
    return 1;
}

// We don't wrap gas_constant as we do R in its place.

extern(C) int intEnergy(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 2, Q);
    if ( Q.T <= 0.0 || isNaN(Q.T) ) {
        string errMsg = "ERROR: when calling 'intEnergy'\n";
        errMsg ~= "        The supplied temperature value is negative, 0 or has not been set.\n";
        errMsg ~= "        Check that the 'T' field is set with a valid value.\n";
        errMsg ~= "        The complete gas state is:\n";
        errMsg ~= Q.toString();
        errMsg ~= "\nBailing out\n";
        throw new Error(errMsg);
    }
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(Q); }
    int narg = lua_gettop(L);
    auto e = gm.internal_energy(Q);
    lua_pushnumber(L, e);
    return 1;
}

extern(C) int enthalpy(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 2, Q);
    if ( Q.T <= 0.0 || isNaN(Q.T) ) {
        string errMsg = "ERROR: when calling 'enthalpy'\n";
        errMsg ~= "        The supplied temperature value is negative, 0 or has not been set.\n";
        errMsg ~= "        Check that the 'T' field is set with a valid value.\n";
        errMsg ~= "        The complete gas state is:\n";
        errMsg ~= Q.toString();
        errMsg ~= "\nBailing out\n";
        throw new Error(errMsg);
    }
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(Q); }
    int narg = lua_gettop(L);
    number h;
    if ( narg >= 3 ) { // Call species-specific version
        int isp = luaL_checkint(L, 3);
        h = gm.enthalpy(Q, isp);
    } else { // Call total gas mix version
        h = gm.enthalpy(Q);
    }
    lua_pushnumber(L, h);
    return 1;
}

extern(C) int entropy(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 2, Q);
    if ( Q.T <= 0.0 || isNaN(Q.T) ) {
        string errMsg = "ERROR: when calling 'entropy'\n";
        errMsg ~= "        The supplied temperature value is negative, 0 or has not been set.\n";
        errMsg ~= "        Check that the 'T' field is set with a valid value.\n";
        errMsg ~= "        The complete gas state is:\n";
        errMsg ~= Q.toString();
        errMsg ~= "\nBailing out\n";
        throw new Error(errMsg);
    }
    if (cast(CEAGas) gm !is null) { gm.update_thermo_from_pT(Q); }
    int narg = lua_gettop(L);
    number s;
    if ( narg >= 3 ) { // Call species-specific version
        int isp = luaL_checkint(L, 3);
        s = gm.entropy(Q, isp);
    } else {
        s = gm.entropy(Q);
    }
    lua_pushnumber(L, s);
    return 1;
}

extern(C) int gibbsFreeEnergy(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 2, Q);
    if ( Q.T <= 0.0 || isNaN(Q.T) ) {
        string errMsg = "ERROR: when calling 'gibbsFreeEnergy'\n";
        errMsg ~= "        The supplied temperature value is negative, 0 or has not been set.\n";
        errMsg ~= "        Check that the 'T' field is set with a valid value.\n";
        errMsg ~= "        The complete gas state is:\n";
        errMsg ~= Q.toString();
        errMsg ~= "\nBailing out\n";
        throw new Error(errMsg);
    }
    int narg = lua_gettop(L);
    number G;
    int isp = luaL_checkint(L, 3);
    G = gm.gibbs_free_energy(Q, isp);
    lua_pushnumber(L, G);
    return 1;
}


extern(C) int Cv(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 2, Q);
    if ( Q.T <= 0.0 || isNaN(Q.T) ) {
        string errMsg = "ERROR: when calling 'Cv'\n";
        errMsg ~= "        The supplied temperature value is negative, 0 or has not been set.\n";
        errMsg ~= "        Check that the 'T' field is set with a valid value.\n";
        errMsg ~= "        The complete gas state is:\n";
        errMsg ~= Q.toString();
        errMsg ~= "\nBailing out\n";
        throw new Error(errMsg);
    }
    auto Cv = gm.Cv(Q);
    lua_pushnumber(L, Cv);
    return 1;
}

extern(C) int Cp(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 2, Q);
    if ( Q.T <= 0.0 || isNaN(Q.T) ) {
        string errMsg = "ERROR: when calling 'Cp'\n";
        errMsg ~= "        The supplied temperature value is negative, 0 or has not been set.\n";
        errMsg ~= "        Check that the 'T' field is set with a valid value.\n";
        errMsg ~= "        The complete gas state is:\n";
        errMsg ~= Q.toString();
        errMsg ~= "\nBailing out\n";
        throw new Error(errMsg);
    }
    auto Cp = gm.Cp(Q);
    lua_pushnumber(L, Cp);
    return 1;
}

extern(C) int R(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 2, Q);
    if ( Q.T <= 0.0 || isNaN(Q.T) ) {
        string errMsg = "ERROR: when calling 'R'\n";
        errMsg ~= "        The supplied temperature value is negative, 0 or has not been set.\n";
        errMsg ~= "        Check that the 'T' field is set with a valid value.\n";
        errMsg ~= "        The complete gas state is:\n";
        errMsg ~= Q.toString();
        errMsg ~= "\nBailing out\n";
        throw new Error(errMsg);
    }
    auto R = gm.R(Q);
    lua_pushnumber(L, R);
    return 1;
}

extern(C) int gamma(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 2, Q);
    if ( Q.T <= 0.0 || isNaN(Q.T) ) {
        string errMsg = "ERROR: when calling 'gamma'\n";
        errMsg ~= "        The supplied temperature value is negative, 0 or has not been set.\n";
        errMsg ~= "        Check that the 'T' field is set with a valid value.\n";
        errMsg ~= "        The complete gas state is:\n";
        errMsg ~= Q.toString();
        errMsg ~= "\nBailing out\n";
        throw new Error(errMsg);
    }
    auto gamma = gm.gamma(Q);
    lua_pushnumber(L, gamma);
    return 1;
}

extern(C) int molMass(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 2, Q);
    auto molMass = gm.molecular_mass(Q);
    lua_pushnumber(L, molMass);
    return 1;
}

extern(C) int massf2molef(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 2, Q);
    number[] molef; molef.length = gm.n_species;
    gm.massf2molef(Q, molef);
    // Place molef in an array and leave at at
    // top-of-stack as a return to the caller.
    lua_newtable(L);
    foreach (i, mf; molef) {
        lua_pushnumber(L, mf);
        lua_setfield(L, -2, toStringz(gm.species_name(i)));
    }
    return 1;
}

extern(C) int molef2massf(lua_State* L)
{
    int narg = lua_gettop(L);
    auto gm = checkGasModel(L, 1);
    number[] molef;
    molef.length = gm.n_species;
    if ( lua_istable(L, 2) ) {
        getSpeciesValsFromTable(L, gm, 2, molef, "molef");
    } else {
        string errMsg = "Error in call to molef2massf():\n";
        errMsg ~= "The value for 'molef' is not a table of key-val pairs.\n";
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    auto Q = GasState(gm);
    if ( narg == 3 ) {
        getGasStateFromTable(L, gm, 3, Q);
    }
    gm.molef2massf(molef, Q);
    // Update table with new mass fractions
    if ( narg == 3 ) {
        setGasStateInTable(L, gm, 3, Q);
    }
    // and return a table to the caller.
    lua_newtable(L);
    foreach (i, mf; Q.massf) {
        lua_pushnumber(L, mf);
        lua_setfield(L, -2, toStringz(gm.species_name(i)));
    }
    return 1;
}

/**
 * Two forms are supported for this method.
 * 1. A GasState table as input:
 *
 *    gmodel:massf2conc(Q)
 *
 * 2. Density and a mass fractions table as input.
 *
 *    gmodel:massf2conc(rho, massf)
 */
extern(C) int massf2conc(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = GasState(gm);
    number[] conc; conc.length = gm.n_species;
    int nargs = lua_gettop(L);
    if ( nargs == 2 ) { // Support for option 1.
        getGasStateFromTable(L, gm, 2, Q);
    } else {
        number rho = luaL_checknumber(L, 2);
        number[] massf; massf.length = gm.n_species;
        getSpeciesValsFromTable(L, gm, 3, massf, "massf");
        Q.rho = rho;
        Q.massf[] = massf[];
    }
    gm.massf2conc(Q, conc);
    // Place conc in an array and leave at at
    // top-of-stack as a return to the caller.
    lua_newtable(L);
    foreach (i, c; conc) {
        lua_pushnumber(L, c);
        lua_setfield(L, -2, toStringz(gm.species_name(i)));
    }
    return 1;
}

extern(C) int conc2massf(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    number[] conc;
    conc.length = gm.n_species();
    if ( lua_istable(L, 2) ) {
        getSpeciesValsFromTable(L, gm, 2, conc, "conc");
    } else {
        string errMsg = "Error in call to conc2massf():\n";
        errMsg ~= "The value for 'conc' is not a table of key-val pairs.\n";
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 3, Q);
    gm.conc2massf(conc, Q);
    // Update table with new mass fractions
    setGasStateInTable(L, gm, 3, Q);
    // and return a table to the caller.
    lua_newtable(L);
    foreach (i, mf; Q.massf) {
        lua_pushnumber(L, mf);
        lua_setfield(L, -2, toStringz(gm.species_name(i)));
    }
    return 1;
}

extern(C) int binary_diffusion_coefficients(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    size_t nsp = gm.n_species;
    number[][] D;
    D.length = nsp;
    foreach(ref Di; D) Di.length=nsp;

    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 2, Q);
    gm.binary_diffusion_coefficients(Q, D);

    // Return the Diffusion_coefficients to the caller
    size_t idx = 1;
    lua_newtable(L);
    foreach (i, Di; D) {
        lua_newtable(L);
        foreach (j, Dij; Di){
            lua_pushnumber(L, Dij);
            lua_rawseti(L, -2, to!int(j+1));
        }
        lua_rawseti(L, -2, to!int(i+1));
    }
    return 1;
}

void getSpeciesValsFromTable(lua_State* L, GasModel gm, int idx,
                             ref double[] vals, string tabName)
{
    int old_top = lua_gettop(L);
    assert(vals.length == gm.n_species(), "Array for species fractions wrong length.");
    // 1. Check all keys are valid species names.
    lua_pushnil(L);
    while (lua_next(L, idx) != 0) {
        string key = to!string(lua_tostring(L, -2));
        auto isp = gm.species_index(key);
        if (isp == -1) {
            string errMsg = format("Species name used in %s table does not exist: %s\n",
                                   tabName, key);
            lua_pop(L, 1);
            throw new LuaInputException(errMsg);
        }
        lua_pop(L, 1);
    }
    // 2. Now set all values to 0.0
    //    (then we'll correct that in the next step)
    //    [Or to 1.0 in the case of a n_species = 1]
    if ( gm.n_species() == 1 ) {
        vals[0] = 1.0;
    } else {
        vals[] = 0.0;
    }
    // 3. Now find those values that we have explicitly set
    foreach (isp; 0 .. gm.n_species()) {
        lua_getfield(L, -1, toStringz(gm.species_name(isp)));
        if ( lua_isnumber(L, -1) ) {
            vals[isp] = lua_tonumber(L, -1);
        } else if ( lua_isnil(L, -1) ) {
            vals[isp] = 0.0;
        } else {
            string errMsg = format("The value for species '%s' in the %s table is not a number.\n",
                                   gm.species_name(isp), tabName);
            lua_pop(L, 1);
            throw new LuaInputException(errMsg);
        }
        lua_pop(L, 1);
    }
    warnOnStackChange(L, old_top);
} // end getSpeciesValsFromTable()

void getSpeciesValsFromTable(lua_State* L, GasModel gm, int idx,
                             ref Complex!double[] vals, string tabName)
{
    assert(vals.length == gm.n_species(), "Array for species fractions wrong length.");
    double[] dvals; dvals.length = gm.n_species();
    getSpeciesValsFromTable(L, gm, idx, dvals, tabName);
    foreach (i, v; dvals) { vals[i] = Complex!double(v, 0.0); }
} // end getSpeciesValsFromTable()

void checkAndScaleMassFractions(number[] massf, double tol)
{
    auto massfSum = sum(massf);
    if ( fabs(massfSum - 1.0) > tol ) {
        string errMsg = "The given mass fraction values do not sum to 1.0.\n";
        errMsg ~= format("The sum value is: %e\n", massfSum);
        errMsg ~= format("The error in the sum is larger than tolerance= %e \n", tol);
        throw new Error(errMsg);
    }
    // otherwise, do scaling based on normalising by sum
    massf[] /= massfSum;
}

void pushNewGasTable(lua_State* L, ref const(GasState) Q, GasModel gm)
{
    lua_newtable(L);
    int idx = lua_gettop(L);
    setGasStateInTable(L, gm, idx, Q);
    // It's convenient in the Lua code to put a reference
    // to the gas model in the GasState table.
    pushObj!(GasModel, GasModelMT)(L, gm);
    lua_setfield(L, idx, "gasmodel");
}

GasState makeNewGasState(GasModel gm)
{
    GasState Q = GasState(gm);
    // For the special case of n_species = 1, it's quite likely
    // that the user might never consider the massf array.
    // In which case, we'll just set that value to 1.0.
    // For all other cases, we take the default action of
    // setting all mass fractions to 0.0.
    // For a multi-component gas, we expect the user to
    // take care with setting mass fraction values.
    if ( gm.n_species == 1 ) {
        Q.massf[0] = 1.0;
    } else {
        foreach (ref mf; Q.massf) { mf = 0.0; }
    }
    return Q;
}

extern(C) int createTableForGasState(lua_State* L)
{
    GasModel gm = checkGasModel(L, 1);
    GasState Q = makeNewGasState(gm);
    pushNewGasTable(L, Q, gm);
    return 1;
}

extern(C) int newTableForGasState(lua_State* L)
{
    lua_remove(L, 1); // Remove first argument this
    lua_rawgeti(L, 1, 1);
    GasModel gm = checkGasModel(L, -1);
    lua_pop(L, 1);
    GasState Q = makeNewGasState(gm);
    pushNewGasTable(L, Q, gm);
    return 1;
}

void getGasStateFromTable(lua_State* L, GasModel gm, int idx, ref GasState Q)
{
    int old_top = lua_gettop(L);
    lua_getfield(L, idx, "rho");
    if ( lua_isnumber(L, -1) ) {
        Q.rho = lua_tonumber(L, -1);
    } else if ( lua_isnil(L, -1) ) {
        // leave untouched
    } else {
        string errMsg = "The value for 'rho' is not a number.\n";
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "p");
    if ( lua_isnumber(L, -1) ) {
        Q.p = lua_tonumber(L, -1);
    } else if ( lua_isnil(L, -1) ) {
        // leave untouched
    } else {
        string errMsg = "The value for 'p' is not a number.\n";
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "T");
    if ( lua_isnumber(L, -1) ) {
        Q.T = lua_tonumber(L, -1);
    } else if ( lua_isnil(L, -1) ) {
        // leave untouched
    } else {
        string errMsg = "The value for 'T' is not a number.\n";
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "p_e");
    if ( lua_isnumber(L, -1) ) {
        Q.p_e = lua_tonumber(L, -1);
    } else if ( lua_isnil(L, -1) ) {
        // leave untouched
    } else {
        string errMsg = "The value for 'p_e' is not a number.\n";
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "a");
    if ( lua_isnumber(L, -1) ) {
        Q.a = lua_tonumber(L, -1);
    } else if ( lua_isnil(L, -1) ) {
        // leave untouched
    } else {
        string errMsg = "The value for 'a' is not a number.\n";
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "u");
    if ( lua_isnumber(L, -1) ) {
        Q.u = lua_tonumber(L, -1);
    } else if ( lua_isnil(L, -1) ) {
        // leave untouched
    } else {
        string errMsg = "The value for 'u' is not a number.\n";
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "u_modes");
    if ( lua_istable(L, -1) ) {
        auto n = to!int(lua_objlen(L, -1));
        if ( n != Q.u_modes.length ) {
            string errMsg = format("Wrong number of internal energy values " ~
                                   "in GasState table.\n Expected: %d, Given: %d\n",
                                   Q.u_modes.length, n);
            lua_pop(L, 1);
            throw new Error(errMsg);
        }
        foreach ( i; 1..n+1 ) {
            lua_rawgeti(L, -1, i);
            if ( lua_isnumber(L, -1) ) {
                Q.u_modes[i-1] = lua_tonumber(L, -1);
            } else {
                string errMsg = format("The value for 'u_modes[%d]' is not a number.\n", i);
                lua_pop(L, 1);
                throw new Error(errMsg);
            }
            lua_pop(L, 1);
        }
    } else if ( lua_isnil(L, -1) ) {
        // leave untouched
    } else {
        string errMsg = "The value for 'u_modes' is not an array of numbers.\n";
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "T_modes");
    if ( lua_istable(L, -1) ) {
        auto n = to!int(lua_objlen(L, -1));
        if ( n != Q.T_modes.length ) {
            string errMsg = format("Wrong number of internal temperature values " ~
                                   " in GasState table.\n Expected: %d, Given: %d\n",
                                   Q.T_modes.length, n);
            lua_pop(L, 1);
            throw new Error(errMsg);
        }
        foreach ( i; 1..n+1 ) {
            lua_rawgeti(L, -1, i);
            if ( lua_isnumber(L, -1) ) {
                Q.T_modes[i-1] = lua_tonumber(L, -1);
            } else {
                string errMsg = format("The value for 'T_modes[%d]' is not a number.\n", i);
                lua_pop(L, 1);
                throw new Error(errMsg);
            }
            lua_pop(L, 1);
        }
    } else if ( lua_isnil(L, -1) ) {
        // leave untouched
    } else {
        string errMsg = "The value for 'T_modes' is not an array of numbers.\n";
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "mu");
    if ( lua_isnumber(L, -1) ) {
        Q.mu = lua_tonumber(L, -1);
    } else if ( lua_isnil(L, -1) ) {
        // leave untouched
    } else {
        string errMsg = "The value for 'mu' is not a number.\n";
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "k");
    if ( lua_isnumber(L, -1) ) {
        Q.k = lua_tonumber(L, -1);
    } else if ( lua_isnil(L, -1) ) {
        // leave untouched
    } else {
        string errMsg = "The value for 'k' is not a number.\n";
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "k_modes");
    if ( lua_istable(L, -1) ) {
        auto n = to!int(lua_objlen(L, -1));
        if ( n != Q.k_modes.length ) {
            string errMsg = format("Wrong number of internal thermal conductivity values " ~
                                   "in GasState table.\n Expected: %d, Given: %d\n",
                                   Q.k_modes.length, n);
            lua_pop(L, 1);
            throw new Error(errMsg);
        }
        foreach ( i; 1..n+1 ) {
            lua_rawgeti(L, -1, i);
            if ( lua_isnumber(L, -1) ) {
                Q.k_modes[i-1] = lua_tonumber(L, -1);
            } else {
                string errMsg = format("The value for 'k_modes[%d]' is not a number.\n", i);
                lua_pop(L, 1);
                throw new Error(errMsg);
            }
            lua_pop(L, 1);
        }
    } else if ( lua_isnil(L, -1) ) {
        // leave untouched
    } else {
        string errMsg = "The value for 'k_modes' is not an array of numbers.\n";
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "sigma");
    if ( lua_isnumber(L, -1) ) {
        Q.sigma = lua_tonumber(L, -1);
    } else if ( lua_isnil(L, -1) ) {
        // leave untouched
    } else {
        string errMsg = "The value for 'sigma' is not a number.\n";
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "massf");
    if ( lua_istable(L, -1) ) {
        int massfIdx = lua_gettop(L);
        getSpeciesValsFromTable(L, gm, massfIdx, Q.massf, "massf");
        checkAndScaleMassFractions(Q.massf, MASSF_ERROR_TOL);
    } else if ( lua_isnil(L, -1) ) {
        // leave untouched
    } else {
        string errMsg = "The value for 'massf' is not a table of key-val pairs.\n";
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "quality");
    if ( lua_isnumber(L, -1) ) {
        Q.quality = lua_tonumber(L, -1);
    } else if ( lua_isnil(L, -1) ) {
        // leave untouched
    } else {
        string errMsg = "The value for 'quality' is not a number.\n";
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "ceaSavedData");
    if ( lua_istable(L, -1) && Q.ceaSavedData ) {
        int ceaSavedDataIdx = lua_gettop(L);
        lua_getfield(L, ceaSavedDataIdx, "p");
        if ( lua_isnumber(L, -1) ) { Q.ceaSavedData.p = lua_tonumber(L, -1); }
        lua_pop(L, 1);
        lua_getfield(L, ceaSavedDataIdx, "rho");
        if ( lua_isnumber(L, -1) ) { Q.ceaSavedData.rho = lua_tonumber(L, -1); }
        lua_pop(L, 1);
        lua_getfield(L, ceaSavedDataIdx, "u");
        if ( lua_isnumber(L, -1) ) { Q.ceaSavedData.u = lua_tonumber(L, -1); }
        lua_pop(L, 1);
        lua_getfield(L, ceaSavedDataIdx, "h");
        if ( lua_isnumber(L, -1) ) { Q.ceaSavedData.h = lua_tonumber(L, -1); }
        lua_pop(L, 1);
        lua_getfield(L, ceaSavedDataIdx, "T");
        if ( lua_isnumber(L, -1) ) { Q.ceaSavedData.T = lua_tonumber(L, -1); }
        lua_pop(L, 1);
        lua_getfield(L, ceaSavedDataIdx, "a");
        if ( lua_isnumber(L, -1) ) { Q.ceaSavedData.a = lua_tonumber(L, -1); }
        lua_pop(L, 1);
        lua_getfield(L, ceaSavedDataIdx, "Mmass");
        if ( lua_isnumber(L, -1) ) { Q.ceaSavedData.Mmass = lua_tonumber(L, -1); }
        lua_pop(L, 1);
        lua_getfield(L, ceaSavedDataIdx, "Rgas");
        if ( lua_isnumber(L, -1) ) { Q.ceaSavedData.Rgas = lua_tonumber(L, -1); }
        lua_pop(L, 1);
        lua_getfield(L, ceaSavedDataIdx, "gamma");
        if ( lua_isnumber(L, -1) ) { Q.ceaSavedData.gamma = lua_tonumber(L, -1); }
        lua_pop(L, 1);
        lua_getfield(L, ceaSavedDataIdx, "Cp");
        if ( lua_isnumber(L, -1) ) { Q.ceaSavedData.Cp = lua_tonumber(L, -1); }
        lua_pop(L, 1);
        lua_getfield(L, ceaSavedDataIdx, "s");
        if ( lua_isnumber(L, -1) ) { Q.ceaSavedData.s = lua_tonumber(L, -1); }
        lua_pop(L, 1);
        lua_getfield(L, ceaSavedDataIdx, "k");
        if ( lua_isnumber(L, -1) ) { Q.ceaSavedData.k = lua_tonumber(L, -1); }
        lua_pop(L, 1);
        lua_getfield(L, ceaSavedDataIdx, "mu");
        if ( lua_isnumber(L, -1) ) { Q.ceaSavedData.mu = lua_tonumber(L, -1); }
        lua_pop(L, 1);
        // We skip the table of mass fractions within the ceaSavedData.
    }
    lua_pop(L, 1);
    warnOnStackChange(L, old_top);
} // end getGasStateFromTable()

void setGasStateInTable(lua_State* L, GasModel gm, int idx, const(GasState) Q)
{
    int old_top = lua_gettop(L);
    lua_pushnumber(L, Q.rho);
    lua_setfield(L, idx, "rho");

    lua_pushnumber(L, Q.p);
    lua_setfield(L, idx, "p");

    lua_pushnumber(L, Q.T);
    lua_setfield(L, idx, "T");

    lua_pushnumber(L, Q.p_e);
    lua_setfield(L, idx, "p_e");

    lua_pushnumber(L, Q.a);
    lua_setfield(L, idx, "a");

    lua_pushnumber(L, Q.u);
    lua_setfield(L, idx, "u");

    lua_newtable(L);
    foreach (i, e; Q.u_modes) {
        lua_pushnumber(L, e); lua_rawseti(L, -2, to!int(i)+1);
    }
    lua_setfield(L, idx, "u_modes");

    lua_newtable(L);
    foreach (i, T; Q.T_modes) {
        lua_pushnumber(L, T); lua_rawseti(L, -2, to!int(i)+1);
    }
    lua_setfield(L, idx, "T_modes");

    lua_pushnumber(L, Q.mu);
    lua_setfield(L, idx, "mu");

    lua_pushnumber(L, Q.k);
    lua_setfield(L, idx, "k");

    lua_newtable(L);
    foreach (i, k; Q.k_modes) {
        lua_pushnumber(L, k); lua_rawseti(L, -2, to!int(i)+1);
    }
    lua_setfield(L, idx, "k_modes");

    lua_pushnumber(L, Q.sigma);
    lua_setfield(L, idx, "sigma");

    lua_newtable(L);
    foreach (i, mf; Q.massf) {
        lua_pushnumber(L, mf);
        string spName = gm.species_name(i);
        lua_setfield(L, -2, toStringz(spName));
    }
    lua_setfield(L, idx, "massf");

    lua_pushnumber(L, Q.quality);
    lua_setfield(L, idx, "quality");

    // We will also add a special field for the ceaSavedData
    // if we find at non-null pointer to the saved data struct.
    // Yes, it will appear that we are duplicating some of
    // the information stored directly in the GasState table
    // but this will give us direct access to the data that
    // we've scanned from the CEA output.
    if ( Q.ceaSavedData ) {
        lua_newtable(L);

        lua_pushnumber(L, Q.ceaSavedData.p);
        lua_setfield(L, -2, "p");

        lua_pushnumber(L, Q.ceaSavedData.rho);
        lua_setfield(L, -2, "rho");

        lua_pushnumber(L, Q.ceaSavedData.u);
        lua_setfield(L, -2, "u");

        lua_pushnumber(L, Q.ceaSavedData.h);
        lua_setfield(L, -2, "h");

        lua_pushnumber(L, Q.ceaSavedData.T);
        lua_setfield(L, -2, "T");

        lua_pushnumber(L, Q.ceaSavedData.a);
        lua_setfield(L, -2, "a");

        lua_pushnumber(L, Q.ceaSavedData.Mmass);
        lua_setfield(L, -2, "Mmass");

        lua_pushnumber(L, Q.ceaSavedData.Rgas);
        lua_setfield(L, -2, "Rgas");

        lua_pushnumber(L, Q.ceaSavedData.gamma);
        lua_setfield(L, -2, "gamma");

        lua_pushnumber(L, Q.ceaSavedData.Cp);
        lua_setfield(L, -2, "Cp");

        lua_pushnumber(L, Q.ceaSavedData.s);
        lua_setfield(L, -2, "s");

        lua_pushnumber(L, Q.ceaSavedData.k);
        lua_setfield(L, -2, "k");

        lua_pushnumber(L, Q.ceaSavedData.mu);
        lua_setfield(L, -2, "mu");

        lua_newtable(L);
        foreach (sp, mf; Q.ceaSavedData.massf) {
            lua_pushnumber(L, mf);
            lua_setfield(L, -2, toStringz(sp));
        }
        lua_setfield(L, -2, "massf");

        lua_setfield(L, idx, "ceaSavedData");
    }
    warnOnStackChange(L, old_top);
} // end setGasStateInTable()

extern(C) int printValues(lua_State* L)
{
    lua_getfield(L, 1, "gasmodel");
    auto gm = checkGasModel(L, -1);
    lua_pop(L, 1);
    auto Q = GasState(gm);
    getGasStateFromTable(L, gm, 1, Q);
    writeln(Q.toString);
    return 0;
}

extern(C) int copyValues(lua_State* L)
{
    lua_getfield(L, 1, "gasmodel");
    auto gm = checkGasModel(L, -1);
    lua_pop(L, 1);

    // Source GasState table
    auto Q0 = GasState(gm);
    getGasStateFromTable(L, gm, 1, Q0);

    // Target GasState table
    // Assume table is available at index 2
    setGasStateInTable(L, gm, 2, Q0);

    return 0;
}

void registerGasModel(lua_State* L)
{
    luaL_newmetatable(L, GasModelMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates current metatable
    lua_setfield(L, -2, "__index");
    /* Register methods for use. */
    lua_pushcfunction(L, &newGasModel);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &toStringObj!(GasModel, GasModelMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &nSpecies);
    lua_setfield(L, -2, "nSpecies");
    lua_pushcfunction(L, &nModes);
    lua_setfield(L, -2, "nModes");
    lua_pushcfunction(L, &molMasses);
    lua_setfield(L, -2, "molMasses");
    lua_pushcfunction(L, &charge);
    lua_setfield(L, -2, "charge");
    lua_pushcfunction(L, &speciesName);
    lua_setfield(L, -2, "speciesName");
    lua_pushcfunction(L, &speciesIndex);
    lua_setfield(L, -2, "speciesIndex");
    lua_pushcfunction(L, &createTableForGasState);
    lua_setfield(L, -2, "createGasState");
    lua_pushcfunction(L, &thermoPU);
    lua_setfield(L, -2, "updateThermoFromPU");
    lua_pushcfunction(L, &thermoPT);
    lua_setfield(L, -2, "updateThermoFromPT");
    lua_pushcfunction(L, &thermoRHOU);
    lua_setfield(L, -2, "updateThermoFromRHOU");
    lua_pushcfunction(L, &thermoRHOU);
    lua_setfield(L, -2, "updateThermoFromRHOE"); // keep the old name, as well
    lua_pushcfunction(L, &thermoRHOT);
    lua_setfield(L, -2, "updateThermoFromRHOT");
    lua_pushcfunction(L, &thermoRHOP);
    lua_setfield(L, -2, "updateThermoFromRHOP");
    lua_pushcfunction(L, &thermoPS);
    lua_setfield(L, -2, "updateThermoFromPS");
    lua_pushcfunction(L, &thermoHS);
    lua_setfield(L, -2, "updateThermoFromHS");
    lua_pushcfunction(L, &soundSpeed);
    lua_setfield(L, -2, "updateSoundSpeed");
    lua_pushcfunction(L, &transCoeffs);
    lua_setfield(L, -2, "updateTransCoeffs");
    lua_pushcfunction(L, &dpdrhoConstT);
    lua_setfield(L, -2, "dpdrhoConstT");
    lua_pushcfunction(L, &intEnergy);
    lua_setfield(L, -2, "intEnergy");
    lua_pushcfunction(L, &enthalpy);
    lua_setfield(L, -2, "enthalpy");
    lua_pushcfunction(L, &entropy);
    lua_setfield(L, -2, "entropy");
    lua_pushcfunction(L, &gibbsFreeEnergy);
    lua_setfield(L, -2, "gibbsFreeEnergy");
    lua_pushcfunction(L, &Cv);
    lua_setfield(L, -2, "Cv");
    lua_pushcfunction(L, &Cv);
    lua_setfield(L, -2, "dudTConstV");
    lua_pushcfunction(L, &Cp);
    lua_setfield(L, -2, "Cp");
    lua_pushcfunction(L, &Cp);
    lua_setfield(L, -2, "dhdTConstP");
    lua_pushcfunction(L, &R);
    lua_setfield(L, -2, "R");
    lua_pushcfunction(L, &R);
    lua_setfield(L, -2, "gasConstant");
    lua_pushcfunction(L, &gamma);
    lua_setfield(L, -2, "gamma");
    lua_pushcfunction(L, &molMass);
    lua_setfield(L, -2, "molMass");
    lua_pushcfunction(L, &massf2molef);
    lua_setfield(L, -2, "massf2molef");
    lua_pushcfunction(L, &molef2massf);
    lua_setfield(L, -2, "molef2massf");
    lua_pushcfunction(L, &massf2conc);
    lua_setfield(L, -2, "massf2conc");
    lua_pushcfunction(L, &conc2massf);
    lua_setfield(L, -2, "conc2massf");
    lua_pushcfunction(L, &binary_diffusion_coefficients);
    lua_setfield(L, -2, "binary_diffusion_coefficients");

    // Make class visible
    lua_setglobal(L, GasModelMT.toStringz);

    // Make initialisation of GasState table look like a class constructor
    luaL_newmetatable(L, "GasState");

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates current metatable
    lua_setfield(L, -2, "__index");
    /* Register methods for use. */
    lua_pushcfunction(L, &newTableForGasState);
    lua_setfield(L, -2, "new");

    // Make GasState constructor visible
    lua_setglobal(L, "GasState");

    // Set a global function to print values in GasState
    lua_pushcfunction(L, &printValues);
    lua_setglobal(L, "printValues");
    // Set a global function to copy GasState values from one table to another.
    lua_pushcfunction(L, &copyValues);
    lua_setglobal(L, "copyValues");

    // Set some of the physical constants as global
    lua_pushnumber(L, P_atm);
    lua_setglobal(L, "P_atm");
    lua_pushnumber(L, R_universal);
    lua_setglobal(L, "R_universal");
}
