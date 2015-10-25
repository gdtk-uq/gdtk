/**
 * luagas_model.d
 * Lua interface to the gas model.
 *
 * Author: Rowan G.
 * Version: Initial cut.
 */

module gas.luagas_model;

import std.stdio;
import std.conv;
import std.string;
import util.lua;
import util.lua_service;

import gas.gas_model;
import gas.gas_model_util;


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
 * Construction of GasModel object in Lua will accept
 * a filename.
 * -------------
 * gmodel = GasModel:new{fname}
 *--------------
 */
extern(C) int newGasModel(lua_State* L)
{
    lua_remove(L, 1); // Remove first argument this
    
    string fname = to!string(luaL_checkstring(L, 1));
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
    foreach ( int isp, mMass; gm.mol_masses ) {
	lua_pushnumber(L, mMass);
	lua_rawseti(L, -2, isp+1);
    }
    return 1;
}

extern(C) int speciesName(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    int i = to!int(lua_tointeger(L, 2));
    // Use "i-1" since D indexes from 0
    lua_pushstring(L, toStringz(gm.species_name(i-1)));
    return 1;
}

extern(C) int thermoPT(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, 2, Q);
    gm.update_thermo_from_pT(Q);
    setGasStateInTable(L, 2, Q);
    return 0;
}

extern(C) int thermoRHOE(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, 2, Q);
    gm.update_thermo_from_rhoe(Q);
    setGasStateInTable(L, 2, Q);
    return 0;
}

extern(C) int thermoRHOT(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, 2, Q);
    gm.update_thermo_from_rhoT(Q);
    setGasStateInTable(L, 2, Q);
    return 0;
}

extern(C) int thermoRHOP(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, 2, Q);
    gm.update_thermo_from_rhop(Q);
    setGasStateInTable(L, 2, Q);
    return 0;
}

extern(C) int thermoPS(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, 2, Q);
    auto s = lua_tonumber(L, 3);
    gm.update_thermo_from_ps(Q, s);
    setGasStateInTable(L, 2, Q);
    return 0;
}

extern(C) int thermoHS(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, 2, Q);
    auto h = lua_tonumber(L, 3);
    auto s = lua_tonumber(L, 4);
    gm.update_thermo_from_hs(Q, h, s);
    setGasStateInTable(L, 2, Q);
    return 0;
}

extern(C) int soundSpeed(lua_State* L)
{
   auto gm = checkGasModel(L, 1);
   auto Q = new GasState(gm.n_species, gm.n_modes);
   getGasStateFromTable(L, 2, Q); 
   gm.update_sound_speed(Q);
   setGasStateInTable(L, 2, Q);
   return 0;
}

// We don't wrap dedT_const_v as we do Cv in its place.
// We don't wrap dhdT_const_p as we do Cp in its place.

extern(C) int dpdrhoConstT(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, 2, Q);
    auto dpdrho = gm.dpdrho_const_T(Q);
    lua_pushnumber(L, dpdrho);
    return 1;
}

// We don't wrap gas_constant as we do R in its place.

extern(C) int intEnergy(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, 2, Q);
    auto e = gm.internal_energy(Q);
    lua_pushnumber(L, e);
    return 1;
}

extern(C) int enthalpy(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, 2, Q);
    auto h = gm.enthalpy(Q);
    lua_pushnumber(L, h);
    return 1;
}

extern(C) int entropy(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, 2, Q);
    auto s = gm.entropy(Q);
    lua_pushnumber(L, s);
    return 1;
}

extern(C) int Cv(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, 2, Q);
    auto Cv = gm.Cv(Q);
    lua_pushnumber(L, Cv);
    return 1;
}

extern(C) int Cp(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, 2, Q);
    auto Cp = gm.Cp(Q);
    lua_pushnumber(L, Cp);
    return 1;
}

extern(C) int R(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, 2, Q);
    auto R = gm.R(Q);
    lua_pushnumber(L, R);
    return 1;
}

extern(C) int gamma(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, 2, Q);
    auto gamma = gm.gamma(Q);
    lua_pushnumber(L, gamma);
    return 1;
}

extern(C) int molMass(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, 2, Q);
    auto molMass = gm.molecular_mass(Q);
    lua_pushnumber(L, molMass);
    return 1;
}

extern(C) int massf2molef(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    double[] molef; molef.length = gm.n_species;
    gm.massf2molef(Q, molef);
    // Place molef in an array and leave at at
    // top-of-stack as a return to the caller.
    lua_newtable(L);
    foreach ( int i, mf; molef ) {
	lua_pushnumber(L, mf); lua_rawseti(L, -2, i+1);
    }
    return 1;
}

extern(C) int molef2massf(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    double[] molef;
    auto n = to!int(lua_objlen(L, 2));
    foreach ( isp; 1 .. n+1 ) {
	lua_rawgeti(L, 2, isp);
	molef ~= lua_tonumber(L, -1);
	lua_pop(L, 1);
    }
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, 3, Q);
    gm.molef2massf(molef, Q);
    // Update table with new mass fractions
    setGasStateInTable(L, 3, Q);
    // and return a table to the caller.
    lua_newtable(L);
    foreach ( int i, mf; Q.massf ) {
	lua_pushnumber(L, mf); lua_rawseti(L, -2, i+1);
    }
    return 1;
}

extern(C) int massf2conc(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    double[] conc; conc.length = gm.n_species;
    gm.massf2conc(Q, conc);
    // Place conc in an array and leave at at
    // top-of-stack as a return to the caller.
    lua_newtable(L);
    foreach ( int i, c; conc ) {
	lua_pushnumber(L, c); lua_rawseti(L, -2, i+1);
    }
    return 1;
}

extern(C) int conc2massf(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    double[] conc;
    auto n = to!int(lua_objlen(L, 2));
    foreach ( isp; 1 .. n+1 ) {
	lua_rawgeti(L, 2, isp);
	conc ~= lua_tonumber(L, -1);
	lua_pop(L, 1);
    }
    auto Q = new GasState(gm.n_species, gm.n_modes);
    getGasStateFromTable(L, 3, Q);
    gm.conc2massf(conc, Q);
    // Update table with new mass fractions
    setGasStateInTable(L, 3, Q);
    // and return a table to the caller.
    lua_newtable(L);
    foreach ( int i, mf; Q.massf ) {
	lua_pushnumber(L, mf); lua_rawseti(L, -2, i+1);
    }
    return 1;
}

extern(C) int createTableForGasState(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    lua_newtable(L);
    int idx = lua_gettop(L);
    setGasStateInTable(L, idx, Q);
    return 1;
}

extern(C) int newTableForGasState(lua_State* L)
{
    lua_remove(L, 1); // Remove first argument this
    lua_rawgeti(L, 1, 1);
    auto gm = checkGasModel(L, -1);
    lua_pop(L, 1);
    auto Q = new GasState(gm.n_species, gm.n_modes);
    lua_newtable(L);
    int idx = lua_gettop(L);
    setGasStateInTable(L, idx, Q);
    return 1;
}

void getGasStateFromTable(lua_State* L, int idx, GasState Q)
{
    lua_getfield(L, idx, "rho");
    if ( lua_isnumber(L, -1) ) {
	Q.rho = lua_tonumber(L, -1);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "p");
    if ( lua_isnumber(L, -1) ) {
	Q.p = lua_tonumber(L, -1);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "p_e");
    if ( lua_isnumber(L, -1) ) {
	Q.p_e = lua_tonumber(L, -1);
    }
    lua_pop(L, 1);
    
    lua_getfield(L, idx, "e");
    if ( lua_istable(L, -1) ) {
	auto n = to!int(lua_objlen(L, -1));
	if ( n != Q.e.length ) {
	    string errMsg = format("Wrong number of energy values in GasState table.\n Expected: %d, Given: %d\n", Q.e.length, n);
	    lua_pop(L, 1);
	    throw new Error(errMsg);
	}
	foreach ( i; 1..n+1 ) {
	    lua_rawgeti(L, -1, i); Q.e[i-1] = lua_tonumber(L, -1);
	    lua_pop(L, 1);
	}
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "T");
    if ( lua_istable(L, -1) ) {
	auto n = to!int(lua_objlen(L, -1));
	if ( n != Q.T.length ) {
	    string errMsg = format("Wrong number of temperature values in GasState table.\n Expected: %d, Given: %d\n", Q.T.length, n);
	    lua_pop(L, 1);
	    throw new Error(errMsg);
	}
	foreach ( i; 1..n+1 ) {
	    lua_rawgeti(L, -1, i); Q.T[i-1] = lua_tonumber(L, -1);
	    lua_pop(L, 1);
	}
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "mu");
    if ( lua_isnumber(L, -1) ) {
	Q.mu = lua_tonumber(L, -1);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "k");
    if ( lua_istable(L, -1) ) {
	auto n = to!int(lua_objlen(L, -1));
	if ( n != Q.T.length ) {
	    string errMsg = format("Wrong number of thermal conductivity values in GasState table.\n Expected: %d, Given: %d\n", Q.k.length, n);
	    lua_pop(L, 1);
	    throw new Error(errMsg);
	}
	foreach ( i; 1..n+1 ) {
	    lua_rawgeti(L, -1, i); Q.k[i-1] = lua_tonumber(L, -1);
	    lua_pop(L, 1);
	}
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "sigma");
    if ( lua_isnumber(L, -1) ) {
	Q.sigma = lua_tonumber(L, -1);
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "massf");
    if ( lua_istable(L, -1) ) {
	auto n = to!int(lua_objlen(L, -1));
	if ( n != Q.massf.length ) {
	    string errMsg = format("Wrong number of mass fraction values in GasState table.\n Expected: %d, Given: %d\n", Q.massf.length, n);
	    lua_pop(L, 1);
	    throw new Error(errMsg);
	}
	foreach ( i; 1..n+1 ) {
	    lua_rawgeti(L, -1, i); Q.massf[i-1] = lua_tonumber(L, -1);
	    lua_pop(L, 1);
	}
    }
    lua_pop(L, 1);

    lua_getfield(L, idx, "quality");
    if ( lua_isnumber(L, -1) ) {
	Q.quality = lua_tonumber(L, -1);
    }
    lua_pop(L, 1);
}

void setGasStateInTable(lua_State* L, int idx, GasState Q)
{
    lua_pushnumber(L, Q.rho);
    lua_setfield(L, idx, "rho");

    lua_pushnumber(L, Q.p);
    lua_setfield(L, idx, "p");

    lua_pushnumber(L, Q.p_e);
    lua_setfield(L, idx, "p_e");

    lua_pushnumber(L, Q.a);
    lua_setfield(L, idx, "a");

    lua_newtable(L);
    foreach ( int i, e; Q.e ) {
	lua_pushnumber(L, e); lua_rawseti(L, -2, i+1);
    }
    lua_setfield(L, idx, "e");

    lua_newtable(L);
    foreach ( int i, T; Q.T ) {
	lua_pushnumber(L, T); lua_rawseti(L, -2, i+1);
    }
    lua_setfield(L, idx, "T");
    
    lua_pushnumber(L, Q.mu);
    lua_setfield(L, idx, "mu");

    lua_newtable(L);
    foreach ( int i, k; Q.k ) {
	lua_pushnumber(L, k); lua_rawseti(L, -2, i+1);
    }
    lua_setfield(L, idx, "k");

    lua_pushnumber(L, Q.sigma);
    lua_setfield(L, idx, "sigma");

    lua_newtable(L);
    foreach ( int i, mf; Q.massf ) {
	lua_pushnumber(L, mf); lua_rawseti(L, -2, i+1);
    }
    lua_setfield(L, idx, "massf");

    lua_pushnumber(L, Q.quality);
    lua_setfield(L, idx, "quality");
}

void registerGasModel(lua_State* L, int tblIdx)
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
    lua_pushcfunction(L, &speciesName);
    lua_setfield(L, -2, "speciesName");
    lua_pushcfunction(L, &createTableForGasState);
    lua_setfield(L, -2, "createGasState");
    lua_pushcfunction(L, &thermoPT);
    lua_setfield(L, -2, "updateThermoFromPT");
    lua_pushcfunction(L, &thermoRHOE);
    lua_setfield(L, -2, "updateThermoFromRHOE");
    lua_pushcfunction(L, &thermoRHOT);
    lua_setfield(L, -2, "updateThermoFromRHOT");
    lua_pushcfunction(L, &thermoPT);
    lua_setfield(L, -2, "updateThermoFromRHOP");
    lua_pushcfunction(L, &thermoPS);
    lua_setfield(L, -2, "updateThermoFromPS");
    lua_pushcfunction(L, &thermoHS);
    lua_setfield(L, -2, "updateThermoFromHS");
    lua_pushcfunction(L, &soundSpeed);
    lua_setfield(L, -2, "updateSoundSpeed");
    lua_pushcfunction(L, &dpdrhoConstT);
    lua_setfield(L, -2, "dpdrhoConstT");
    lua_pushcfunction(L, &intEnergy);
    lua_setfield(L, -2, "intEnergy");
    lua_pushcfunction(L, &enthalpy);
    lua_setfield(L, -2, "enthalpy");
    lua_pushcfunction(L, &entropy);
    lua_setfield(L, -2, "entropy");
    lua_pushcfunction(L, &Cv);
    lua_setfield(L, -2, "Cv");
    lua_pushcfunction(L, &Cv);
    lua_setfield(L, -2, "dedTConstV");
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

    // Make class visible
    lua_setfield(L, tblIdx, GasModelMT.toStringz);

    // Make initialisation of GasState table look like a class constructor
    luaL_newmetatable(L, "GasState");

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates current metatable
    lua_setfield(L, -2, "__index");
    /* Register methods for use. */
    lua_pushcfunction(L, &newTableForGasState);
    lua_setfield(L, -2, "new");

    // Make GasState constructor visible
    lua_setfield(L, tblIdx, "GasState");

}

version(gas_calc) {
    int main(string[] args) {
	if ( args.length != 2 ) {
	    writeln("ERROR: Wrong number of arguments.");
	    writeln("");
	    writeln("Usage: ");
	    writeln("  gas-calc inputFile");
	    return 1;
	}
	auto L = luaL_newstate();
	luaL_openlibs(L);
	registerGasModel(L, LUA_GLOBALSINDEX);

	if ( luaL_dofile(L, toStringz(args[1])) != 0 ) {
	    writeln(to!string(lua_tostring(L, -1)));
	    return 1;
	}
	return 0;
    }
}

