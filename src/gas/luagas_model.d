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

extern(C) int createTableForGasState(lua_State* L)
{
    auto gm = checkGasModel(L, 1);
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
    lua_pushcfunction(L, &speciesName);
    lua_setfield(L, -2, "speciesName");
    lua_pushcfunction(L, &createTableForGasState);
    lua_setfield(L, -2, "createGasState");
    lua_pushcfunction(L, &thermoPT);
    lua_setfield(L, -2, "updateThermoFromPT");
    // Make class visible
    lua_setglobal(L, GasModelMT.toStringz);

    // Global functions
    lua_pushcfunction(L, &createTableForGasState);
    lua_setglobal(L, "GasState");

}

