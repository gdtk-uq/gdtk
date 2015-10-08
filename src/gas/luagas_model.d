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
    // Make class visible
    lua_setglobal(L, GasModelMT.toStringz);
}

