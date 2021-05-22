// luaglobalconfig.d
// Lua access to the GlobalConfig class data, for use in the preparation script.
//
// Peter J. and Rowan G.
// 2015-03-02: First code adapted from the other lua wrapper modules.

module luaglobalconfig;

import core.stdc.stdlib : exit;

import std.stdio;
import std.string;
import std.conv;
import util.lua;
import util.lua_service;

import gas;
import gas.luagas_model;
import globalconfig;


//------------------------------------------------------------------------
// Functions related to the managed gas model.

extern(C) int setGasModel(lua_State* L)
{
    if (lua_isstring(L, 1)) {
        string fname = to!string(luaL_checkstring(L, 1));
        GlobalConfig.gas_model_file = fname;
        // 2021-05-22 Also, to set config.gas_model_name in Lua world.
        lua_getglobal(L, "config".toStringz);
        if (lua_istable(L, -1)) {
            lua_pushstring(L, fname.toStringz);
            lua_setfield(L, -2, "gas_model_file".toStringz);
        }
        lua_pop(L, 1); // dispose of config
        try {
            GlobalConfig.gmodel_master = init_gas_model(fname);
        } catch (GasModelException e) {
            string msg = "\nThere is a problem in call to setGasModel. Reported errors are:\n";
            msg ~= e.msg;
            msg ~= "\n---------------------------\n";
            msg ~= "The preparation stage cannot proceed. Exiting without completing.\n";
            writeln(msg);
            exit(1);
        }
        lua_settop(L, 0); // clear the stack
        lua_pushinteger(L, GlobalConfig.gmodel_master.n_species);
        lua_pushinteger(L, GlobalConfig.gmodel_master.n_modes);
        GasModelStore ~= pushObj!(GasModel, GasModelMT)(L, GlobalConfig.gmodel_master);
        return 3;
    } else {
        string msg = "setGasModel expects a string as the name of the gas model file.";
        luaL_error(L, msg.toStringz);
        return 0;
    }
}

extern(C) int getGasModel(lua_State* L)
{
    if (GlobalConfig.gmodel_master is null) {
        string msg = "The master gas model appears to be uninitialized. "~
            "You should initialize it with setGasModel().";
        luaL_error(L, msg.toStringz);
    }
    lua_settop(L, 0); // clear the stack
    pushObj!(GasModel, GasModelMT)(L, GlobalConfig.gmodel_master);
    return 1;
}

//-----------------------------------------------------------------------
// Call the following function from the main program to get the
// functions appearing in the Lua interpreter.

void registerGlobalConfig(lua_State* L)
{
    // Register global functions related to the managed gas model.
    lua_pushcfunction(L, &setGasModel);
    lua_setglobal(L, "setGasModel");
    lua_pushcfunction(L, &getGasModel);
    lua_setglobal(L, "getGasModel");
} // end registerGlobalConfig()
