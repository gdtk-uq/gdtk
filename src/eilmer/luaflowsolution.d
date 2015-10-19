/**
 * luaflowsolution.d
 * Lua interface to the FlowSolution object for use in postprocessing.
 *
 * Authors: Peter Ja dna Rowan G.
 * Version: 2015-10-20 Initial cut.
 */

module luaflowsolution;

import std.array;
import std.format;
import std.stdio;
import std.string;
import std.conv;
import std.traits;
import util.lua;
import util.lua_service;
import globalconfig;
import luaglobalconfig;
import flowsolution;

/// name for FlowSolution object in Lua scripts.
immutable string FlowSolutionMT = "FlowSolution"; 

static const(FlowSolution)[] flowSolutionStore;

FlowSolution checkFlowSolution(lua_State* L, int index)
{
    return checkObj!(FlowSolution, FlowSolutionMT)(L, index);
}

/** 
 * This function implements our constructor for the Lua interface.
 *
 * Construction of a FlowSolution object from in Lua will accept:
 * fs = FlowSolution:new{jobName=<string>, dir=<string>, tindx=<int>, nBlocks=<int>}
 */
extern(C) int newFlowSolution(lua_State* L)
{
    auto managedGasModel = GlobalConfig.gmodel_master;
    if ( managedGasModel is null ) {
	string errMsg = `Error in call to FlowSolution:new.
It appears that you have not yet set the GasModel.
Be sure to call setGasModel(fname) before using a FlowSolution object.`;
	luaL_error(L, errMsg.toStringz);
    }

    lua_remove(L, 1); // Remove first argument "this".
    if ( !lua_istable(L, 1) ) {
	string errMsg = "Error in call to FlowSolution:new. A table is expected as first argument.";
	luaL_error(L, errMsg.toStringz);
    }

    string jobName;
    lua_getfield(L, 1, "jobName");
    if ( lua_isnil(L, -1) ) {
	jobName = "";
	string errMsg = "Error in call to FlowSolution:new.";
	errMsg ~= " A field for jobName was not found.";
	throw new LuaInputException(errMsg);
    } else if ( lua_isstring(L, -1) ) {
	jobName = to!string(luaL_checkstring(L, -1));
    } else {
	string errMsg = "Error in call to FlowSolution:new.";
	errMsg ~= " A field for jobName was found, but the content was not valid.";
	errMsg ~= " The jobName should be given as a string.";
	throw new LuaInputException(errMsg);
    }
    lua_pop(L, 1);

    string dir;
    lua_getfield(L, 1, "dir");
    if ( lua_isnil(L, -1) ) {
	dir = ".";
    } else if ( lua_isstring(L, -1) ) {
	dir = to!string(luaL_checkstring(L, -1));
    } else {
	string errMsg = "Error in call to FlowSolution:new.";
	errMsg ~= " A field for dir was found, but the content was not valid.";
	errMsg ~= " The dir field, if given, should be a string.";
	throw new LuaInputException(errMsg);
    }
    lua_pop(L, 1);

    int tindx;
    lua_getfield(L, 1, "tindx");
    if ( lua_isnil(L, -1) ) {
	tindx = 0;
    } else if ( lua_isnumber(L, -1) ) {
	tindx = to!int(luaL_checknumber(L, -1));
    } else {
	string errMsg = "Error in call to FlowSolution:new.";
	errMsg ~= " A field for tindx was found, but the content was not valid.";
	errMsg ~= " The tindx field, if given, should be an integer.";
	throw new LuaInputException(errMsg);
    }
    lua_pop(L, 1);

    int nBlocks;
    lua_getfield(L, 1, "nBlocks");
    if ( lua_isnil(L, -1) ) {
	nBlocks = 1;
    } else if ( lua_isnumber(L, -1) ) {
	nBlocks = to!int(luaL_checknumber(L, -1));
    } else {
	string errMsg = "Error in call to FlowSolution:new.";
	errMsg ~= " A field for nBlocks was found, but the content was not valid.";
	errMsg ~= " The nBlocks field, if given, should be an integer.";
	throw new LuaInputException(errMsg);
    }
    lua_pop(L, 1);

    auto fsol = new FlowSolution(jobName, dir, tindx, nBlocks);
    flowSolutionStore ~= pushObj!(FlowSolution, FlowSolutionMT)(L, fsol);
    return 1;
} // end newFlowSolution()


void registerFlowSolution(lua_State* L)
{
    luaL_newmetatable(L, FlowSolutionMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");
    /* Register methods for use. */
    lua_pushcfunction(L, &newFlowSolution);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &toStringObj!(FlowSolution, FlowSolutionMT));
    lua_setfield(L, -2, "__tostring");
    // Make class visible
    lua_setglobal(L, FlowSolutionMT.toStringz);

    // [TODO]
    // lua_pushcfunction(L, &write_initial_flow_file_from_lua);
    // lua_setglobal(L, "write_initial_flow_file");
} // end registerFlowSolution()
