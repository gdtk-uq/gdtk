/**
 * luaflowsolution.d
 * Lua interface to the FlowSolution object for use in postprocessing.
 *
 * Authors: Peter J and Rowan G.
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
import geom;
import luageom;
import sgrid;
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
    lua_remove(L, 1); // Remove first argument "this".
    if ( !lua_istable(L, 1) ) {
	string errMsg = "Error in call to FlowSolution:new.";
	errMsg ~= " A table is expected as first (and only) argument.";
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


extern(C) int find_nearest_cell_centre_from_lua(lua_State* L)
{
    auto fsol = checkFlowSolution(L, 1);
    if ( !lua_istable(L, 2) ) {
	string errMsg = "Error in call to FlowSolution:find_nearest_cell_centre.";
	errMsg ~= " A table is expected as first (and only) argument to the method.";
	luaL_error(L, errMsg.toStringz);
    }

    double x;
    lua_getfield(L, 2, "x");
    if ( lua_isnil(L, -1) ) {
	x = 0.0;
    } else if ( lua_isnumber(L, -1) ) {
	x = to!double(luaL_checknumber(L, -1));
    } else {
	string errMsg = "Error in call to FlowSolution:find_nearest_cell_centre.";
	errMsg ~= " A field for x was found, but the content was not valid.";
	errMsg ~= " The x field, if given, should be a double.";
	throw new LuaInputException(errMsg);
    }
    lua_pop(L, 1);

    double y;
    lua_getfield(L, 2, "y");
    if ( lua_isnil(L, -1) ) {
	y = 0.0;
    } else if ( lua_isnumber(L, -1) ) {
	y = to!double(luaL_checknumber(L, -1));
    } else {
	string errMsg = "Error in call to FlowSolution:find_nearest_cell_centre.";
	errMsg ~= " A field for y was found, but the content was not valid.";
	errMsg ~= " The y field, if given, should be a double.";
	throw new LuaInputException(errMsg);
    }
    lua_pop(L, 1);

    double z;
    lua_getfield(L, 2, "z");
    if ( lua_isnil(L, -1) ) {
	z = 0.0;
    } else if ( lua_isnumber(L, -1) ) {
	z = to!double(luaL_checknumber(L, -1));
    } else {
	string errMsg = "Error in call to FlowSolution:find_nearest_cell_centre.";
	errMsg ~= " A field for z was found, but the content was not valid.";
	errMsg ~= " The z field, if given, should be a double.";
	throw new LuaInputException(errMsg);
    }
    lua_pop(L, 1);

    auto indices = fsol.find_nearest_cell_centre(x, y, z);

    lua_settop(L, 0); // clear stack
    lua_newtable(L); // anonymous table { }
    auto tblIndx = lua_gettop(L);
    lua_pushnumber(L, indices[0]);
    lua_setfield(L, tblIndx, "ib");
    lua_pushnumber(L, indices[1]);
    lua_setfield(L, tblIndx, "i");
    lua_pushnumber(L, indices[2]);
    lua_setfield(L, tblIndx, "j");
    lua_pushnumber(L, indices[3]);
    lua_setfield(L, tblIndx, "k");
    return 1; // Just the table of indices is left on the stack.
} // end find_nearest_cell_centre_from_lua()


extern(C) int get_nic(lua_State* L)
{
    auto fsol = checkFlowSolution(L, 1);
    int ib = to!int(luaL_checknumber(L, 2));
    lua_settop(L, 0);
    lua_pushnumber(L, fsol.flowBlocks[ib].nic);
    return 1;
}

extern(C) int get_njc(lua_State* L)
{
    auto fsol = checkFlowSolution(L, 1);
    int ib = to!int(luaL_checknumber(L, 2));
    lua_settop(L, 0);
    lua_pushnumber(L, fsol.flowBlocks[ib].njc);
    return 1;
}

extern(C) int get_nkc(lua_State* L)
{
    auto fsol = checkFlowSolution(L, 1);
    int ib = to!int(luaL_checknumber(L, 2));
    lua_settop(L, 0);
    lua_pushnumber(L, fsol.flowBlocks[ib].nkc);
    return 1;
}


extern(C) int get_vtx_pos(lua_State* L)
{
    auto fsol = checkFlowSolution(L, 1);
    if ( !lua_istable(L, 2) ) {
	string errMsg = "Error in call to FlowSolution:get_vtx_pos.";
	errMsg ~= " A table is expected as first (and only) argument to the method.";
	luaL_error(L, errMsg.toStringz);
    }

    int ib;
    lua_getfield(L, 2, "ib");
    if ( lua_isnil(L, -1) ) {
	ib = 0;
    } else if ( lua_isnumber(L, -1) ) {
	ib = to!int(luaL_checknumber(L, -1));
    } else {
	string errMsg = "Error in call to FlowSolution:get_vtx_pos.";
	errMsg ~= " A field for ib was found, but the content was not valid.";
	errMsg ~= " The ib field, if given, should be an integer.";
	throw new LuaInputException(errMsg);
    }
    lua_pop(L, 1);

    int i;
    lua_getfield(L, 2, "i");
    if ( lua_isnil(L, -1) ) {
	i = 0;
    } else if ( lua_isnumber(L, -1) ) {
	i = to!int(luaL_checknumber(L, -1));
    } else {
	string errMsg = "Error in call to FlowSolution:get_vtx_pos.";
	errMsg ~= " A field for i was found, but the content was not valid.";
	errMsg ~= " The i field, if given, should be an integer.";
	throw new LuaInputException(errMsg);
    }
    lua_pop(L, 1);

    int j;
    lua_getfield(L, 2, "j");
    if ( lua_isnil(L, -1) ) {
	j = 0;
    } else if ( lua_isnumber(L, -1) ) {
	j = to!int(luaL_checknumber(L, -1));
    } else {
	string errMsg = "Error in call to FlowSolution:get_vtx_pos.";
	errMsg ~= " A field for j was found, but the content was not valid.";
	errMsg ~= " The j field, if given, should be an integer.";
	throw new LuaInputException(errMsg);
    }
    lua_pop(L, 1);

    int k;
    lua_getfield(L, 2, "k");
    if ( lua_isnil(L, -1) ) {
	k = 0;
    } else if ( lua_isnumber(L, -1) ) {
	k = to!int(luaL_checknumber(L, -1));
    } else {
	string errMsg = "Error in call to FlowSolution:get_vtx_pos.";
	errMsg ~= " A field for k was found, but the content was not valid.";
	errMsg ~= " The k field, if given, should be an integer.";
	throw new LuaInputException(errMsg);
    }
    lua_pop(L, 1);
    Vector3 vtx = fsol.gridBlocks[ib][i, j, k];
    lua_settop(L, 0);
    return pushVector3(L, vtx);
} // end get_vtx_pos()


extern(C) int get_var_names(lua_State* L)
{
    auto fsol = checkFlowSolution(L, 1);
    lua_settop(L, 0); // clear stack
    lua_newtable(L); // anonymous table { }
    auto tblIndx = lua_gettop(L);
    foreach (ivar; 0 .. fsol.flowBlocks[0].variableNames.length) {
	lua_pushstring(L, fsol.flowBlocks[0].variableNames[ivar].toStringz);
	lua_rawseti(L, tblIndx, to!int(ivar+1)); // Lua table indexing starts from 1
    }
    return 1; // Just the table of indices is left on the stack.
} // end get_var_names()


extern(C) int get_cell_data(lua_State* L)
{
    auto fsol = checkFlowSolution(L, 1);
    if ( !lua_istable(L, 2) ) {
	string errMsg = "Error in call to FlowSolution:get_cell_data.";
	errMsg ~= " A table is expected as first (and only) argument to the method.";
	luaL_error(L, errMsg.toStringz);
    }

    int ib;
    lua_getfield(L, 2, "ib");
    if ( lua_isnil(L, -1) ) {
	ib = 0;
    } else if ( lua_isnumber(L, -1) ) {
	ib = to!int(luaL_checknumber(L, -1));
    } else {
	string errMsg = "Error in call to FlowSolution:get_cell_data.";
	errMsg ~= " A field for ib was found, but the content was not valid.";
	errMsg ~= " The ib field, if given, should be an integer.";
	throw new LuaInputException(errMsg);
    }
    lua_pop(L, 1);

    int i;
    lua_getfield(L, 2, "i");
    if ( lua_isnil(L, -1) ) {
	i = 0;
    } else if ( lua_isnumber(L, -1) ) {
	i = to!int(luaL_checknumber(L, -1));
    } else {
	string errMsg = "Error in call to FlowSolution:get_cell_data.";
	errMsg ~= " A field for i was found, but the content was not valid.";
	errMsg ~= " The i field, if given, should be an integer.";
	throw new LuaInputException(errMsg);
    }
    lua_pop(L, 1);

    int j;
    lua_getfield(L, 2, "j");
    if ( lua_isnil(L, -1) ) {
	j = 0;
    } else if ( lua_isnumber(L, -1) ) {
	j = to!int(luaL_checknumber(L, -1));
    } else {
	string errMsg = "Error in call to FlowSolution:get_cell_data.";
	errMsg ~= " A field for j was found, but the content was not valid.";
	errMsg ~= " The j field, if given, should be an integer.";
	throw new LuaInputException(errMsg);
    }
    lua_pop(L, 1);

    int k;
    lua_getfield(L, 2, "k");
    if ( lua_isnil(L, -1) ) {
	k = 0;
    } else if ( lua_isnumber(L, -1) ) {
	k = to!int(luaL_checknumber(L, -1));
    } else {
	string errMsg = "Error in call to FlowSolution:get_cell_data.";
	errMsg ~= " A field for k was found, but the content was not valid.";
	errMsg ~= " The k field, if given, should be an integer.";
	throw new LuaInputException(errMsg);
    }
    lua_pop(L, 1);

    lua_settop(L, 0); // clear stack
    lua_newtable(L); // anonymous table { }
    auto tblIndx = lua_gettop(L);
    foreach (varName; fsol.flowBlocks[ib].variableNames) {
	lua_pushnumber(L, fsol.flowBlocks[ib][varName,i,j,k]);
	lua_setfield(L, tblIndx, varName.toStringz);
    }
    return 1; // Just the table of indices is left on the stack.
} // end get_cell_data()


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
    lua_pushcfunction(L, &find_nearest_cell_centre_from_lua);
    lua_setfield(L, -2, "find_nearest_cell_centre");
    lua_pushcfunction(L, &get_nic);
    lua_setfield(L, -2, "get_nic");
    lua_pushcfunction(L, &get_njc);
    lua_setfield(L, -2, "get_njc");
    lua_pushcfunction(L, &get_nkc);
    lua_setfield(L, -2, "get_nkc");
    lua_pushcfunction(L, &get_vtx_pos);
    lua_setfield(L, -2, "get_vtx");
    lua_pushcfunction(L, &get_var_names);
    lua_setfield(L, -2, "get_var_names");
    lua_pushcfunction(L, &get_cell_data);
    lua_setfield(L, -2, "get_cell_data");
    // Make class visible
    lua_setglobal(L, FlowSolutionMT.toStringz);
} // end registerFlowSolution()
