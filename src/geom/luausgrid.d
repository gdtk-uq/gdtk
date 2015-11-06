/**
 * A Lua interface for the D usgrid (UnstructuredGrid) module.
 *
 * Authors: Peter J. and Rowan G.
 * Date: 2015-11-06 adapted from luasgrid.d
 */

module luausgrid;

import std.stdio;
import std.string;
import std.conv;
import util.lua;
import util.lua_service;
import univariatefunctions;
import geom;
import gpath;
import surface;
import volume;
import sgrid;
import usgrid;
import luasgrid;
import luageom;
// Might put the following back later when we generate our own 
// UnstructuredGrid objects from parametric surfaces and volumes.
// import luasurface;
// import luavolume;
// import luaunifunction;

// Name of metatables
immutable string UnstructuredGridMT = "UnstructuredGrid";

static const(UnstructuredGrid)[] unstructuredGridStore;

UnstructuredGrid checkUnstructuredGrid(lua_State* L, int index) {
    if ( isObjType(L, index, UnstructuredGridMT) ) {
	return checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, index);
    }
    // all else fails
    return null;
}

extern(C) int copyUnstructuredGrid(T, string MTname)(lua_State* L)
{
    // Sometimes it's convenient to get a copy of a grid.
    auto grid = checkObj!(T, MTname)(L, 1);
    unstructuredGridStore ~= pushObj!(T, MTname)(L, grid);
    return 1;
}

extern(C) int get_dimensions(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(T, MTname)(L, 1);
    lua_pushnumber(L, grid.dimensions);
    return 1;
}

extern(C) int get_nvertices(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(T, MTname)(L, 1);
    lua_pushnumber(L, grid.nvertices);
    return 1;
}

extern(C) int get_ncells(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(T, MTname)(L, 1);
    lua_pushnumber(L, grid.ncells);
    return 1;
}

extern(C) int get_nfaces(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(T, MTname)(L, 1);
    lua_pushnumber(L, grid.nfaces);
    return 1;
}

extern(C) int get_vtx(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L);
    auto grid = checkObj!(T, MTname)(L, 1);
    size_t i = to!size_t(luaL_checkint(L, 2)); // Note that we expect 0 <= i < niv
    Vector3 vtx = grid.vertices[i];
    return pushVector3(L, vtx);
}

extern(C) int write_to_gzip_file(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 2;
    auto grid = checkObj!(T, MTname)(L, 1);
    auto fileName = to!string(luaL_checkstring(L, 2));
    grid.write_to_gzip_file(fileName);
    return 0;
}

extern(C) int write_to_vtk_file(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 2;
    auto grid = checkObj!(T, MTname)(L, 1);
    auto fileName = to!string(luaL_checkstring(L, 2));
    grid.write_to_vtk_file(fileName);
    return 0;
}

/**
 * The Lua constructor for a UnstructuredGrid.
 *
 * Example construction in Lua:
 * grid = UnstructuredGrid:new{sgrid=sgrid_object, new_label="a-new-label"}
 * grid = UnstructuredGrid:new{filename="a-file-name", new_label="a-new-label"}
 */
extern(C) int newUnstructuredGrid(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = `Error in call to UnstructuredGrid:new{}.;
A table containing arguments is expected, but no table was found.`;
	luaL_error(L, errMsg.toStringz);
    }
    UnstructuredGrid usgrid;
    StructuredGrid sgrid;
    // First, look for a StructuredGrid field.
    lua_getfield(L, 1, "sgrid".toStringz);
    if ( !lua_isnil(L, -1) ) {
	sgrid = checkStructuredGrid(L, -1);
	if (!sgrid) {
	    string errMsg = "Error in UnstructuredGrid:new{}. sgrid not a StructuredGrid.";
	    luaL_error(L, errMsg.toStringz);
	}
	usgrid = new UnstructuredGrid(sgrid); // [TODO] new_label
    }
    lua_pop(L, 1);

    if (!usgrid) {
	// We didn't find a sgrid entry, so try for a file name.
	lua_getfield(L, 1, "filename".toStringz);
	if ( !lua_isnil(L, -1) ) {
	    string filename = to!string(luaL_checkstring(L, -1));
	    if (filename.length > 0) {
		usgrid = new UnstructuredGrid(filename, "gziptext"); // [TODO] new_label
	    }
	} else {
	    string errMsg = "Error in UnstructuredGrid:new{}. expected a string for filename.";
	    luaL_error(L, errMsg.toStringz);
	}
	lua_pop(L, 1);
    }

    if (usgrid) {
	unstructuredGridStore ~= pushObj!(UnstructuredGrid, UnstructuredGridMT)(L, usgrid);
	return 1;
    } else {
	string errMsg = "Error in UnstructuredGrid:new{}. Failed to construct a grid.";
	luaL_error(L, errMsg.toStringz);
	return 0;
    }
} // end newUnstructuredGrid()


void registerUnstructuredGrid(lua_State* L)
{
    // Register the UnstructuredGrid object
    luaL_newmetatable(L, UnstructuredGridMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newUnstructuredGrid);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &toStringObj!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyUnstructuredGrid!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &get_dimensions!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "get_dimensions");
    lua_pushcfunction(L, &get_nvertices!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "get_nvertices");
    lua_pushcfunction(L, &get_ncells!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "get_ncells");
    lua_pushcfunction(L, &get_nfaces!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "get_nfaces");
    lua_pushcfunction(L, &get_vtx!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "get_vtx");
    lua_pushcfunction(L, &write_to_gzip_file!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "write_to_gzip_file");
    lua_pushcfunction(L, &write_to_vtk_file!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "write_to_vtk_file");

    lua_setglobal(L, UnstructuredGridMT.toStringz);

} // end registerUnstructuredGrid()
    






