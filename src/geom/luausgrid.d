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
import grid;
import sgrid;
import usgrid;
import luagrid;
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

extern(C) int get_vtx(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L);
    auto grid = checkObj!(T, MTname)(L, 1);
    size_t i = to!size_t(luaL_checkint(L, 2)); // Note that we expect 0 <= i < niv
    Vector3 vtx = grid.vertices[i];
    return pushVector3(L, vtx);
}

extern(C) int get_nfaces(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 1);
    lua_pushnumber(L, grid.nfaces);
    return 1;
}

extern(C) int get_nboundaries(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 1);
    lua_pushnumber(L, grid.nboundaries);
    return 1;
}

extern(C) int get_boundaryset_tag(lua_State* L)
{
    int narg = lua_gettop(L);
    auto grid = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 1);
    size_t i = to!size_t(luaL_checkint(L, 2)); // Note that we expect 0 <= i < nboundaries
    if (i < grid.boundaries.length) {
	lua_pushstring(L, grid.boundaries[i].tag.toStringz());
    } else {
	lua_pushnil(L);
    }
    return 1;
}

extern(C) int find_nearest_cell_centre_usg(lua_State *L)
{
    double x, y, z;
    int narg = lua_gettop(L); // assume narg == 2;
    auto grid = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 1);
    lua_getfield(L, 2, "x");
    if ( !lua_isnil(L, -1) ) {
	x = luaL_checknumber(L, -1);
    } else {
	x = 0.0;
    }
    lua_pop(L, 1);
    lua_getfield(L, 2, "y");
    if ( !lua_isnil(L, -1) ) {
	y = luaL_checknumber(L, -1);
    } else {
	y = 0.0;
    }
    lua_pop(L, 1);
    lua_getfield(L, 2, "z");
    if ( !lua_isnil(L, -1) ) {
	z = luaL_checknumber(L, -1);
    } else {
	z = 0.0;
    }
    lua_pop(L, 1);
    size_t indx;
    double dist;
    Vector3 p = Vector3(x, y, z);
    grid.find_nearest_cell_centre(p, indx, dist);
    lua_pushinteger(L, indx);
    lua_pushnumber(L, dist);
    return 2;
} // end find_nearest_cell_centre_usg()

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
    if (!checkAllowedNames(L, 1, ["sgrid","new_label","filename","fileName","scale","fmt",
				  "expect_gmsh_order_for_wedges", "label"])) {
	string errMsg = "Error in call to UnstructuredGrid:new{}. Invalid name in table.";
	luaL_error(L, errMsg.toStringz);
    }
    //
    UnstructuredGrid usgrid;
    string label = "";
    lua_getfield(L, 1, "label".toStringz);
    if (lua_isstring(L, -1)) {
	label = to!string(luaL_checkstring(L, -1));
    }
    lua_pop(L, 1); // dispose of label entry, even if it is a nil
    //
    // First, look for a StructuredGrid field.
    lua_getfield(L, 1, "sgrid".toStringz);
    if ( !lua_isnil(L, -1) ) {
	StructuredGrid sgrid = checkStructuredGrid(L, -1);
	if (!sgrid) {
	    string errMsg = "Error in UnstructuredGrid:new{}. sgrid not a StructuredGrid.";
	    luaL_error(L, errMsg.toStringz);
	}
	usgrid = new UnstructuredGrid(sgrid, label);
    }
    lua_pop(L, 1);
    //
    if (!usgrid) {
	// We didn't find a sgrid entry, so try for a file name.
	lua_getfield(L, 1, "filename".toStringz);
	if (lua_isnil(L, -1)) {
	    lua_pop(L, 1);
	    lua_getfield(L, 1, "fileName".toStringz);
	}
	if ( !lua_isnil(L, -1) ) {
	    string filename = to!string(luaL_checkstring(L, -1));
	    double scale = 1.0;
	    lua_getfield(L, 1, "scale".toStringz);
	    if ( !lua_isnil(L, -1) ) { scale = to!double(luaL_checknumber(L, -1)); }
	    lua_pop(L, 1); // dispose of scale item
	    string fmt = "gziptext";
	    lua_getfield(L, 1, "fmt".toStringz);
	    if ( !lua_isnil(L, -1) ) { fmt = to!string(luaL_checkstring(L, -1)); }
	    lua_pop(L, 1); // dispose of fmt item
	    bool expect_gmsh_order_for_wedges = true;
	    lua_getfield(L, 1, "expect_gmsh_order_for_wedges".toStringz);
	    if ( !lua_isnil(L, -1) ) {
		expect_gmsh_order_for_wedges = to!bool(lua_toboolean(L, -1));
	    }
	    lua_pop(L, 1); // dispose of expect_gmsh_order_for_wedges item
	    if (filename.length > 0) {
		usgrid = new UnstructuredGrid(filename, fmt, scale,
					      expect_gmsh_order_for_wedges,
					      label);
	    }
	} else {
	    string errMsg = "Error in UnstructuredGrid:new{}. expected a string for filename.";
	    luaL_error(L, errMsg.toStringz);
	}
	lua_pop(L, 1); // dispose of filename item
    }
    //
    if (usgrid) {
	unstructuredGridStore ~= pushObj!(UnstructuredGrid, UnstructuredGridMT)(L, usgrid);
	return 1;
    } else {
	string errMsg = "Error in UnstructuredGrid:new{}. Failed to construct a grid.";
	luaL_error(L, errMsg.toStringz);
	return 0;
    }
} // end newUnstructuredGrid()


extern(C) int usg_joinGrid(lua_State* L)
{
    int narg = lua_gettop(L);
    UnstructuredGrid master = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 1);
    if (!master) {
	string errMsg = "Error in UnstructuredGrid:joinGrid(). master grid not an UnstructuredGrid.";
	luaL_error(L, errMsg.toStringz);
    }
    UnstructuredGrid other = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 2);
    if (!other) {
	string errMsg = "Error in UnstructuredGrid:joinGrid(). other grid not an UnstructuredGrid.";
	luaL_error(L, errMsg.toStringz);
    }
    master.joinGrid(other);
    lua_settop(L, 1); // We leave the master grid at stack position 1
    return 1;
} // end usg_joinGrid()

extern(C) int usg_writeStats(lua_State* L)
{
    int narg = lua_gettop(L);
    UnstructuredGrid usgrid = checkObj!(UnstructuredGrid, UnstructuredGridMT)(L, 1);
    if (!usgrid) {
	string errMsg = "Error in UnstructuredGrid:writeStats(). grid not an UnstructuredGrid.";
	luaL_error(L, errMsg.toStringz);
    }
    usgrid.writeStats();
    lua_settop(L, 0);
    return 0;
} // end usg_writeStats()


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
    lua_pushcfunction(L, &get_type!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "get_type");
    lua_pushcfunction(L, &get_nvertices!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "get_nvertices");
    lua_pushcfunction(L, &get_ncells!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "get_ncells");
    lua_pushcfunction(L, &get_nfaces);
    lua_setfield(L, -2, "get_nfaces");
    lua_pushcfunction(L, &get_nboundaries);
    lua_setfield(L, -2, "get_nboundaries");
    lua_pushcfunction(L, &get_boundaryset_tag);
    lua_setfield(L, -2, "get_boundaryset_tag");
    lua_pushcfunction(L, &get_vtx!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "get_vtx");
    lua_pushcfunction(L, &cellVolume!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "cellVolume");
    lua_pushcfunction(L, &write_to_gzip_file!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "write_to_gzip_file");
    lua_pushcfunction(L, &write_to_vtk_file!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "write_to_vtk_file");
    lua_pushcfunction(L, &write_to_su2_file!(UnstructuredGrid, UnstructuredGridMT));
    lua_setfield(L, -2, "write_to_su2_file");
    lua_pushcfunction(L, &find_nearest_cell_centre_usg);
    lua_setfield(L, -2, "find_nearest_cell_centre");
    lua_pushcfunction(L, &usg_joinGrid);
    lua_setfield(L, -2, "joinGrid");
    lua_pushcfunction(L, &usg_writeStats);
    lua_setfield(L, -2, "writeStats");

    lua_setglobal(L, UnstructuredGridMT.toStringz);

} // end registerUnstructuredGrid()
    






