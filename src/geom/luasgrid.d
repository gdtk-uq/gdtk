/**
 * A Lua interface for the D sgrid (StructuredGrid) module.
 *
 * Authors: Rowan G. and Peter J.
 * Date: 2015-02-26
 */

module luasgrid;

// We cheat to get the C Lua headers by using LuaD.
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
import luageom;
import luasurface;
import luavolume;
import luaunifunction;

// Name of metatables
immutable string StructuredGridMT = "StructuredGrid";

static const(StructuredGrid)[] structuredGridStore;

StructuredGrid checkStructuredGrid(lua_State* L, int index) {
    if ( isObjType(L, index, StructuredGridMT) ) {
	return checkObj!(StructuredGrid, StructuredGridMT)(L, index);
    }
    // if all else fails
    return null;
}

extern(C) int copyStructuredGrid(T, string MTname)(lua_State* L)
{
    // Sometimes it's convenient to get a copy of a function.
    auto grid = checkObj!(T, MTname)(L, 1);
    structuredGridStore ~= pushObj!(T, MTname)(L, grid);
    return 1;
}

extern(C) int get_niv(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(T, MTname)(L, 1);
    lua_pushnumber(L, grid.niv);
    return 1;
}

extern(C) int get_njv(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(T, MTname)(L, 1);
    lua_pushnumber(L, grid.njv);
    return 1;
}

extern(C) int get_nkv(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(T, MTname)(L, 1);
    lua_pushnumber(L, grid.nkv);
    return 1;
}

extern(C) int get_vtx(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L);
    auto grid = checkObj!(T, MTname)(L, 1);
    size_t i = to!size_t(luaL_checkint(L, 2)); // Note that we expect 0 <= i < niv
    size_t j = to!size_t(luaL_checkint(L, 3));
    size_t k;
    if (narg > 3) { 
	k = to!size_t(luaL_checkint(L, 4));
    } else {
	k = 0; // Assume 2D grid
    }
    Vector3 vtx = grid[i,j,k];
    return pushVector3(L, vtx);
}

extern(C) int subgrid(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L);
    auto grid = checkObj!(T, MTname)(L, 1);
    size_t i0 = to!size_t(luaL_checkint(L, 2)); // Note that we expect 0 <= i0 < niv
    size_t ni = to!size_t(luaL_checkint(L, 3));
    size_t j0 = to!size_t(luaL_checkint(L, 4)); // Note that we expect 0 <= j0 < njv
    size_t nj = to!size_t(luaL_checkint(L, 5));
    size_t k0;
    if (narg > 5) { 
	k0 = to!size_t(luaL_checkint(L, 6));
    } else {
	k0 = 0; // Assume 2D grid
    }
    size_t nk;
    if (narg > 6) { 
	nk = to!size_t(luaL_checkint(L, 7));
    } else {
	nk = 1; // Assume 2D grid
    }
    auto new_subgrid = grid.subgrid(i0,ni,j0,nj,k0,nk);
    structuredGridStore ~= pushObj!(T, MTname)(L, new_subgrid);
    return 1;
} // end subgrid()

extern(C) int write_to_text_file(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 2;
    auto grid = checkObj!(T, MTname)(L, 1);
    auto fileName = to!string(luaL_checkstring(L, 2));
    grid.write_to_text_file(fileName, false);
    return 0;
}

/**
 * The Lua constructor for a StructuredGrid.
 *
 * Example construction in Lua:
 * grid = StructuredGrid:new{psurface=someParametricSurface, niv=10, njv=10,
 *                           cfList={cf_north, cf_east, cf_south, cf_west},
 *                           label="A-2D-Grid"}
 * grid3D = StructuredGrid:new{pvolume=someParametricVolume,
 *                             niv=11, njv=21, nkv=11,
 *                             cfList={},
 *                             label="A-3D-Grid"}
 */
extern(C) int newStructuredGrid(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = `Error in call to StructuredGrid:new{}.;
A table containing arguments is expected, but no table was found.`;
	luaL_error(L, errMsg.toStringz);
    }
    ParametricSurface psurface;
    ParametricVolume pvolume;
    int dimensions;
    // First, look for a ParametricSurface field, for a 2D grid.
    lua_getfield(L, 1, "psurface".toStringz);
    if ( !lua_isnil(L, -1) ) {
	dimensions = 2;
	psurface = checkSurface(L, -1);
	if (!psurface) {
	    string errMsg = "Error in StructuredGrid:new{}. psurface not a ParametricSurface.";
	    luaL_error(L, errMsg.toStringz);
	}
    }
    lua_pop(L, 1);
    // We didn't find a psurface entry, so try for a pvolume, for a 3D grid.
    if (!psurface) {
	lua_getfield(L, 1, "pvolume".toStringz);
	if ( !lua_isnil(L, -1) ) {
	    dimensions = 3;
	    pvolume = checkVolume(L, -1);
	    if (!pvolume) {
		string errMsg = "Error in StructuredGrid:new{}. pvolume not a ParametricVolume.";
		luaL_error(L, errMsg.toStringz);
	    }
	} else {
	    string errMsg = "Error in StructuredGrid:new{}. neither psurface nor pvolume found.";
	    luaL_error(L, errMsg.toStringz);
	}
	lua_pop(L, 1);
    }

    // Get clustering functions, 
    // filling in nil or invalid entries with LinearFunction.
    UnivariateFunction[] cfList;
    int number_of_edges = (dimensions == 2) ? 4 : 12;
    lua_getfield(L, 1, "cfList".toStringz);
    if ( lua_istable(L, -1) ) {
	// Extract the cluster functions from the table at top of stack. 
	foreach (i; 0 .. number_of_edges) {
	    lua_rawgeti(L, -1, i+1);
	    if ( lua_isnil(L, -1) ) {
		cfList ~= new LinearFunction();
	    } else {
		auto mycf = checkUnivariateFunction(L, -1);
		if ( mycf ) {
		    cfList ~= mycf;
		} else {
		    cfList ~= new LinearFunction();
		}
	    }
	    lua_pop(L, 1);
	}
    } else {
	// Didn't find a table of cluster functions.
	foreach (i; 0 .. number_of_edges) {
	    cfList ~= new LinearFunction();
	}
    }
    lua_pop(L, 1);

    string errMsgTmplt = `Error in StructuredGrid:new{}.
A valid value for '%s' was not found in list of arguments.
The value, if present, should be a number.`;
    int niv = getIntegerFromTable(L, 1, "niv", true, 0, true, format(errMsgTmplt, "niv"));
    int njv = getIntegerFromTable(L, 1, "njv", true, 0, true, format(errMsgTmplt, "njv"));
    int nkv = 1;
    if ( dimensions == 3 ) {
	nkv = getIntegerFromTable(L, 1, "nkv", true, 0, true, format(errMsgTmplt, "nkv"));
    }

    StructuredGrid grid;
    if ( dimensions == 2 ) {
	grid = new StructuredGrid(psurface, niv, njv, cfList);
    } else {
	grid = new StructuredGrid(pvolume, niv, njv, nkv, cfList);
    }
    structuredGridStore ~= pushObj!(StructuredGrid, StructuredGridMT)(L, grid);
    return 1;
} // end newStructuredGrid()


extern(C) int importGridproGrid(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg == 0 ) {
	string errMsg = `Error in call to importGridproGrid().
At least one argument is required: the name of the Gridpro file.`;
	luaL_error(L, errMsg.toStringz);
    }
    auto fname = to!string(luaL_checkstring(L, 1));
    double scale = 1.0;
    if ( narg >= 2 ) {
	scale = luaL_checknumber(L, 2);
    }
    auto sgrids = import_gridpro_grid(fname, scale);
    lua_newtable(L);
    foreach ( int i, grid; sgrids ) {
	structuredGridStore ~= pushObj!(StructuredGrid, StructuredGridMT)(L, grid);
	lua_rawseti(L, -2, i+1);
    }
    return 1;
}


void registerStructuredGrid(lua_State* L)
{
    // Register the StructuredGrid object
    luaL_newmetatable(L, StructuredGridMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newStructuredGrid);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &toStringObj!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyStructuredGrid!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &get_niv!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "get_niv");
    lua_pushcfunction(L, &get_njv!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "get_njv");
    lua_pushcfunction(L, &get_nkv!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "get_nkv");
    lua_pushcfunction(L, &get_vtx!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "get_vtx");
    lua_pushcfunction(L, &subgrid!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "subgrid");
    lua_pushcfunction(L, &write_to_text_file!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "write_to_text_file");

    lua_setglobal(L, StructuredGridMT.toStringz);

    // Global functions available for use
    lua_pushcfunction(L, &importGridproGrid);
    lua_setglobal(L, "importGridproGrid");

} // end registerStructuredGrid()
    






