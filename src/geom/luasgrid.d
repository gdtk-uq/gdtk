/**
 * A Lua interface for the D sgrid (StructuredGrid) module.
 *
 * Authors: Rowan G. and Peter J.
 * Date: 2015-02-26
 */

module luasgrid;

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
import luageom;
import luagpath;
import luasurface;
import luavolume;
import luaunifunction;
import luagrid;

// Name of metatables
immutable string StructuredGridMT = "StructuredGrid";

static const(StructuredGrid)[] structuredGridStore;

StructuredGrid checkStructuredGrid(lua_State* L, int index) {
    if ( isObjType(L, index, StructuredGridMT) ) {
	return checkObj!(StructuredGrid, StructuredGridMT)(L, index);
    }
    return null; // on fail
}

extern(C) int copyStructuredGrid(T, string MTname)(lua_State* L)
{
    // Sometimes it's convenient to get a copy of a grid.
    auto grid = checkObj!(T, MTname)(L, 1);
    structuredGridStore ~= pushObj!(T, MTname)(L, grid);
    return 1;
}

extern(C) int get_niv(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    lua_pushnumber(L, grid.niv);
    return 1;
}

extern(C) int get_njv(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    lua_pushnumber(L, grid.njv);
    return 1;
}

extern(C) int get_nkv(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 1; This is a getter
    auto grid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    lua_pushnumber(L, grid.nkv);
    return 1;
}

extern(C) int get_vtx(T, string MTname)(lua_State* L)
{
    int narg = lua_gettop(L);
    auto grid = checkObj!(T, MTname)(L, 1);
    size_t i = to!size_t(luaL_checkint(L, 2));
    // Note that we expect 0 <= i < niv
    size_t j = 0; if (narg > 2) { j = to!size_t(luaL_checkint(L, 3)); }
    size_t k = 0; if (narg > 3) { k = to!size_t(luaL_checkint(L, 4)); }
    Vector3* vtx = grid[i,j,k];
    return pushVector3(L, *vtx);
}

extern(C) int subgrid(lua_State* L)
{
    int narg = lua_gettop(L);
    auto grid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    // Note that we expect 0 <= i0 < niv
    size_t i0 = to!size_t(luaL_checkint(L, 2));
    size_t ni = to!size_t(luaL_checkint(L, 3));
    size_t j0 = 0; if (narg > 3) { j0 = to!size_t(luaL_checkint(L, 4)); }
    size_t nj = 1; if (narg > 4) { nj = to!size_t(luaL_checkint(L, 5)); }
    size_t k0 = 0; if (narg > 5) { k0 = to!size_t(luaL_checkint(L, 6)); }
    size_t nk = 1; if (narg > 6) { nk = to!size_t(luaL_checkint(L, 7)); }
    auto new_subgrid = grid.subgrid(i0,ni,j0,nj,k0,nk);
    structuredGridStore ~= pushObj!(StructuredGrid, StructuredGridMT)(L, new_subgrid);
    return 1;
} // end subgrid()

extern(C) int joinGrid(lua_State* L)
{
    int narg = lua_gettop(L); // assume narg == 3;
    auto selfGrid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
    auto otherGrid = checkObj!(StructuredGrid, StructuredGridMT)(L, 2);
    auto joinLocation = to!string(luaL_checkstring(L, 3));
    selfGrid.joinGrid(otherGrid, joinLocation);
    return 0;
}

extern(C) int find_nearest_cell_centre_sg(lua_State *L)
{
    double x, y, z;
    int narg = lua_gettop(L); // assume narg == 2;
    auto grid = checkObj!(StructuredGrid, StructuredGridMT)(L, 1);
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
} // end find_nearest_cell_centre_sg()

/**
 * The Lua constructor for a StructuredGrid.
 *
 * Example construction in Lua:
 * grid1D = StructuredGrid:new{path=somePath, niv=10, cf=my_cf01}
 * grid2D = StructuredGrid:new{psurface=someParametricSurface, niv=10, njv=10,
 *                             cfList={north=cfn, east=cfe, south=cfs, west=cfw},
 *                             label="A-2D-Grid"}
 * grid3D = StructuredGrid:new{pvolume=someParametricVolume,
 *                             niv=11, njv=21, nkv=11,
 *                             cfList={edge01=cf01, edge32=cf32, edge45=cf45, edge76=cf76},
 *                             label="A-3D-Grid"}
 */
extern(C) int newStructuredGrid(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = "Error in call to StructuredGrid:new{}.; " ~
	    "A table containing arguments is expected, but no table was found.";
	luaL_error(L, errMsg.toStringz);
    }
    Path mypath;
    ParametricSurface psurface;
    ParametricVolume pvolume;
    int dimensions = 0;
    // First, look for a Path field, for a 1D grid.
    lua_getfield(L, 1, "path".toStringz);
    if ( !lua_isnil(L, -1) ) {
	dimensions = 1;
	mypath = checkPath(L, -1);
	if (!mypath) {
	    string errMsg = "Error in StructuredGrid:new{}. path not a Path.";
	    luaL_error(L, errMsg.toStringz);
	}
    }
    lua_pop(L, 1);
    if (dimensions == 0) {
	// Next, look for a ParametricSurface field, for a 2D grid.
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
    }
    if (dimensions == 0) {
	// No path or surface, so try for a pvolume, for a 3D grid.
	lua_getfield(L, 1, "pvolume".toStringz);
	if ( !lua_isnil(L, -1) ) {
	    dimensions = 3;
	    pvolume = checkVolume(L, -1);
	    if (!pvolume) {
		string errMsg = "Error in StructuredGrid:new{}. pvolume not a ParametricVolume.";
		luaL_error(L, errMsg.toStringz);
	    }
	}
	lua_pop(L, 1);
    }
    if (dimensions == 0) {
	string errMsg = "Error in StructuredGrid:new{}. no path, psurface or pvolume found.";
	luaL_error(L, errMsg.toStringz);
    }
 
    // Get clustering functions, filling in nil or invalid entries with LinearFunction.
    UnivariateFunction cf;
    UnivariateFunction[] cfList;
    if (dimensions == 1) {
	lua_getfield(L, 1, "cf".toStringz);
	cf = checkUnivariateFunction(L, -1);
	if (!cf) { cf = new LinearFunction(); }
	lua_pop(L, 1);
    } else {
	// There will be a list of cluster functions for 2D or 3D.
	int number_of_edges = (dimensions == 2) ? 4 : 12;
	lua_getfield(L, 1, "cfList".toStringz);
	if ( lua_istable(L, -1) ) {
	    // Before going on, check the caller gave us named parameters.
	    if ( lua_objlen(L, -1) != 0 ) {
		// It appears that the caller has tried to set arguments as an array
		string errMsg = "Error in StructuredGrid:new{}. " ~
		    "A table of named parameters is expected for the cfList. " ~
		    "An array style of parameters was found.";
		luaL_error(L, errMsg.toStringz);
	    }

	    string[] edges;
	    // Extract the cluster functions from the table at top of stack.
	    if (dimensions == 2) {
		edges = ["north", "east", "south", "west"];
	    } else {
		edges = ["edge01", "edge12", "edge32", "edge03",
			 "edge45", "edge56", "edge76", "edge47",
			 "edge04", "edge15", "edge26", "edge37"];
	    }
	    foreach (edge_name; edges) {
		lua_getfield(L, -1, edge_name.toStringz);
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
    }

    string errMsgTmplt = "Error in StructuredGrid:new{}. " ~
	"A valid value for '%s' was not found in list of arguments. " ~
	"The value, if present, should be a number.";
    int niv = getIntegerFromTable(L, 1, "niv", true, 0, true, format(errMsgTmplt, "niv"));
    int njv = 1;
    if (dimensions > 1) {
	njv = getIntegerFromTable(L, 1, "njv", true, 0, true, format(errMsgTmplt, "njv"));
    }
    int nkv = 1;
    if (dimensions == 3) {
	nkv = getIntegerFromTable(L, 1, "nkv", true, 0, true, format(errMsgTmplt, "nkv"));
    }

    StructuredGrid grid;
    switch (dimensions) {
    case 1: grid = new StructuredGrid(mypath, niv, cf); break;
    case 2: grid = new StructuredGrid(psurface, niv, njv, cfList); break;
    case 3: grid = new StructuredGrid(pvolume, niv, njv, nkv, cfList); break;
    default: assert(0);
    }
    structuredGridStore ~= pushObj!(StructuredGrid, StructuredGridMT)(L, grid);
    return 1;
} // end newStructuredGrid()


extern(C) int importGridproGrid(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg == 0 ) {
	string errMsg = "Error in call to importGridproGrid(). " ~
	    "At least one argument is required: the name of the Gridpro file.";
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
} // end importGridproGrid()

extern(C) int writeGridsAsPlot3D(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( narg != 3 ) {
	string errMsg = "Error in call to writeGridsAsPlot3D(). " ~
	    "3 arguments are required: name of file, array of grids, and dimension/";
	luaL_error(L, errMsg.toStringz);
    }
    auto fname = to!string(luaL_checkstring(L, 1));
    if (!lua_istable(L, 2)) {
	string errMsg = "Error in call to writeGridsAsPlot3D(). " ~
	    "The second argument should be an array of grids.";
	luaL_error(L, errMsg.toStringz);
    }
    int n = to!int(lua_objlen(L, 2));
    StructuredGrid[] grids;
    foreach (i; 1 .. n+1) {
	lua_rawgeti(L, 2, i);
	grids ~= checkObj!(StructuredGrid, StructuredGridMT)(L, -1);
	lua_pop(L, 1);
    }
    int dim = to!int(luaL_checkinteger(L, 3));
    sgrid.writeGridsAsPlot3D(fname, grids, dim);
    return 0;
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
    lua_pushcfunction(L, &get_ncells!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "get_ncells");
    lua_pushcfunction(L, &get_niv);
    lua_setfield(L, -2, "get_niv");
    lua_pushcfunction(L, &get_njv);
    lua_setfield(L, -2, "get_njv");
    lua_pushcfunction(L, &get_nkv);
    lua_setfield(L, -2, "get_nkv");
    lua_pushcfunction(L, &get_vtx!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "get_vtx");
    lua_pushcfunction(L, &subgrid);
    lua_setfield(L, -2, "subgrid");
    lua_pushcfunction(L, &cellVolume!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "cellVolume");
    lua_pushcfunction(L, &write_to_vtk_file!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "write_to_vtk_file");
    lua_pushcfunction(L, &write_to_gzip_file!(StructuredGrid, StructuredGridMT));
    lua_setfield(L, -2, "write_to_gzip_file");
    lua_pushcfunction(L, &joinGrid);
    lua_setfield(L, -2, "joinGrid");
    lua_pushcfunction(L, &find_nearest_cell_centre_sg);
    lua_setfield(L, -2, "find_nearest_cell_centre");

    lua_setglobal(L, StructuredGridMT.toStringz);

    // Global functions available for use
    lua_pushcfunction(L, &importGridproGrid);
    lua_setglobal(L, "importGridproGrid");
    lua_pushcfunction(L, &writeGridsAsPlot3D);
    lua_setglobal(L, "writeGridsAsPlot3D");


} // end registerStructuredGrid()
    






