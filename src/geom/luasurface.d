/**
 * luasurface.d
 * Lua interface to ParametricSurface objects.
 *
 * Authors: Rowan G. and Peter J.
 * Version: 2015-02-24 first code
 *          2015-07-05 SubRangedSurface
 */

module luasurface;

import std.stdio;
import std.string;
import std.conv;
import std.uni;
import util.lua;
import util.lua_service;
import geom;
import gpath;
import surface;
import luageom;
import luagpath;

// Name of metatables -- these are the Lua access names.
immutable string CoonsPatchMT = "CoonsPatch";
immutable string AOPatchMT = "AOPatch";
immutable string ChannelPatchMT = "ChannelPatch";
immutable string SubRangedSurfaceMT = "SubRangedSurface";

static const(ParametricSurface)[] surfaceStore;

ParametricSurface checkSurface(lua_State* L, int index) {
    // We have to do a brute force test for each object
    // type in turn.
    if ( isObjType(L, index, CoonsPatchMT) )
	return checkObj!(CoonsPatch, CoonsPatchMT)(L, index);
    if ( isObjType(L, index, AOPatchMT ) )
	return checkObj!(AOPatch, AOPatchMT)(L, index);
    // if no match found then
    return null;
}

extern(C) int isSurface(lua_State* L)
{
    if ( checkSurface(L, 1) )
	lua_pushboolean(L, 1);
    else
	lua_pushboolean(L, 0);
    return 1;
}

extern(C) int opCallSurface(T, string MTname)(lua_State* L)
{
    auto surface = checkObj!(T, MTname)(L, 1);
    auto r = luaL_checknumber(L, 2);
    auto s = luaL_checknumber(L, 3);
    auto pt = surface(r, s);
    return pushVector3(L, pt);
}


// Helper functions for constructors.

string getPathFromTable(string pName)
{
    return `lua_getfield(L, index, "`~pName~`");
    paths["`~pName~`"] = checkPath(L, -1);
    if ( paths["`~pName~`"] is null ) luaL_error(L, toStringz(format(errMsgTmplt, ctorName, "`~pName~`")));
    lua_pop(L, 1);`;
}

void getPaths(lua_State *L, string ctorName, out Path[string] paths)
{
    // Assume that table is at index 1.
    string errMsgTmplt = `Error in call to %s:new.
The value set for the '%s' path was not of Path type.`; 
    int index = 1;
    mixin(getPathFromTable("north"));
    mixin(getPathFromTable("east"));
    mixin(getPathFromTable("south"));
    mixin(getPathFromTable("west"));
}

string getVecFromTable(string vName)
{
    return `lua_getfield(L, index, "`~vName~`");
auto `~vName~`Ptr = checkVector3(L, -1);
if ( `~vName~`Ptr is null ) luaL_error(L, toStringz(format(errMsgTmplt, ctorName, "`~vName~`")));
    else `~vName~` = *`~vName~`Ptr;
    lua_pop(L, 1);`;
}

void getVector3s(lua_State *L, string ctorName,
		 out Vector3 p00, out Vector3 p10, out Vector3 p11, out Vector3 p01)
{
    // Assume that table is at index 1.
    string errMsgTmplt = `Error in call to %s:new.
The value set for the '%s' corner was not of Vector3 type.`; 
    int index = 1;
    mixin(getVecFromTable("p00"));
    mixin(getVecFromTable("p10"));
    mixin(getVecFromTable("p11"));
    mixin(getVecFromTable("p01"));
}

void getRandS(lua_State* L, string ctorName,
	      out double r0, out double r1, out double s0, out double s1)
{
    // Table is at index 1
    string errMsgTmplt = format("Error in call to %s:new.", ctorName);
    errMsgTmplt ~= "The value for variable '%s' is not valid. It should be a number value.";
    r0 = getNumberFromTable(L, 1, "r0", false, 0.0, true, format(errMsgTmplt, "r0"));
    r1 = getNumberFromTable(L, 1, "r1", false, 1.0, true, format(errMsgTmplt, "r1"));
    s0 = getNumberFromTable(L, 1, "s0", false, 0.0, true, format(errMsgTmplt, "s0"));
    s1 = getNumberFromTable(L, 1, "s1", false, 1.0, true, format(errMsgTmplt, "s1"));
}

// --------------- Specific constructors for each Surface type -------------------

/**
 * This is the constructor for a CoonsPatch to be used from the Lua interface.
 *
 * At successful completion of this function, a new CoonsPatch object
 * is pushed onto the Lua stack.
 *
 * Supported constructions are:
 * -------------------------
 * patch0 = CoonsPatch:new{north=nPath, east=ePath, south=sPath, west=wPath}
 * patch2 = CoonsPatch:new{p00=a, p10=b, p11=c, p01=d}
 * --------------------------
 * Notes:
 * 1. See PJs diagram at top of geom.surface.d for ordering and labelling of
 *    paths and corners.
 * 2. No mix-n-match of constructors allowed. It is one of: 4 paths OR 4 corner points. If a path is found first, that constructor wins.
 *
 */

extern(C) int newCoonsPatch(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    
    if ( !lua_istable(L, 1) ) {
	string errMsg = `Error in constructor CoonPatch:new.
A table with input parameters is expected as the first argument.`;
	luaL_error(L, errMsg.toStringz);
    }
    // Look for a path. If found, proceed with construction from paths.
    lua_getfield(L, 1, "north");
    if ( ! lua_isnil(L, -1) ) { 
	lua_pop(L, 1);
	Path[string] paths;
	getPaths(L, "CoonsPatch", paths);
	auto cpatch = new CoonsPatch(paths["south"], paths["north"],
				     paths["west"], paths["east"]);
	surfaceStore ~= pushObj!(CoonsPatch, CoonsPatchMT)(L, cpatch);
	return 1;
    } else {
	lua_pop(L, 1);
    }
    // Instead, look for Vector3 objects.
    lua_getfield(L, 1, "p00");
    if ( ! lua_isnil(L, -1) ) {
	lua_pop(L, 1);
	Vector3 p00, p10, p11, p01;
	getVector3s(L, "CoonsPatch", p00, p10, p11, p01);
	writeln("p00= ", p00, " p10= ", p10, " p11= ", p11, " p01= ", p01);
	auto cpatch = new CoonsPatch(p00, p10, p11, p01);
	surfaceStore ~= pushObj!(CoonsPatch, CoonsPatchMT)(L, cpatch);
	return 1;
    }
    lua_pop(L, 1);
    // If we make it here, there's been an error in construction.
    string errMsg = `There's a problem in call to CoonsPatch.new.
Neither a list of named paths ('north', 'east', 'south', 'west')
nor a list of named corners ('p00', 'p10', 'p11', 'p01') were found.`;
    luaL_error(L, errMsg.toStringz);
    return 0;
} // end newCoonsPatch()


/**
 * This is constructor for an AOPatch object to be used from the Lua interface.
 *
 * At successful completion of this function, a new AOPatch object
 * is pushed onto the Lua stack.
 *
 * Supported constructions are:
 * -------------------------
 * patch0 = AOPatch:new{north=nPath, east=ePath, south=sPath, west=wPath}
 * patch2 = AOPatch:new{p00=a, p10=b, p11=c, p01=d}
 * --------------------------
 * Notes:
 * 1. See PJs diagram at top of geom.surface.d for ordering and labelling of
 *    paths and corners.
 * 2. No mix-n-match of constructors allowed. It is one of: 4 paths OR 4 corner points. If a path is found first, that constructor wins.
 *
 */
extern(C) int newAOPatch(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    
    if ( !lua_istable(L, 1) ) {
	string errMsg = `Error in constructor AOPatch:new.
A table with input parameters is expected as the first argument.`;
	luaL_error(L, errMsg.toStringz);
    }
    string errMsgTmplt = "Error in call to AOPatch:new.\n";
    errMsgTmplt ~= "A valid value for '%s' is not found in arguments.\n";
    errMsgTmplt ~= "The value, if present, should be a number.";
    int nx = to!int(getNumberFromTable(L, 1, "nx", false, 10.0, true, format(errMsgTmplt, "nx")));
    int ny = to!int(getNumberFromTable(L, 1, "ny", false, 10.0, true, format(errMsgTmplt, "ny")));
    // Look for a path. If found, proceed with construction from paths.
    lua_getfield(L, 1, "north");
    if ( ! lua_isnil(L, -1) ) { 
	lua_pop(L, 1);
	Path[string] paths;
	getPaths(L, "AOPatch", paths);
	auto aopatch = new AOPatch(paths["south"], paths["north"],
				   paths["west"], paths["east"],
				   nx, ny);
	surfaceStore ~= pushObj!(AOPatch, AOPatchMT)(L, aopatch);
	return 1;
    } else {
	lua_pop(L, 1);
    }
    // Instead, look for Vector3 objects.
    lua_getfield(L, 1, "p00");
    if ( ! lua_isnil(L, -1) ) {
	lua_pop(L, 1);
	Vector3 p00, p10, p11, p01;
	getVector3s(L, "AOPatch", p00, p10, p11, p01);
	writeln("p00= ", p00, " p10= ", p10, " p11= ", p11, " p01= ", p01);
	auto aopatch = new AOPatch(p00, p10, p11, p01, nx, ny);
	surfaceStore ~= pushObj!(AOPatch, AOPatchMT)(L, aopatch);
	return 1;
    }
    lua_pop(L, 1);
    // If we make it here, there's been an error in construction.
    string errMsg = `There's a problem in call to AOPatch.new.
Neither a list of named paths ('north', 'east', 'south', 'west')
nor a list of named corners ('p00', 'p10', 'p11', 'p01') were found.`;
    luaL_error(L, errMsg.toStringz);
    return 0;
} // end newAOPatch()


/**
 * This is the constructor for a ChannelPatch to be used from the Lua interface.
 *
 * At successful completion of this function, a new ChannelPatch object
 * is pushed onto the Lua stack.
 *
 * Supported constructions are:
 * -------------------------
 * patch0 = ChannelPatch:new{south=sPath, north=nPath}
 * patch1 = ChannelPatch:new{south=sPath, north=nPath, ruled=false, pure2D=false}
 * --------------------------
 */

extern(C) int newChannelPatch(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    
    if ( !lua_istable(L, 1) ) {
	string errMsg = `Error in constructor ChannelPatch:new.
A table with input parameters is expected as the first argument.`;
	luaL_error(L, errMsg.toStringz);
    }
    // Look for south and north paths.
    lua_getfield(L, 1, "south");
    auto south = checkPath(L, -1);
    if ( south is null ) {
	string errMsg = `Error in constructor ChannelPatch:new.
Couldn't find south Path.`;
	luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);
    lua_getfield(L, 1, "north");
    auto north = checkPath(L, -1);
    if ( north is null ) {
	string errMsg = `Error in constructor ChannelPatch:new.
Couldn't find north Path.`;
	luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);
    // Optional parameters
    string errMsgTmpltBool = `Error in call to ChannelPatch:new{}.
A valid value for '%s' was not found in list of arguments.
The value, if present, should be boolean (true or false).`;
    bool ruled = getBooleanFromTable(L, 1, "ruled", false, false, true, format(errMsgTmpltBool, "ruled"));
    bool pure2D = getBooleanFromTable(L, 1, "pure2D", false, false, true, format(errMsgTmpltBool, "pure2D"));
    // Construct the actual surface.
    auto cpatch = new ChannelPatch(south, north, ruled, pure2D);
    surfaceStore ~= pushObj!(ChannelPatch, ChannelPatchMT)(L, cpatch);
    return 1;
} // end newChannelPatch()


/**
 * This is constructor for a SubRangedSurface object
 * to be used from the Lua interface.
 *
 * At successful completion of this function, a new SubRangedSurface object
 * is pushed onto the Lua stack.
 *
 * Supported constructions are:
 * -------------------------
 * srs = SubRangedSurface:new{psurf, r0=0.0, r1=1.0, s0=0.0, s1=1.0}
 * --------------------------
 */
extern(C) int newSubRangedSurface(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    
    if ( !lua_istable(L, 1) ) {
	string errMsg = `Error in constructor SubRangedSurface:new.
A table with input parameters is expected as the first argument.`;
	luaL_error(L, errMsg.toStringz);
    }
    // Expect a Surface object at the first array position in the table.
    lua_rawgeti(L, 1, 1);
    if ( lua_isnil(L, -1) ) {
	string errMsg = `Error in call to SubRangedSurface:new{}. No table entry found.`;
	luaL_error(L, errMsg.toStringz());
    }
    auto psurf = checkSurface(L, -1);
    lua_pop(L, 1);
    if ( psurf is null ) {
	string errMsg = `Error in call to SubRangedSurface:new{}. No valid Surface object found.`;
	luaL_error(L, errMsg.toStringz());
    }
    double r0, r1, s0, s1;
    getRandS(L, "SubRangedSurface", r0, r1, s0, s1);
    auto srs = new SubRangedSurface(psurf, r0, r1, s0, s1);
    surfaceStore ~= pushObj!(SubRangedSurface, SubRangedSurfaceMT)(L, srs);
    return 1;
} // end newSubRangedSurface()

/* ---------- convenience functions -------------- */

string getPath(string path, string pos)
{
    return `lua_rawgeti(L, 1, `~pos~`);
auto `~path~` = checkPath(L, -1);
if ( `~path~` is null ) luaL_error(L, toStringz(format(errMsgTmpl, "`~path~`")));
lua_pop(L, 1);`;
}
extern(C) int makePatch(lua_State* L)
{
    if ( !lua_istable(L, 1) ) {
	string errMsg = `Error in call to makePatch.
A table is expected as the first argument. No table was found.`;
	luaL_error(L, errMsg.toStringz);
    }
    size_t n = lua_objlen(L, 1);
    if ( n != 4 ) {
	string errMsg = `Error in call to makePatch.
There should be four values listed first in the table.
These are Paths listed in the orded North, East, South, West.`;
	luaL_error(L, errMsg.toStringz);
    }
    string errMsgTmpl = `Error in call to makePatch.
The %s path is not a valid Path object.`;
    mixin(getPath("north", "1"));
    mixin(getPath("east", "2"));
    mixin(getPath("south", "3"));
    mixin(getPath("west", "4"));

    string gridType = "TFI";
    lua_getfield(L, 1, "gridType");
    if ( !lua_isnil(L, -1) ) {
	gridType = to!string(luaL_checkstring(L, -1));
    }
    lua_pop(L, 1);
    
    gridType = toUpper(gridType);
    if ( gridType == "AO" ) {
	auto patch = new AOPatch(south, north, west, east);
	surfaceStore ~= pushObj!(AOPatch, AOPatchMT)(L, patch);
	return 1;
    }
    // else
    auto patch = new CoonsPatch(south, north, west, east);
    surfaceStore ~= pushObj!(CoonsPatch, CoonsPatchMT)(L, patch);
    return 1;
} // end makePatch()

void registerSurfaces(lua_State* L)
{
    // Register the CoonsPatch object
    luaL_newmetatable(L, CoonsPatchMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newCoonsPatch);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallSurface!(CoonsPatch, CoonsPatchMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallSurface!(CoonsPatch, CoonsPatchMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(CoonsPatch, CoonsPatchMT));
    lua_setfield(L, -2, "__tostring");

    lua_setglobal(L, CoonsPatchMT.toStringz);

    // Register the AOPatch object
    luaL_newmetatable(L, AOPatchMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newAOPatch);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallSurface!(AOPatch, AOPatchMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallSurface!(AOPatch, AOPatchMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(AOPatch, AOPatchMT));
    lua_setfield(L, -2, "__tostring");

    lua_setglobal(L, AOPatchMT.toStringz);

    // Register the ChannelPatch object
    luaL_newmetatable(L, ChannelPatchMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newChannelPatch);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallSurface!(ChannelPatch, ChannelPatchMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallSurface!(ChannelPatch, ChannelPatchMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(ChannelPatch, ChannelPatchMT));
    lua_setfield(L, -2, "__tostring");

    lua_setglobal(L, ChannelPatchMT.toStringz);

    // Register the SubRangedSurface object
    luaL_newmetatable(L, SubRangedSurfaceMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newSubRangedSurface);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallSurface!(SubRangedSurface, SubRangedSurfaceMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallSurface!(SubRangedSurface, SubRangedSurfaceMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(SubRangedSurface, SubRangedSurfaceMT));
    lua_setfield(L, -2, "__tostring");

    lua_setglobal(L, SubRangedSurfaceMT.toStringz);

    // Register utility functions.
    lua_pushcfunction(L, &isSurface);
    lua_setglobal(L, "isSurface");
    lua_pushcfunction(L, &makePatch);
    lua_setglobal(L, "makePatch");
} // end registerSurfaces()
