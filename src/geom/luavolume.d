/**
 * luavolume.d
 * Lua interface to ParametricVolume objects.
 *
 * Authors: Peter J. and Rowan G.
 * Version: 2015-04-07, following Rowan's code in luasurface.d
 *          2015-07-05, introduced SubRangedVolume
 */

module luavolume;

import std.stdio;
import std.string;
import std.conv;
import std.uni;
import util.lua;
import util.lua_service;
import geom;
import gpath;
import surface;
import volume;
import luageom;
import luagpath;
import luasurface;

/// Name of the metatables -- these are the Lua access name.
immutable string TFIVolumeMT = "TFIVolume";
immutable string SubRangedVolumeMT = "SubRangedVolume";

static const(ParametricVolume)[] volumeStore;

ParametricVolume checkVolume(lua_State* L, int index) {
    // We have to do a brute force test for each object
    // type in turn.
    if ( isObjType(L, index, TFIVolumeMT) )
	return checkObj!(TFIVolume, TFIVolumeMT)(L, index);
    // TODO MeshVolume...
    // if no match found then
    return null;
}

extern(C) int isVolume(lua_State* L)
{
    if ( checkVolume(L, 1) )
	lua_pushboolean(L, 1);
    else
	lua_pushboolean(L, 0);
    return 1;
}

extern(C) int opCallVolume(T, string MTname)(lua_State* L)
{
    auto volume = checkObj!(T, MTname)(L, 1);
    auto r = luaL_checknumber(L, 2);
    auto s = luaL_checknumber(L, 3);
    auto t = luaL_checknumber(L, 4);
    auto pt = volume(r, s, t);
    return pushVector3(L, pt);
}

// Helper functions for constructors.

ParametricSurface[] get6Surfaces(lua_State *L, string ctorName)
{
    // Assume that table containing the Surfaces is at top of stack.
    string errMsgTmplt = `Error in call to %s:new.
The value set for the face[%d] was not of ParametricSurface type.`;
    ParametricSurface[] faces;
    foreach(i; 0 .. 6) {
	lua_rawgeti(L, -1, i+1);
	faces ~= checkSurface(L, -1);
	if ( faces[i] is null ) {
	    luaL_error(L, toStringz(format(errMsgTmplt, ctorName, i)));
	}
	lua_pop(L, 1);
    }
    return faces;
} // end get6surfaces()

Vector3[] get8Vector3s(lua_State *L, string ctorName)
{
    // Assume that table containing the Vector3 objects is at top of stack.
    string errMsgTmplt = `Error in call to %s:new.
The value set for the corner[%d] was not of Vector3 type.`;
    Vector3[] corners;
    foreach(i; 0 .. 8) {
	lua_rawgeti(L, -1, i+1);
	auto ptr = checkVector3(L, -1);
	if ( ptr is null ) {
	    luaL_error(L, toStringz(format(errMsgTmplt, ctorName, i)));
	} else {
	    corners ~= *ptr;
	}
	lua_pop(L, 1);
    }
    return corners;
} // end get8Vector3s()


// Constructor for the TFIVolume, to be used from the Lua domain.
//
// Supported forms are:
// vol0 = TFIVolume:new{faces={nFace, eFace, sFace, wFace, tFace, bFace}}
// vol1 = TFIVolume:new{vertices={p0, p1, p2, p3, p4, p5, p6, p7}}
//
// Notes:
// 1. See PJs diagram at top of geom.volume.d for ordering and labelling of
//    paths and corners.
// 2. No mix-n-match of constructors allowed.
//    It is one of: 6 faces OR 8 corner points.
//    If the faces table is found first, that constructor wins.

extern(C) int newTFIVolume(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    
    if ( !lua_istable(L, 1) ) {
	string errMsg = `Error in constructor TFIVolume:new.
A table with input parameters is expected as the first argument.`;
	luaL_error(L, errMsg.toStringz);
    }
    // Look for an array of ParametricSurfaces. 
    // If found, proceed with construction from these faces.
    lua_getfield(L, 1, "faces");
    if ( lua_istable(L, -1) ) {
	auto faces = get6Surfaces(L, "TFIVolume");
	lua_pop(L, 1);
	auto tfivolume = new TFIVolume(faces);
	volumeStore ~= pushObj!(TFIVolume, TFIVolumeMT)(L, tfivolume);
	return 1;
    } else {
	lua_pop(L, 1); // get rid of the non-table item
    }
    // Instead, look for an array of Vector3 objects.
    lua_getfield(L, 1, "vertices");
    if ( lua_istable(L, -1) ) {
	auto corners = get8Vector3s(L, "TFIVolume");
	lua_pop(L, 1);
	auto tfivolume = new TFIVolume(corners);
	volumeStore ~= pushObj!(TFIVolume, TFIVolumeMT)(L, tfivolume);
	return 1;
    }
    lua_pop(L, 1);
    // If we make it here, there's been an error in construction.
    string errMsg = `There's a problem in call to TFIVolume.new.
Neither the array of 6 surfaces nor a list of 8 vertices was found.`;
    luaL_error(L, errMsg.toStringz);
    return 0;
} // end newTFIVolume()


void getRSandT(lua_State* L, string ctorName,
	       out double r0, out double r1,
	       out double s0, out double s1,
	       out double t0, out double t1)
{
    // Table is at index 1
    string errMsgTmplt = format("Error in call to %s:new.", ctorName);
    errMsgTmplt ~= "The value for variable '%s' is not valid. It should be a number value.";
    r0 = getNumberFromTable(L, 1, "r0", false, 0.0, true, format(errMsgTmplt, "r0"));
    r1 = getNumberFromTable(L, 1, "r1", false, 1.0, true, format(errMsgTmplt, "r1"));
    s0 = getNumberFromTable(L, 1, "s0", false, 0.0, true, format(errMsgTmplt, "s0"));
    s1 = getNumberFromTable(L, 1, "s1", false, 1.0, true, format(errMsgTmplt, "s1"));
    t0 = getNumberFromTable(L, 1, "t0", false, 0.0, true, format(errMsgTmplt, "t0"));
    t1 = getNumberFromTable(L, 1, "t1", false, 1.0, true, format(errMsgTmplt, "t1"));
} // end getRSandT()

// Constructor for the SubRangedVolume, to be used from the Lua domain.
//
// Supported form is:
// vol1 = SubRangedVolume:new{vol0}

extern(C) int newSubRangedVolume(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    
    if ( !lua_istable(L, 1) ) {
	string errMsg = `Error in constructor SubRangeVolume:new.
A table with input parameters is expected as the first argument.`;
	luaL_error(L, errMsg.toStringz);
    }
    // Look for the original ParametricVolume in the first array position. 
    lua_rawgeti(L, 1, 1);
    if ( lua_isnil(L, -1) ) {
	string errMsg = `Error in call to SubRangedVolume:new{}. No table entry found.`;
	luaL_error(L, errMsg.toStringz());
    }
    auto pvolume = checkVolume(L, -1);
    lua_pop(L, 1);
    if ( pvolume is null ) {
	string errMsg = `Error in call to SubRangedVolume:new{}. No valid Surface object found.`;
	luaL_error(L, errMsg.toStringz());
    }
    double r0, r1, s0, s1, t0, t1;
    getRSandT(L, "SubRangedVolume", r0, r1, s0, s1, t0, t1);
    auto srvolume = new SubRangedVolume(pvolume, r0, r1, s0, s1, t0, t1);
    volumeStore ~= pushObj!(SubRangedVolume, SubRangedVolumeMT)(L, srvolume);
    return 1;
} // end newSubRangedVolume()


void registerVolumes(lua_State* L)
{
    // Register the TFIVolume object
    luaL_newmetatable(L, TFIVolumeMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newTFIVolume);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallVolume!(TFIVolume, TFIVolumeMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallVolume!(TFIVolume, TFIVolumeMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(TFIVolume, TFIVolumeMT));
    lua_setfield(L, -2, "__tostring");

    lua_setglobal(L, TFIVolumeMT.toStringz);

    // Register the SubRangedVolume object
    luaL_newmetatable(L, SubRangedVolumeMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newSubRangedVolume);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallVolume!(SubRangedVolume, SubRangedVolumeMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallVolume!(SubRangedVolume, SubRangedVolumeMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(SubRangedVolume, SubRangedVolumeMT));
    lua_setfield(L, -2, "__tostring");

    lua_setglobal(L, SubRangedVolumeMT.toStringz);

    // Utility functions...
    lua_pushcfunction(L, &isVolume);
    lua_setglobal(L, "isVolume");
} // end registerVolumes()
