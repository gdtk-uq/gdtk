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
immutable string LuaFnVolumeMT = "LuaFnVolume";
immutable string SubRangedVolumeMT = "SubRangedVolume";
// TODO MeshVolume...

static const(ParametricVolume)[] volumeStore;

ParametricVolume checkVolume(lua_State* L, int index) {
    // We have to do a brute force test for each object type, in turn.
    if ( isObjType(L, index, TFIVolumeMT) )
	return checkObj!(TFIVolume, TFIVolumeMT)(L, index);
    if ( isObjType(L, index, TFIVolumeMT) )
	return checkObj!(LuaFnVolume, LuaFnVolumeMT)(L, index);
    if ( isObjType(L, index, TFIVolumeMT) )
	return checkObj!(SubRangedVolume, SubRangedVolumeMT)(L, index);
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
    string errMsgTmplt = "Error in call to %s:new. " ~
	"The value set for the face[%s] was not of ParametricSurface type.";
    string[] face_names = ["north", "east", "south", "west", "top", "bottom"];
    ParametricSurface[] faces;
    foreach(i, name; face_names) {
	lua_getfield(L, -1, name.toStringz());
	faces ~= checkSurface(L, -1);
	if ( faces[i] is null ) {
	    luaL_error(L, toStringz(format(errMsgTmplt, ctorName, name)));
	}
	lua_pop(L, 1);
    }
    return faces;
} // end get6surfaces()

Vector3[] get8Vector3s(lua_State *L, string ctorName)
{
    // Assume that table containing the Vector3 objects is at top of stack.
    string errMsgTmplt = "Error in call to %s:new. " ~
	"The value set for the corner[%d] was not of Vector3 type.";
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
// vol0 = TFIVolume:new{north=nFace, east=eFace, south=sFace, west=wFace,
//                      top=tFace, bottom=bFace}
// vol1 = TFIVolume:new{vertices={p0, p1, p2, p3, p4, p5, p6, p7}}
//
// Notes:
// 1. See PJs diagram at top of geom.volume.d for ordering and labelling of
//    paths and corners.
// 2. No mix-n-match of constructors allowed.
//    It is one of: 6 named faces OR 8 corner points.
//    If the north face is found first, that constructor wins.

extern(C) int newTFIVolume(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    
    if ( !lua_istable(L, 1) ) {
	string errMsg = "Error in constructor TFIVolume:new. " ~
	    "A table with input parameters is expected as the first argument.";
	luaL_error(L, errMsg.toStringz);
    }
    // Look for named ParametricSurfaces. 
    // If found, proceed with construction from these faces.
    lua_getfield(L, 1, "north"); // test item
    if ( !lua_isnil(L, -1) ) {
	lua_pop(L, 1); // discard the test item
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
    string errMsg = "There's a problem in call to TFIVolume.new. " ~
	"Neither the set of 6 named surfaces nor a list of 8 vertices was found.";
    luaL_error(L, errMsg.toStringz);
    return 0;
} // end newTFIVolume()


/**
 * LuaFnVolume class and it's Lua constructor.
 *
 * This is hangs onto a Lua call-back function that is invoked from the D domain.
 *
 * Example:
 * function myLuaFunction(r, s, t)
 *    -- Simple cube
 *    return {x=r, y=s, z=t}
 * end
 * myPath = LuaFnVolume:new{luaFnName="myLuaFunction"}
 */

class LuaFnVolume : ParametricVolume {
public:
    lua_State* L; // a pointer to the Lua interpreter's state.
    // Even though some of the class methods claim that they don't change
    // the object state, we have to get the Lua interpreter to evaluate
    // things and that diddles with the Lua interpreter's internal state.
    // So the const on the lua_State pointer is more a statement that
    // "I'm not going to switch interpreters on you."
    // Hence the ugly but (hopefully safe) casts where ever we get 
    // the Lua interpreter to do something.
    // This is the best I can do for the moment.  PJ, 2014-04-22, 2015-02-27
    string luaFnName;
    this(const lua_State* L, string luaFnName)
    {
	this.L = cast(lua_State*)L;
	this.luaFnName = luaFnName;
    }
    this(ref const(LuaFnVolume) other)
    {
	L = cast(lua_State*)other.L;
	luaFnName = other.luaFnName;
    }
    override LuaFnVolume dup() const
    {
	return new LuaFnVolume(L, luaFnName);
    }
    override Vector3 opCall(double r, double s, double t) const 
    {
	// Call back to the Lua function.
	lua_getglobal(cast(lua_State*)L, luaFnName.toStringz);
	lua_pushnumber(cast(lua_State*)L, r);
	lua_pushnumber(cast(lua_State*)L, s);
	lua_pushnumber(cast(lua_State*)L, t);
	if ( lua_pcall(cast(lua_State*)L, 3, 1, 0) != 0 ) {
	    string errMsg = "Error in call to " ~ luaFnName ~ 
		" from LuaFnVolume:opCall(): " ~ 
		to!string(lua_tostring(cast(lua_State*)L, -1));
	    luaL_error(cast(lua_State*)L, errMsg.toStringz);
	}
	// We are expecting a table to be returned, containing three numbers.
	if ( !lua_istable(cast(lua_State*)L, -1) ) {
	    string errMsg = "Error in call to LuaFnVolume:opCall().; " ~
		"A table containing arguments is expected, but no table was found.";
	    luaL_error(cast(lua_State*)L, errMsg.toStringz);
	}
	double x = 0.0; // default value
	lua_getfield(cast(lua_State*)L, -1, "x".toStringz());
	if ( lua_isnumber(cast(lua_State*)L, -1) ) {
	    x = to!double(lua_tonumber(cast(lua_State*)L, -1));
	}
	lua_pop(cast(lua_State*)L, 1);
	double y = 0.0; // default value
	lua_getfield(cast(lua_State*)L, -1, "y".toStringz());
	if ( lua_isnumber(cast(lua_State*)L, -1) ) {
	    y = to!double(lua_tonumber(cast(lua_State*)L, -1));
	}
	lua_pop(cast(lua_State*)L, 1);
	double z = 0.0; // default value
	lua_getfield(cast(lua_State*)L, -1, "z".toStringz());
	if ( lua_isnumber(cast(lua_State*)L, -1) ) {
	    z = to!double(lua_tonumber(cast(lua_State*)L, -1));
	}
	lua_pop(cast(lua_State*)L, 1);
	//
	lua_settop(cast(lua_State*)L, 0); // clear the stack
	return Vector3(x, y, z);
    } // end opCall()
    override string toString() const
    {
	return "LuaFnVolume(luaFnName=\"" ~ luaFnName ~ "\")";
    }
} // end class LuaFnVolume

extern(C) int newLuaFnVolume(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = "Error in call to LuaFnVolume:new{}.; " ~
	    "A table containing arguments is expected, but no table was found.";
	luaL_error(L, errMsg.toStringz);
    }
    // Expect function name in table.
    string fnName = "";
    lua_getfield(L, 1, "luaFnName".toStringz());
    if ( lua_isnil(L, -1) ) {
	string errMsg = "Error in call to LuaFnVolume:new{}. No luaFnName entry found.";
	luaL_error(L, errMsg.toStringz());
    }
    if ( lua_isstring(L, -1) ) {
	fnName ~= to!string(lua_tostring(L, -1));
    }
    lua_pop(L, 1);
    if ( fnName == "" ) {
	string errMsg = "Error in call to LuaFnVolume:new{}. No function name found.";
	luaL_error(L, errMsg.toStringz());
    }
    auto lfv = new LuaFnVolume(L, fnName);
    volumeStore ~= pushObj!(LuaFnVolume, LuaFnVolumeMT)(L, lfv);
    return 1;
} // end newLuaFnVolume()


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
// vol1 = SubRangedVolume:new{underlying_pvolume=vol0, r0=r0value, r1=r1value,
//                            s0=s0value, s1=s1vlaue, t0=t0value, t1=t1value}

extern(C) int newSubRangedVolume(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    
    if ( !lua_istable(L, 1) ) {
	string errMsg = "Error in constructor SubRangeVolume:new. " ~
	    "A table with input parameters is expected as the first argument.";
	luaL_error(L, errMsg.toStringz);
    }
    // Look for the underlying ParametricVolume.
    lua_getfield(L, 1, "underlying_pvolume");
    if ( lua_isnil(L, -1) ) {
	string errMsg = "Error in call to SubRangedVolume:new{}. No underlying_pvolume field found.";
	luaL_error(L, errMsg.toStringz());
    }
    auto pvolume = checkVolume(L, -1);
    lua_pop(L, 1);
    if ( pvolume is null ) {
	string errMsg = "Error in call to SubRangedVolume:new{}. No valid Volume object found.";
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

    // Register the LuaFnVolume object
    luaL_newmetatable(L, LuaFnVolumeMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newLuaFnVolume);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallVolume!(LuaFnVolume, LuaFnVolumeMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallVolume!(LuaFnVolume, LuaFnVolumeMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(LuaFnVolume, LuaFnVolumeMT));
    lua_setfield(L, -2, "__tostring");

    lua_setglobal(L, LuaFnVolumeMT.toStringz);

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
