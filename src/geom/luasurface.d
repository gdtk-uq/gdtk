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
import luasgrid;

// Name of metatables -- these are the Lua access names.
immutable string CoonsPatchMT = "CoonsPatch";
immutable string AOPatchMT = "AOPatch";
immutable string ChannelPatchMT = "ChannelPatch";
immutable string SweptPathPatchMT = "SweptPathPatch";
immutable string MeshPatchMT = "MeshPatch";
immutable string LuaFnSurfaceMT = "LuaFnSurface";
immutable string SubRangedSurfaceMT = "SubRangedSurface";

static const(ParametricSurface)[] surfaceStore;

ParametricSurface checkSurface(lua_State* L, int index) {
    // We have to do a brute force test for each object type, in turn.
    if ( isObjType(L, index, CoonsPatchMT) )
	return checkObj!(CoonsPatch, CoonsPatchMT)(L, index);
    if ( isObjType(L, index, AOPatchMT ) )
	return checkObj!(AOPatch, AOPatchMT)(L, index);
    if ( isObjType(L, index, ChannelPatchMT ) )
	return checkObj!(ChannelPatch, ChannelPatchMT)(L, index);
    if ( isObjType(L, index, SweptPathPatchMT ) )
	return checkObj!(SweptPathPatch, SweptPathPatchMT)(L, index);
    if ( isObjType(L, index, MeshPatchMT ) )
	return checkObj!(MeshPatch, MeshPatchMT)(L, index);
    if ( isObjType(L, index, LuaFnSurfaceMT ) )
	return checkObj!(LuaFnSurface, LuaFnSurfaceMT)(L, index);
    if ( isObjType(L, index, SubRangedSurfaceMT ) )
	return checkObj!(SubRangedSurface, SubRangedSurfaceMT)(L, index);
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

void getPaths(lua_State *L, string ctorName, out Path[string] paths)
{
    string errMsgTmplt = "Error in call to %s:new. " ~
	"The value set for the '%s' path was not of Path type."; 
    int index = 1;  // Assume that table is at index 1.
    string[] path_names = ["north", "east", "south", "west"];
    foreach (name; path_names) {
	lua_getfield(L, index, name.toStringz());
	paths[name] = checkPath(L, -1);
	if (paths[name] is null) luaL_error(L, toStringz(format(errMsgTmplt, ctorName, name)));
	lua_pop(L, 1);
    }
}

void getVector3s(lua_State *L, string ctorName, out Vector3[string] corners)
{
    string errMsgTmplt = "Error in call to %s:new. " ~
	"The value set for the '%s' corner was not of Vector3 type."; 
    int index = 1;  // Assume that table is at index 1.
    string[] corner_names = ["p00", "p10", "p11", "p01"];
    foreach (name; corner_names) {
	lua_getfield(L, index, name.toStringz());
	auto p = checkVector3(L, -1);
	if (p is null) {
	    luaL_error(L, toStringz(format(errMsgTmplt, ctorName, name)));
	} else {
	    corners[name] = *p;
	}
	lua_pop(L, 1);
    }
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
 * 2. No mix-n-match of constructors allowed. It is one of: 4 paths OR 4 corner points. 
 *    If a path is found first, that constructor wins.
 *
 */

extern(C) int newCoonsPatch(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    
    if ( !lua_istable(L, 1) ) {
	string errMsg = "Error in constructor CoonPatch:new. " ~
	    "A table with input parameters is expected as the first argument.";
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
	Vector3[string] corners;
	getVector3s(L, "CoonsPatch", corners);
	auto cpatch = new CoonsPatch(corners["p00"], corners["p10"],
				     corners["p11"], corners["p01"]);
	surfaceStore ~= pushObj!(CoonsPatch, CoonsPatchMT)(L, cpatch);
	return 1;
    }
    lua_pop(L, 1);
    // If we make it here, there's been an error in construction.
    string errMsg = "There's a problem in call to CoonsPatch.new. " ~
	"Neither a list of named paths ('north', 'east', 'south', 'west') " ~
	"nor a list of named corners ('p00', 'p10', 'p11', 'p01') were found.";
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
 * 2. No mix-n-match of constructors allowed. It is one of: 4 paths OR 4 corner points.
 *    If a path is found first, that constructor wins.
 *
 */
extern(C) int newAOPatch(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    
    if ( !lua_istable(L, 1) ) {
	string errMsg = "Error in constructor AOPatch:new. " ~
	    "A table with input parameters is expected as the first argument.";
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
	auto aopatch = new AOPatch(paths["south"], paths["north"], paths["west"], paths["east"], nx, ny);
	surfaceStore ~= pushObj!(AOPatch, AOPatchMT)(L, aopatch);
	return 1;
    } else {
	lua_pop(L, 1);
    }
    // Instead, look for Vector3 objects.
    lua_getfield(L, 1, "p00");
    if ( ! lua_isnil(L, -1) ) {
	lua_pop(L, 1);
	Vector3[string] corners;
	getVector3s(L, "AOPatch", corners);
	auto aopatch = new AOPatch(corners["p00"], corners["p10"], corners["p11"], corners["p01"], nx, ny);
	surfaceStore ~= pushObj!(AOPatch, AOPatchMT)(L, aopatch);
	return 1;
    }
    lua_pop(L, 1);
    // If we make it here, there's been an error in construction.
    string errMsg = "There's a problem in call to AOPatch.new. " ~
	"Neither a list of named paths ('north', 'east', 'south', 'west') " ~
	"nor a list of named corners ('p00', 'p10', 'p11', 'p01') were found.";
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
	string errMsg = "Error in constructor ChannelPatch:new. " ~
	    "A table with input parameters is expected as the first argument.";
	luaL_error(L, errMsg.toStringz);
    }
    // Look for south and north paths.
    lua_getfield(L, 1, "south");
    auto south = checkPath(L, -1);
    if ( south is null ) {
	string errMsg = "Error in constructor ChannelPatch:new. Couldn't find south Path.";
	luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);
    lua_getfield(L, 1, "north");
    auto north = checkPath(L, -1);
    if ( north is null ) {
	string errMsg = "Error in constructor ChannelPatch:new. Couldn't find north Path.";
	luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);
    // Optional parameters
    string errMsgTmpltBool = "Error in call to ChannelPatch:new{}. " ~
	"A valid value for '%s' was not found in list of arguments. " ~
	"The value, if present, should be boolean (true or false).";
    bool ruled = getBooleanFromTable(L, 1, "ruled", false, false, true, format(errMsgTmpltBool, "ruled"));
    bool pure2D = getBooleanFromTable(L, 1, "pure2D", false, false, true, format(errMsgTmpltBool, "pure2D"));
    // Construct the actual surface.
    auto cpatch = new ChannelPatch(south, north, ruled, pure2D);
    surfaceStore ~= pushObj!(ChannelPatch, ChannelPatchMT)(L, cpatch);
    return 1;
} // end newChannelPatch()


/**
 * This is the constructor for a SweptPathPatch to be used from the Lua interface.
 *
 * At successful completion of this function, a new SweptPathPatch object
 * is pushed onto the Lua stack.
 *
 * Supported constructions are:
 * -------------------------
 * patch0 = SweptPathPatch:new{west=wPath, south=sPath}
 * --------------------------
 */

extern(C) int newSweptPathPatch(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    
    if (!lua_istable(L, 1)) {
	string errMsg = "Error in constructor SweptPathPatch:new. " ~
	    "A table with input parameters is expected as the first argument.";
	luaL_error(L, errMsg.toStringz);
    }
    // Look for west and south paths.
    lua_getfield(L, 1, "west");
    auto west = checkPath(L, -1);
    if (west is null) {
	string errMsg = "Error in constructor SweptPathPatch:new. Couldn't find west Path.";
	luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);
    lua_getfield(L, 1, "south");
    auto south = checkPath(L, -1);
    if (south is null) {
	string errMsg = "Error in constructor SweptPathPatch:new. Couldn't find south Path.";
	luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);
    // Construct the actual surface.
    auto spatch = new SweptPathPatch(west, south);
    surfaceStore ~= pushObj!(SweptPathPatch, SweptPathPatchMT)(L, spatch);
    return 1;
} // end newSweptPathPatch()


/**
 * This is the constructor for a MeshPatch to be used from the Lua interface.
 *
 * At successful completion of this function, a new MeshPatch object
 * is pushed onto the Lua stack.
 *
 * Supported constructions are:
 * -------------------------
 * patch0 = MeshPatch:new{sgrid=myStructuredGrid}
 * --------------------------
 */

extern(C) int newMeshPatch(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    
    if ( !lua_istable(L, 1) ) {
	string errMsg = "Error in constructor MeshPatch:new. " ~
	    "A table with input parameters is expected as the first argument.";
	luaL_error(L, errMsg.toStringz);
    }
    // Look for the StructuredGrid object.
    lua_getfield(L, 1, "sgrid");
    auto grid = checkStructuredGrid(L, -1);
    if ( grid is null ) {
	string errMsg = "Error in constructor MeshPatch:new. " ~
	    "Couldn't find StructuredGrid object in sgrid field.";
	luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);
    // Construct the actual surface.
    auto mpatch = new MeshPatch(grid);
    surfaceStore ~= pushObj!(MeshPatch, MeshPatchMT)(L, mpatch);
    return 1;
} // end newMeshPatch()


/**
 * LuaFnSurface class and it's Lua constructor.
 *
 * This is hangs onto a Lua call-back function that is invoked from the D domain.
 *
 * Example:
 * function myLuaFunction(r, s)
 *    -- Simple plane
 *    return {x=r, y=s, z=0.0}
 * end
 * mySurf = LuaFnSurface:new{luaFnName="myLuaFunction"}
 */

class LuaFnSurface : ParametricSurface {
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
    this(ref const(LuaFnSurface) other)
    {
	L = cast(lua_State*)other.L;
	luaFnName = other.luaFnName;
    }
    override LuaFnSurface dup() const
    {
	return new LuaFnSurface(L, luaFnName);
    }
    override Vector3 opCall(double r, double s) const 
    {
	// Call back to the Lua function.
	lua_getglobal(cast(lua_State*)L, luaFnName.toStringz);
	lua_pushnumber(cast(lua_State*)L, r);
	lua_pushnumber(cast(lua_State*)L, s);
	if ( lua_pcall(cast(lua_State*)L, 2, 1, 0) != 0 ) {
	    string errMsg = "Error in call to " ~ luaFnName ~ 
		" from LuaFnSurface:opCall(): " ~ 
		to!string(lua_tostring(cast(lua_State*)L, -1));
	    luaL_error(cast(lua_State*)L, errMsg.toStringz);
	}
	// We are expecting a table to be returned, containing three numbers.
	if ( !lua_istable(cast(lua_State*)L, -1) ) {
	    string errMsg = "Error in call to LuaFnSurface:opCall().; " ~
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
	return "LuaFnSurface(luaFnName=\"" ~ luaFnName ~ "\")";
    }
} // end class LuaFnSurface

extern(C) int newLuaFnSurface(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = "Error in call to LuaFnSurface:new{}.; " ~
	    "A table containing arguments is expected, but no table was found.";
	luaL_error(L, errMsg.toStringz);
    }
    // Expect function name in table.
    string fnName = "";
    lua_getfield(L, 1, "luaFnName".toStringz());
    if ( lua_isnil(L, -1) ) {
	string errMsg = "Error in call to LuaFnSurface:new{}. No luaFnName entry found.";
	luaL_error(L, errMsg.toStringz());
    }
    if ( lua_isstring(L, -1) ) {
	fnName ~= to!string(lua_tostring(L, -1));
    }
    lua_pop(L, 1);
    if ( fnName == "" ) {
	string errMsg = "Error in call to LuaFnSurface:new{}. No function name found.";
	luaL_error(L, errMsg.toStringz());
    }
    auto lfs = new LuaFnSurface(L, fnName);
    surfaceStore ~= pushObj!(LuaFnSurface, LuaFnSurfaceMT)(L, lfs);
    return 1;
} // end newLuaFnSurface()


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
    lua_getfield(L, 1, "underlying_psurface");
    if ( lua_isnil(L, -1) ) {
	string errMsg = "Error in call to SubRangedSurface:new{}. No underlying_psurface field found.";
	luaL_error(L, errMsg.toStringz());
    }
    auto psurf = checkSurface(L, -1);
    lua_pop(L, 1);
    if ( psurf is null ) {
	string errMsg = "Error in call to SubRangedSurface:new{}. No valid Surface object found.";
	luaL_error(L, errMsg.toStringz());
    }
    double r0, r1, s0, s1;
    getRandS(L, "SubRangedSurface", r0, r1, s0, s1);
    auto srs = new SubRangedSurface(psurf, r0, r1, s0, s1);
    surfaceStore ~= pushObj!(SubRangedSurface, SubRangedSurfaceMT)(L, srs);
    return 1;
} // end newSubRangedSurface()

/* ---------- convenience functions -------------- */

extern(C) int makePatch(lua_State* L)
{
    if ( !lua_istable(L, 1) ) {
	string errMsg = "Error in call to makePatch. " ~
	    "A table is expected as the first argument. No table was found.";
	luaL_error(L, errMsg.toStringz);
    }
    // Do some error checking of arguments.
    if ( lua_objlen(L, 1) != 0 ) {
	// It appears that the caller has tried to set arguments as an array
	string errMsg = "Error in call to makePatch. " ~
	    "Name parameters are expected. An array style of parameters was found.";
	luaL_error(L, errMsg.toStringz);
    }
    // Get boundary paths.
    string[] edges = ["north", "east", "south", "west"];
    Path[string] paths;
    string errMsgTmpl = "Error in call to makePatch. The %s path is not a valid Path object.";
    foreach (edge_name; edges) {
	lua_getfield(L, 1, edge_name.toStringz());
	auto my_path = checkPath(L, -1);
	if (my_path is null) {
	    luaL_error(L, toStringz(format(errMsgTmpl, edge_name)));
	} else {
	    paths[edge_name] = my_path;
	}
	lua_pop(L, 1);
    }
    string gridType = "TFI";
    lua_getfield(L, 1, "gridType");
    if ( !lua_isnil(L, -1) ) {
	gridType = to!string(luaL_checkstring(L, -1));
    }
    lua_pop(L, 1);
    gridType = toUpper(gridType);
    if ( gridType == "AO" ) {
	auto patch = new AOPatch(paths["south"], paths["north"], paths["west"], paths["east"]);
	surfaceStore ~= pushObj!(AOPatch, AOPatchMT)(L, patch);
    } else { // Default to CoonsPatch
	auto patch = new CoonsPatch(paths["south"], paths["north"], paths["west"], paths["east"]);
	surfaceStore ~= pushObj!(CoonsPatch, CoonsPatchMT)(L, patch);
    }
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

    // Register the SweptPathPatch object
    luaL_newmetatable(L, SweptPathPatchMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newSweptPathPatch);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallSurface!(SweptPathPatch, SweptPathPatchMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallSurface!(SweptPathPatch, SweptPathPatchMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(SweptPathPatch, SweptPathPatchMT));
    lua_setfield(L, -2, "__tostring");

    lua_setglobal(L, SweptPathPatchMT.toStringz);

    // Register the MeshPatch object
    luaL_newmetatable(L, MeshPatchMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newMeshPatch);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallSurface!(MeshPatch, MeshPatchMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallSurface!(MeshPatch, MeshPatchMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(MeshPatch, MeshPatchMT));
    lua_setfield(L, -2, "__tostring");

    lua_setglobal(L, MeshPatchMT.toStringz);

    // Register the LuaFnSurface object
    luaL_newmetatable(L, LuaFnSurfaceMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newLuaFnSurface);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallSurface!(LuaFnSurface, LuaFnSurfaceMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallSurface!(LuaFnSurface, LuaFnSurfaceMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(LuaFnSurface, LuaFnSurfaceMT));
    lua_setfield(L, -2, "__tostring");

    lua_setglobal(L, LuaFnSurfaceMT.toStringz);

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
