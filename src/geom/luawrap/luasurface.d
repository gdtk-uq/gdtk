/**
 * luasurface.d
 * Lua interface to ParametricSurface objects.
 *
 * Authors: Rowan G. and Peter J.
 * Version: 2015-02-24 first code
 *          2015-07-05 SubRangedSurface
 */

module geom.luawrap.luasurface;

import std.stdio;
import std.string;
import std.conv;
import std.uni;
import util.lua;
import util.lua_service;
import ntypes.complex;
import nm.number;
import geom;
import geom.luawrap.luageom;
import geom.luawrap.luagpath;
import geom.luawrap.luasgrid;

// Name of metatables -- these are the Lua access names.
immutable string CoonsPatchMT = "CoonsPatch";
immutable string AOPatchMT = "AOPatch";
immutable string ChannelPatchMT = "ChannelPatch";
immutable string RuledSurfaceMT = "RuledSurface";
immutable string NozzleExpansionPatchMT = "NozzleExpansionPatch";
immutable string SweptPathPatchMT = "SweptPathPatch";
immutable string SpherePatchMT = "SpherePatch";
immutable string CubePatchMT = "CubePatch";
immutable string MeshPatchMT = "MeshPatch";
immutable string LuaFnSurfaceMT = "LuaFnSurface";
immutable string SubRangedSurfaceMT = "SubRangedSurface";
immutable string BezierPatchMT = "BezierPatch";
immutable string ControlPointPatchMT = "ControlPointPatch";
immutable string BezierTrianglePatchMT = "BezierTrianglePatch";
immutable string NURBSPatchMT = "NURBSPatch";

static const(ParametricSurface)[] surfaceStore;
static const(BezierTrianglePatch)[] triPatchStore;

ParametricSurface checkSurface(lua_State* L, int index) {
    // We have to do a brute force test for each object type, in turn.
    if ( isObjType(L, index, CoonsPatchMT) )
        return checkObj!(CoonsPatch, CoonsPatchMT)(L, index);
    if ( isObjType(L, index, AOPatchMT ) )
        return checkObj!(AOPatch, AOPatchMT)(L, index);
    if ( isObjType(L, index, ChannelPatchMT ) )
        return checkObj!(ChannelPatch, ChannelPatchMT)(L, index);
    if ( isObjType(L, index, RuledSurfaceMT ) )
        return checkObj!(RuledSurface, RuledSurfaceMT)(L, index);
    if ( isObjType(L, index, NozzleExpansionPatchMT ) )
        return checkObj!(NozzleExpansionPatch, NozzleExpansionPatchMT)(L, index);
    if ( isObjType(L, index, SweptPathPatchMT ) )
        return checkObj!(SweptPathPatch, SweptPathPatchMT)(L, index);
    if ( isObjType(L, index, SpherePatchMT) )
        return checkObj!(SpherePatch, SpherePatchMT)(L, index);
    if ( isObjType(L, index, CubePatchMT) )
        return checkObj!(CubePatch, CubePatchMT)(L, index);
    if ( isObjType(L, index, MeshPatchMT ) )
        return checkObj!(MeshPatch, MeshPatchMT)(L, index);
    if ( isObjType(L, index, LuaFnSurfaceMT ) )
        return checkObj!(LuaFnSurface, LuaFnSurfaceMT)(L, index);
    if ( isObjType(L, index, SubRangedSurfaceMT ) )
        return checkObj!(SubRangedSurface, SubRangedSurfaceMT)(L, index);
    if ( isObjType(L, index, BezierPatchMT) )
        return checkObj!(BezierPatch, BezierPatchMT)(L, index);
    if ( isObjType(L, index, ControlPointPatchMT) )
        return checkObj!(ControlPointPatch, ControlPointPatchMT)(L, index);
    if ( isObjType(L, index, NURBSPatchMT) )
        return checkObj!(NURBSSurface, NURBSPatchMT)(L, index);

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

extern(C) int areaOfSurface(T, string MTname)(lua_State* L)
{
    int nargs = lua_gettop(L);
    auto surface = checkObj!(T, MTname)(L, 1);
    int nr = 10;
    if (nargs > 1) { nr = luaL_checkint(L, 2); }
    int ns = 10;
    if (nargs > 2) { ns = luaL_checkint(L, 3); }
    return pushVector3(L, surface.area(nr, ns));
} // end areaOfSurface

extern(C) int luaWriteSurfaceAsVtkXml(lua_State *L)
{
    int narg = lua_gettop(L);
    if (narg < 4) {
        string errMsg = "Not enough arguments passed to writeSurfaceAsVtkXml\n";
        luaL_error(L, errMsg.toStringz);
    }
    auto surface = checkSurface(L, 1);
    string fName = to!string(luaL_checkstring(L, 2));
    int nrPts = luaL_checkint(L, 3);
    int nsPts = luaL_checkint(L, 4);
    writeSurfaceAsVtkXml(surface, fName, nrPts, nsPts);
    return 0;
}


// Helper functions for constructors.

void getPaths(lua_State *L, string ctorName, out Path[string] paths, bool allowNull=false)
{
    string errMsgTmplt = "Error in call to %s:new. " ~
        "The value set for the '%s' path was not of Path type.";
    int index = 1;  // Assume that table is at index 1.
    string[] path_names = ["north", "east", "south", "west"];
    foreach (name; path_names) {
        lua_getfield(L, index, name.toStringz());
        paths[name] = checkPath(L, -1);
        if (paths[name] is null && !allowNull) {
            luaL_error(L, toStringz(format(errMsgTmplt, ctorName, name)));
        }
        lua_pop(L, 1);
    }
}

void getVector3s(lua_State *L, string ctorName, out Vector3[string] corners)
{
    string errMsgTmplt = "Error in call to %s:new. " ~
        "A value for the '%s' corner was not available.";
    int index = 1;  // Assume that table is at index 1.
    string[] corner_names = ["p00", "p10", "p11", "p01"];
    foreach (name; corner_names) {
        lua_getfield(L, index, name.toStringz());
        if (lua_isnil(L, -1)) { luaL_error(L, toStringz(format(errMsgTmplt, ctorName, name))); }
        auto p = toVector3(L, -1);
        corners[name] = p;
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
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected CoonsPatch:new{}; ";
        errMsg ~= "maybe you tried CoonsPatch.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in constructor CoonsPatch:new{}. " ~
            "A table with input parameters is expected as the first argument.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["p00","p10","p11","p01","north","south","west","east"])) {
        string errMsg = "Error in call to CoonsPatch:new{}. Invalid name in table.";
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
    string errMsg = "There's a problem in call to CoonsPatch.new{}. " ~
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
 * patch0 = AOPatch:new{north=nPath, east=ePath, south=sPath, west=wPath, nx=10, ny=10}
 * patch2 = AOPatch:new{p00=a, p10=b, p11=c, p01=d, nx=10, ny=10}
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
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected AOPatch:new{}; ";
        errMsg ~= "maybe you tried AOPatch.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in constructor AOPatch:new{}. " ~
            "A table with input parameters is expected as the first argument.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["p00","p10","p11","p01","north","south","west","east","nx","ny"])) {
        string errMsg = "Error in call to AOPatch:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    string errMsgTmplt = "Error in call to AOPatch:new{}.\n";
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
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected ChannelPatch:new{}; ";
        errMsg ~= "maybe you tried ChannelPatch.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in constructor ChannelPatch:new{}. " ~
            "A table with input parameters is expected as the first argument.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["north","south","ruled","pure2D"])) {
        string errMsg = "Error in call to ChannelPatch:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    // Look for south and north paths.
    lua_getfield(L, 1, "south");
    auto south = checkPath(L, -1);
    if ( south is null ) {
        string errMsg = "Error in constructor ChannelPatch:new{}. Couldn't find south Path.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);
    lua_getfield(L, 1, "north");
    auto north = checkPath(L, -1);
    if ( north is null ) {
        string errMsg = "Error in constructor ChannelPatch:new{}. Couldn't find north Path.";
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

extern(C) int make_bridging_path_ChannelPatch(lua_State* L)
{
    auto surface = checkObj!(ChannelPatch, ChannelPatchMT)(L, 1);
    if (!surface) {
        string errMsg = "Error in ChannelPatch:make_bridging_path. Didn't get a ChannelPatch object.";
        luaL_error(L, errMsg.toStringz);
    }
    auto r = luaL_checknumber(L, 2);
    // Construct the actual path.
    auto bpath = surface.make_bridging_path(r);
    if (auto dpath = cast(Line) bpath) {
        pathStore ~= pushObj!(Line, LineMT)(L, dpath);
    } else if (auto dpath = cast(Bezier) bpath) {
        pathStore ~= pushObj!(Bezier, BezierMT)(L, dpath);
    } else {
        string errMsg = "Error in ChannelPatch:make_bridging_path. Didn't push a Path object.";
        luaL_error(L, errMsg.toStringz);
    }
    return 1;
} // end make_bridging_path_ChannelPatch()


/**
 * This is the constructor for a RuledSurface to be used from the Lua interface.
 *
 * At successful completion of this function, a new RuledSurface object
 * is pushed onto the Lua stack.
 *
 * Supported constructions are:
 * -------------------------
 * patch0 = RuledSurface:new{edge0=sPath, edge1=nPath}
 * patch1 = RuledSurface:new{edge0=sPath, edge1=nPath, ruled_direction='s', pure2D=false}
 * patch2 = RuledSurface:new{edge0=wPath, edge1=ePath, ruled_direction='r'}
 * For the ruled_direction string, 's' acts the same as 'j' and 'r' acts the same as 'i'.
 * --------------------------
 */

extern(C) int newRuledSurface(lua_State* L)
{
    int narg = lua_gettop(L);
    if (!(narg == 2 && lua_istable(L, 1))) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected RuledSurface:new{}; ";
        errMsg ~= "maybe you tried RuledSurface.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if (!lua_istable(L, 1)) {
        string errMsg = "Error in constructor RuledSurface:new{}. " ~
            "A table with input parameters is expected as the first argument.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["edge0","edge1","ruled_direction","pure2D"])) {
        string errMsg = "Error in call to RuledSurface:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    // Look for edge0 and edge1 paths.
    lua_getfield(L, 1, "edge0");
    auto edge0 = checkPath(L, -1);
    if (edge0 is null) {
        string errMsg = "Error in constructor RuledSurface:new{}. Couldn't find edge0 Path.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);
    lua_getfield(L, 1, "edge1");
    auto edge1 = checkPath(L, -1);
    if (edge1 is null) {
        string errMsg = "Error in constructor RuledSurface:new{}. Couldn't find edge1 Path.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);
    // Optional parameters
    string ruled_direction = getStringWithDefault(L, 1, "ruled_direction", "s");
    string errMsgTmpltBool = "Error in call to RuledSurface:new{}. " ~
        "A valid value for '%s' was not found in list of arguments. " ~
        "The value, if present, should be boolean (true or false).";
    bool pure2D = getBooleanFromTable(L, 1, "pure2D", false, false, true, format(errMsgTmpltBool, "pure2D"));
    // Construct the actual surface.
    auto cpatch = new RuledSurface(edge0, edge1, ruled_direction, pure2D);
    surfaceStore ~= pushObj!(RuledSurface, RuledSurfaceMT)(L, cpatch);
    return 1;
} // end newRuledSurface()


/**
 * This is the constructor for a NozzleExpansionPatch to be used from the Lua interface.
 *
 * At successful completion of this function, a new NozzleExpansionPatch object
 * is pushed onto the Lua stack.
 *
 * Supported constructions are:
 * -------------------------
 * patch = NozzleExpansionPatch:new{north=nPath}
 * --------------------------
 */

extern(C) int newNozzleExpansionPatch(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected NozzleExpansionPatch:new{}; ";
        errMsg ~= "maybe you tried NozzleExpansionPatch.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in constructor NozzleExpansionPatch:new{}. " ~
            "A table with input parameters is expected as the first argument.";
        luaL_error(L, errMsg.toStringz);
    }
    // Because of the history of this particular surface patch,
    // we will allow the user to specify all of the edges but we will use only north.
    if (!checkAllowedNames(L, 1, ["north","south","west","east"])) {
        string errMsg = "Error in call to NozzleExpansionPatch:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    // Look for north path.
    lua_getfield(L, 1, "south");
    lua_getfield(L, 1, "north");
    auto north = checkPath(L, -1);
    if (north is null) {
        string errMsg = "Error in constructor NozzleExpansionPatch:new{}. Couldn't find north Path.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);
    auto south = checkPath(L, -1);
    if (south) {
        writeln("Warning in constructor NozzleExpansionPatch:new{}. Ignored south Path.");
    }
    lua_pop(L, 1);
    lua_getfield(L, 1, "west");
    auto west = checkPath(L, -1);
    if (west) {
        writeln("Warning in constructor NozzleExpansionPatch:new{}. Ignored west Path.");
    }
    lua_pop(L, 1);
    lua_getfield(L, 1, "east");
    auto east = checkPath(L, -1);
    if (east) {
        writeln("Warning in constructor NozzleExpansionPatch:new{}. Ignored east Path.");
    }
    lua_pop(L, 1);

    // Construct the actual surface.
    auto cpatch = new NozzleExpansionPatch(north);
    surfaceStore ~= pushObj!(NozzleExpansionPatch, NozzleExpansionPatchMT)(L, cpatch);
    return 1;
} // end newNozzleExpansionPatch()

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
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected SweptPathPatch:new{}; ";
        errMsg ~= "maybe you tried SweptPathPatch.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if (!lua_istable(L, 1)) {
        string errMsg = "Error in constructor SweptPathPatch:new{}. " ~
            "A table with input parameters is expected as the first argument.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["south","west"])) {
        string errMsg = "Error in call to SweptPathPatch:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    // Look for west and south paths.
    lua_getfield(L, 1, "west");
    auto west = checkPath(L, -1);
    if (west is null) {
        string errMsg = "Error in constructor SweptPathPatch:new{}. Couldn't find west Path.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);
    lua_getfield(L, 1, "south");
    auto south = checkPath(L, -1);
    if (south is null) {
        string errMsg = "Error in constructor SweptPathPatch:new{}. Couldn't find south Path.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_pop(L, 1);
    // Construct the actual surface.
    auto spatch = new SweptPathPatch(west, south);
    surfaceStore ~= pushObj!(SweptPathPatch, SweptPathPatchMT)(L, spatch);
    return 1;
} // end newSweptPathPatch()


/**
 * This is the constructor for a SpherePatch to be used from the Lua interface.
 *
 * The model is a cube -1,1 remapped to a sphere of unit radius and then scaled
 * to radius R and moved to center C.
 * The surface starts as a single, named face of the unit cube. Names are the
 * usual north, east, south, west, top, bottom.
 * By default, you will get a full face.  You can request a half face
 * by specifying a name for the adjoining cube face and you can requast
 * a quadrant by specifying the pair of names for the adjoining cube faces.
 *
 * Supported constructions are:
 * -------------------------
 * patch0 = SpherePatch:new{radius=R, centre=C, face_name="east"}
 * patch1 = SpherePatch:new{radius=R, centre=C, face_name="east", which_part="top"}
 * patch1 = SpherePatch:new{radius=R, centre=C, face_name="east", which_part="top-south"}
 * --------------------------
 *
 * At successful completion of this function, a new SpherePatch object
 * is pushed onto the Lua stack.
 */

extern(C) int newSpherePatch(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected SpherePatch:new{}; ";
        errMsg ~= "maybe you tried SpherePatch.new{}.";
        luaL_error(L, errMsg.toStringz);
        return 0;
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in constructor SpherePatch:new{}. " ~
            "A table with input parameters is expected as the first argument.";
        luaL_error(L, errMsg.toStringz);
        return 0;
    }
    if (!checkAllowedNames(L, 1, ["radius","centre","face_name","which_part"])) {
        string errMsg = "Error in call to SpherePatch:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
        return 0;
    }
    // Now, get our table entries.
    // Table is at index 1
    double radius = getDouble(L, 1, "radius");
    string face_name = getString(L, 1, "face_name");
    string which_part = getStringWithDefault(L, 1, "which_part", "");
    lua_getfield(L, 1, "centre");
    if (lua_isnil(L, -1)) {
        string errMsg = "Error in call to SpherePatch:new{}. No centre entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto centre = toVector3(L, -1);
    lua_pop(L, 1);
    //
    auto spatch = new SpherePatch(to!number(radius), centre, face_name, which_part);
    surfaceStore ~= pushObj!(SpherePatch, SpherePatchMT)(L, spatch);
    return 1;
} // end newSpherePatch()


/**
 * This is the constructor for a CubePatch to be used from the Lua interface.
 *
 * The model is a cube -1,1 remapped to a cube of edge length 2a and center C.
 * The surface starts as a single, named face of the unit cube. Names are the
 * usual north, east, south, west, top, bottom.
 * By default, you will get a full face.  You can request a half face
 * by specifying a name for the adjoining cube face and you can requast
 * a quadrant by specifying the pair of names for the adjoining cube faces.
 *
 * Supported constructions are:
 * -------------------------
 * patch0 = CubePatch:new{a=a, centre=C, face_name="east"}
 * patch1 = CubePatch:new{a=a, centre=C, face_name="east", which_part="top"}
 * patch1 = CubePatch:new{a=a, centre=C, face_name="east", which_part="top-south"}
 * --------------------------
 *
 * At successful completion of this function, a new CubePatch object
 * is pushed onto the Lua stack.
 */

extern(C) int newCubePatch(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected CubePatch:new{}; ";
        errMsg ~= "maybe you tried CubePatch.new{}.";
        luaL_error(L, errMsg.toStringz);
        return 0;
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in constructor CubePatch:new{}. " ~
            "A table with input parameters is expected as the first argument.";
        luaL_error(L, errMsg.toStringz);
        return 0;
    }
    if (!checkAllowedNames(L, 1, ["a","centre","face_name","which_part"])) {
        string errMsg = "Error in call to CubePatch:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
        return 0;
    }
    // Now, get our table entries.
    // Table is at index 1
    double a = getDouble(L, 1, "a");
    string face_name = getString(L, 1, "face_name");
    string which_part = getStringWithDefault(L, 1, "which_part", "");
    lua_getfield(L, 1, "centre");
    if (lua_isnil(L, -1)) {
        string errMsg = "Error in call to CubePatch:new{}. No centre entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto centre = toVector3(L, -1);
    lua_pop(L, 1);
    //
    auto spatch = new CubePatch(to!number(a), centre, face_name, which_part);
    surfaceStore ~= pushObj!(CubePatch, CubePatchMT)(L, spatch);
    return 1;
} // end newCubePatch()


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
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected MeshPatch:new{}; ";
        errMsg ~= "maybe you tried MeshPatch.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in constructor MeshPatch:new{}. " ~
            "A table with input parameters is expected as the first argument.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["sgrid"])) {
        string errMsg = "Error in call to MeshPatch:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    // Look for the StructuredGrid object.
    lua_getfield(L, 1, "sgrid");
    auto grid = checkStructuredGrid(L, -1);
    if ( grid is null ) {
        string errMsg = "Error in constructor MeshPatch:new{}. " ~
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
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected LuaFnSurface:new{}; ";
        errMsg ~= "maybe you tried LuaFnSurface.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to LuaFnSurface:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["luaFnName"])) {
        string errMsg = "Error in call to LuaFnSurface:new{}. Invalid name in table.";
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
 * srs = SubRangedSurface:new{underlying_surface=psurf, r0=0.0, r1=1.0, s0=0.0, s1=1.0}
 * --------------------------
 */
extern(C) int newSubRangedSurface(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected SubRangedSurface:new{}; ";
        errMsg ~= "maybe you tried SubRangedSurface.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = `Error in constructor SubRangedSurface:new.
A table with input parameters is expected as the first argument.`;
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["underlying_psurface", "underlying_surface", "r0", "r1", "s0", "s1"])) {
        string errMsg = "Error in call to SubRangedSurface:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect a ParametricSurface object at the first array position in the table.
    // I couldn't decide which field name was clearer and/or more consistent, so I've allowed both.
    lua_getfield(L, 1, "underlying_psurface");
    if (lua_isnil(L, -1)) {
        lua_pop(L, 1); // dispose of the nil value
        lua_getfield(L, 1, "underlying_surface"); // try alternate name
        if (lua_isnil(L, -1)) {
            string errMsg = "Error in call to SubRangedSurface:new{}. No underlying_psurface field found.";
            luaL_error(L, errMsg.toStringz());
        }
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

/**
 * This is the constructor for a BezierPatch to be used from the Lua interface.
 *
 * At successful completion of this function, a new BezierPatch object
 * is pushed onto the Lua stack.
 *
 * Construction is:
 * -------------------------
 * patch = BezierPatch:new{points=Q}
 * --------------------------
 *
 * points  : an array of points Q[n+1][m+1]
 */

extern(C) int newBezierPatch(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 2)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected BezierPatch:new{}; ";
        errMsg ~= "maybe you tried BezierPatch.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in constructor BezierPatch:new{}. " ~
            "A table with input parameters is expected as the first argument.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["points"])) {
        string errMsg = "Error in call to BezierPatch:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_getfield(L, 1, "points");
    if (!lua_istable(L, -1)) {
        string errMsg = "There's a problem in call to BezierPatch:new{}. " ~
            "An array of arrays for 'points' is expected.";
        luaL_error(L, errMsg.toStringz);
    }
    int np1 = to!int(lua_objlen(L, -1));
    int mp1 = 0;
    int n = np1 - 1;
    int m = 0;
    Vector3[][] Q;
    Q.length = np1;
    foreach (i; 0 .. np1) {
        lua_rawgeti(L, -1, i+1); // +1 for Lua offset
        if (i == 0) {
            mp1 = to!int(lua_objlen(L, -1));
            m = mp1 - 1;
        }
        else {
            if (mp1 != lua_objlen(L, -1)) {
                string errMsg = "There's a problem in call to BezierPatch:new{}. " ~
                    "Inconsistent numbers of points in array of 'points'.";
                luaL_error(L, errMsg.toStringz);
            }
        }
        Q[i].length = mp1;
        foreach (j; 0 .. mp1) {
            lua_rawgeti(L, -1, j+1);
            auto p = toVector3(L, -1);
            Q[i][j].set(p);
            lua_pop(L, 1);
        }
        lua_pop(L, 1);
    }
    lua_pop(L, 1);
    // Try to construct object.
    auto bezPatch = new BezierPatch(Q, n, m);
    surfaceStore ~= pushObj!(BezierPatch, BezierPatchMT)(L, bezPatch);
    return 1;
} // end newBezierPatch()

extern(C) int luaWriteCtrlPtsAsVtkXml(lua_State *L)
{
    int narg = lua_gettop(L);
    if (narg < 2) {
        string errMsg = "Not enough arguments passed to writeCtrlPtsVtkXml\n";
        luaL_error(L, errMsg.toStringz);
    }
    auto bPatch = checkObj!(BezierPatch, BezierPatchMT)(L, 1);
    string fName = to!string(luaL_checkstring(L, 2));
    writeCtrlPtsAsVtkXml(bPatch, fName);
    return 0;
}

/**
 * This is the constructor for a ControlPointPatch to be used from the Lua interface.
 *
 * At successful completion of this function, a new ControlPointPatch object
 * is pushed onto the Lua stack.
 *
 * Construction is:
 * -------------------------
 * patch = ControlPointPatch:new{north=nPath, east=ePath, south=sPath, west=wPath, control_points=C}
 * <OR>
 * patch = ControlPointPatch:new{north=nPath, east=ePath, south=sPath, west=wPath,
 *                               ncpi=5, ncpj=3, guide_patch="channel"}
 * --------------------------
 *
 * paths   : see CoonsPatch (above)
 * control_points  : an array of points C[N][M]
 * ncpi : no. control points in i-direction
 * ncpj : no. control points in j-direction
 * guide_patch : a patch type to use for laying down control points
 *
 * @author: Rowan J. Gollan
 */

extern(C) int newControlPointPatch(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 2)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected ControlPointPatch:new{}; ";
        errMsg ~= "maybe you tried ControlPointPatch.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in constructor ControlPointPatch:new{}. " ~
            "A table with input parameters is expected as the first argument.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["north", "east", "south", "west", "control_points",
                                  "ncpi", "ncpj", "guide_patch"])) {
        string errMsg = "Error in call to ControlPointPatch:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }

    Path[string] paths;
    getPaths(L, "ControlPointPatch", paths, true);

    lua_getfield(L, 1, "control_points");
    if (!lua_isnil(L, -1)) {
            // Proceed with constructor based on supplied control points
        if (!lua_istable(L, -1)) {
            string errMsg = "There's a problem in call to ControlPointsPatch:new{}. " ~
                "An array of arrays for 'control_points' is expected.";
            luaL_error(L, errMsg.toStringz);
        }
        int N = to!int(lua_objlen(L, -1));
        int M;
        Vector3[][] C;
        C.length = N;
        foreach (i; 0 .. N) {
            lua_rawgeti(L, -1, i+1); // +1 for Lua offset
            if (i == 0) {
                M = to!int(lua_objlen(L, -1));
            }
            else {
                if (M != lua_objlen(L, -1)) {
                    string errMsg = "There's a problem in call to ControlPointPatch:new{}. " ~
                        "Inconsistent numbers of points in array of 'control_points'.";
                    luaL_error(L, errMsg.toStringz);
                }
            }
            C[i].length = M;
            foreach (j; 0 .. M) {
                lua_rawgeti(L, -1, j+1);
                auto p = toVector3(L, -1);
                C[i][j].set(p);
                lua_pop(L, 1);
            }
            lua_pop(L, 1);
        }
        lua_pop(L, 1);
        // Try to construct object.
        auto ctrlPtPatch = new ControlPointPatch(paths["south"], paths["north"], paths["west"], paths["east"], C);
        surfaceStore ~= pushObj!(ControlPointPatch, ControlPointPatchMT)(L, ctrlPtPatch);
        return 1;
    }
    else {
        // Assume that the constructor based on a guidePatch is at play
        lua_getfield(L, 1, "ncpi");
        int ncpi = to!int(lua_tonumber(L, -1));
        lua_pop(L, 1);
        lua_getfield(L, 1, "ncpj");
        int ncpj = to!int(lua_tonumber(L, -1));
        lua_pop(L, 1);
        lua_getfield(L, 1, "guide_patch");
        string guidePatch = to!string(lua_tostring(L, -1));
        lua_pop(L, 1);
        // Try to construct object.
        auto ctrlPtPatch = new ControlPointPatch(paths["south"], paths["north"], paths["west"], paths["east"], ncpi, ncpj, guidePatch);
        surfaceStore ~= pushObj!(ControlPointPatch, ControlPointPatchMT)(L, ctrlPtPatch);
        return 1;
    }



} // end newControlPointPatch()

extern(C) int luaSetCtrlPtInPatch(lua_State *L)
{
    int narg = lua_gettop(L);
    if (narg < 4) {
        string errMsg = "Not enough arguments passed to ControlPointPatch:setCtrlPt\n";
        luaL_error(L, errMsg.toStringz);
    }
    auto ctrlPatch = checkObj!(ControlPointPatch, ControlPointPatchMT)(L, 1);
    int i = to!int(luaL_checkinteger(L, 2));
    int j = to!int(luaL_checkinteger(L, 3));
    auto p = toVector3(L, 4);
    ctrlPatch.setCtrlPt(i, j, p);
    return 0;
}

extern(C) int luaGetCtrlPts(lua_State *L)
{
    auto ctrlPatch = checkObj!(ControlPointPatch, ControlPointPatchMT)(L, 1);
    int nci, ncj;
    ctrlPatch.nCtrlPts(nci, ncj);
    lua_newtable(L);
    foreach (i; 0 .. nci) {
        lua_newtable(L);
        lua_rawseti(L, -2, i);
    }
    foreach (i; 0 .. nci) {
        lua_rawgeti(L, -1, i);
        foreach (j; 0 .. ncj) {
            pushVector3(L, ctrlPatch.getCtrlPt(i,j));
            lua_rawseti(L, -2, j);
        }
        lua_pop(L, 1);
    }
    return 1;
}

extern(C) int luaWriteCtrlPtsInPatchAsVtkXml(lua_State *L)
{
    int narg = lua_gettop(L);
    if (narg < 2) {
        string errMsg = "Not enough arguments passed to ControlPointPatch:writeCtrlPtsAsVtkXml\n";
        luaL_error(L, errMsg.toStringz);
    }
    auto ctrlPatch = checkObj!(ControlPointPatch, ControlPointPatchMT)(L, 1);
    string bfn = to!string(luaL_checkstring(L, 2));
    ctrlPatch.writeCtrlPtsAsVtkXml(bfn);
    return 0;
}
/**
 * This is the constructor for a NURBSPatch to be used from the Lua interface.
 *
 * At successful completion of this function, a new NURBSPatch object
 * is pushed onto the Lua stack.
 *
 * Construction is:
 * -------------------------
 * patch = NURBSPatch:new{points=P, weights=w,
                          u_knots=U, v_knots=V,
			  u_degree=p, v_degree=q}
 * --------------------------
 *
 * [TODO] -- document inputs
 */

extern(C) int newNURBSPatch(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 2)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected NURBSPatch:new{}; ";
        errMsg ~= "maybe you tried NURBSPatch.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in constructor NURBSPatch:new{}. " ~
            "A table with input parameters is expected as the first argument.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["points", "weights",
				  "u_knots", "v_knots",
				  "u_degree", "v_degree"])) {
        string errMsg = "Error in call to NURBSPatch:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }

    lua_getfield(L, 1, "points");
    if (!lua_istable(L, -1)) {
        string errMsg = "There's a problem in call to NURBSPatch:new{}. " ~
            "An array of arrays for 'points' is expected.";
        luaL_error(L, errMsg.toStringz);
    }
    int np1 = to!int(lua_objlen(L, -1));
    int mp1 = 0;
    int n = np1 - 1;
    int m = 0;
    Vector3[][] P;
    P.length = np1;
    foreach (i; 0 .. np1) {
        lua_rawgeti(L, -1, i+1); // +1 for Lua offset
        if (i == 0) {
            mp1 = to!int(lua_objlen(L, -1));
            m = mp1 - 1;
        }
        else {
            if (mp1 != lua_objlen(L, -1)) {
                string errMsg = "There's a problem in call to NURBSPatch:new{}. " ~
                    "Inconsistent numbers of points in array of 'points'.";
                luaL_error(L, errMsg.toStringz);
            }
        }
        P[i].length = mp1;
        foreach (j; 0 .. mp1) {
            lua_rawgeti(L, -1, j+1);
            auto p = toVector3(L, -1);
            P[i][j].set(p);
            lua_pop(L, 1);
        }
        lua_pop(L, 1);
    }
    lua_pop(L, 1);

    // Get "weights"
    double[4][][] Pw;
    double[] w;
    lua_getfield(L, 1, "weights");
    if (!lua_istable(L, -1)) {
        string errMsg = "There's a problem in call to NURBSPatch:new{}. " ~
            "An array of arrays for 'weights' is expected.";
        luaL_error(L, errMsg.toStringz);
    }
    int nw1 = to!int(lua_objlen(L, -1));
    Pw.length = nw1;
    foreach (i; 0 .. nw1) {
	lua_rawgeti(L, -1, i+1); // +1 for lua offset
	int mw = to!int(lua_objlen(L, -1));
	w.length = 0;
	foreach (j; 0 .. mw) {
	    lua_rawgeti(L, -1, j+1); // +1 for lua offset
	    if (lua_isnumber(L, -1)) w ~= lua_tonumber(L, -1);
	    // Silently ignore anything else
	    lua_pop(L, 1);
	}
	lua_pop(L, 1);
	Pw[i].length = w.length;
	foreach (j, wt; w) {
	    Pw[i][j][0] = P[i][j].x.re * wt;
	    Pw[i][j][1] = P[i][j].y.re * wt;
	    Pw[i][j][2] = P[i][j].z.re * wt;
	    Pw[i][j][3] = wt;
	}
    }
    lua_pop(L, 1);

    // Get "u_knots"
    double[] U;
    getArrayOfDoubles(L, 1, "u_knots", U);

    // Get "v_knots"
    double[] V;
    getArrayOfDoubles(L, 1, "v_knots", V);

    // Get "u_degree"
    int p = getInt(L, 1, "u_degree");

    // Get "v_degree"
    int q = getInt(L, 1, "v_degree");

    // Try to construct object.
    auto nurbsPatch = new NURBSSurface(Pw, U, p, V, q);
    surfaceStore ~= pushObj!(NURBSSurface, NURBSPatchMT)(L, nurbsPatch);
    return 1;
} // end newNURBSPatch()



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
    if (!checkAllowedNames(L, 1, ["north", "south", "west", "east", "gridType"])) {
        string errMsg = "Error in call to makePatch:new{}. Invalid name in table.";
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


/**
 * This is the constructor for a BezierTrianglePatch to be used from the Lua interface.
 *
 * At successful completion of this function, a new BezierTrianglePatch object
 * is pushed onto the Lua stack.
 *
 * Construction is:
 * -------------------------
 * patch = BezierTrianglePatch:new{points=Q, degree=n}
 * --------------------------
 *
 * points  : a linear array of points Q
 * degree  : degree of patch
 *
 * A degree n Bezier triangle has (n+1)(n+2)/2 points.
 * The points are given in linear order. An example for
 * n=3 is:
 *                     + b_(0,3,0)
 *
 *              + b_(0,2,1)   + b_(1,2,0)
 *
 *       + b_(0,1,2)   + b(1,1,1)    + b(2,1,0)
 *
 * + b_(0,0,3)  + b_(1,0,2)  + b_(2,0,1)    + b_(3,0,0)
 *
 * Array entries are:
 * Q = [b_(3,0,0), b_(2,1,0), b_(2,0,1), b_(1,2,0), b_(1,1,1),
 *      b_(1,0,2), b_(0,3,0), b_(0,2,1), b_(0,1,2), b_(0,0,3)]
 *
 * As a loop, the entries are:
 *
 *    for i=n,0,-1 do
 *        for j=n-i,0,-1 do
 *            k = n - (i+j)
 *            pts[#pts+1] = b_(i,j,k)
 *        end
 *    end
 */

extern(C) int newBezierTrianglePatch(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 2)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Excepted BezierTrianglePatch:new{};\n";
        errMsg ~= "Maybe you tried BezierTrianglePatch.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if (!checkAllowedNames(L, 1, ["points", "degree"])) {
        string errMsg = "Error in call to BezierTrianglePatch:new{}.\n";
        errMsg ~= "Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_getfield(L, 1, "degree");
    int n = to!int(luaL_checkint(L, -1));
    lua_pop(L, 1);
    lua_getfield(L, 1, "points");
    if (!lua_istable(L, -1)) {
        string errMsg = "Error in call to BezierTrianglePatch:new{}.\n";
        errMsg ~= "An array is expected for 'points'.\n";
        luaL_error(L, errMsg.toStringz);
    }
    int nPtsExpected = (n+1)*(n+2)/2;
    int nPts = to!int(lua_objlen(L, -1));
    if (nPts != nPtsExpected) {
        string errMsg = "Error in call to BezierTrianglePatch:new{}.\n";
        errMsg ~= "The number of supplied points is not consistent with the supplied degree of patch.\n";
        errMsg ~= format("The degree is %d, so %d points are expected.\n", n, nPtsExpected);
        errMsg ~= format("The number of points found is %d.\n", nPts);
        luaL_error(L, errMsg.toStringz);
    }
    Vector3[] Q;
    foreach (i; 1 .. nPts+1) {
        lua_rawgeti(L, -1, i);
        Q ~= toVector3(L, -1);
        lua_pop(L, 1);
    }
    lua_pop(L, 1);
    // Now try to construct object, and keep a D copy in store.
    auto bezTriPatch = new BezierTrianglePatch(Q, n);
    triPatchStore ~= pushObj!(BezierTrianglePatch, BezierTrianglePatchMT)(L, bezTriPatch);
    return 1;
}

BezierTrianglePatch checkBezierTrianglePatch(lua_State* L, int index)
{
    if (isObjType(L, index, BezierTrianglePatchMT))
        return checkObj!(BezierTrianglePatch, BezierTrianglePatchMT)(L, index);
    // no match found
    return null;
}

extern(C) int opCallBezierTrianglePatch(lua_State* L)
{
    auto bezTri = checkBezierTrianglePatch(L, 1);
    auto u = luaL_checknumber(L, 2);
    auto v = luaL_checknumber(L, 3);
    auto p = bezTri(u, v);
    return pushVector3(L, p);
}

extern(C) int writeBezierTriangleCtrlPtsAsText(lua_State* L)
{
    int narg = lua_gettop(L);
    if (narg < 2) {
        string errMsg = "Not enough arguments passed to writeBezierTriangleCtrlPtsAsText().\n";
        luaL_error(L, errMsg.toStringz);
    }
    auto bezTri = checkBezierTrianglePatch(L, 1);
    string fName = to!string(luaL_checkstring(L, 2));
    geom.writeBezierTriangleCtrlPtsAsText(bezTri, fName);
    return 0;
}

extern(C) int writeBezierTriangleCtrlPtsAsVtkXml(lua_State* L)
{
    int narg = lua_gettop(L);
    if (narg < 2) {
        string errMsg = "Not enough arguments passed to writeBezierTriangleCtrlPtsAsVtkXml().\n";
        luaL_error(L, errMsg.toStringz);
    }
    auto bezTri = checkBezierTrianglePatch(L, 1);
    string fName = to!string(luaL_checkstring(L, 2));
    geom.writeBezierTriangleCtrlPtsAsVtkXml(bezTri, fName);
    return 0;
}

extern(C) int writeBezierTriangleAsVtkXml(lua_State* L)
{
    int narg = lua_gettop(L);
    if (narg < 3) {
        string errMsg = "Not enough arguments passed to writeBezierTriangleCtrlPtsAsVtkXml().\n";
        luaL_error(L, errMsg.toStringz);
    }
    auto bezTri = checkBezierTrianglePatch(L, 1);
    string fName = to!string(luaL_checkstring(L, 2));
    int nEdgePts = luaL_checkint(L, 3);
    geom.writeBezierTriangleAsVtkXml(bezTri, fName, nEdgePts);
    return 0;
}

extern(C) int writeBezierTriangleAsDat(lua_State* L)
{
    int narg = lua_gettop(L);
    if (narg < 3) {
        string errMsg = "Not enough arguments passed to writeBezierTriangleCtrlPtsAsDat().\n";
        luaL_error(L, errMsg.toStringz);
    }
    auto bezTri = checkBezierTrianglePatch(L, 1);
    string fName = to!string(luaL_checkstring(L, 2));
    int nEdgePts = luaL_checkint(L, 3);
    geom.writeBezierTriangleAsDat(bezTri, fName, nEdgePts);
    return 0;
}

extern(C) int bezierTriangleFromPointCloud(lua_State* L)
{
    int narg = lua_gettop(L);
    if (narg < 5) {
        string errMsg = "Not enough arguments passed to bezierTriangleFromPointCloud().\n";
        luaL_error(L, errMsg.toStringz);
    }
    if (!lua_istable(L, 1)) {
        string errMsg = "An array of points is expected as first argument in\n";
        errMsg ~= "bezierTriangleFromPointCloud().\n";
        errMsg ~= "No points array found.\n";
        luaL_error(L, errMsg.toStringz);
    }
    int nPts = to!int(lua_objlen(L, 1));
    Vector3[] pts;
    foreach (i; 1 .. nPts+1) {
        lua_rawgeti(L, 1, i);
        pts ~= toVector3(L, -1);
        lua_pop(L, 1);
    }
    auto b0 = checkObj!(Bezier, BezierMT)(L, 2);
    auto b1 = checkObj!(Bezier, BezierMT)(L, 3);
    auto b2 = checkObj!(Bezier, BezierMT)(L, 4);
    int n = luaL_checkint(L, 5);
    BezierTrianglePatch initGuess;
    if (narg >= 6) {
        // Attempt to grab a BezierTriangle
        initGuess = checkBezierTrianglePatch(L, 6);
    }
    bool fittingSuccess;
    auto bezTri = geom.bezierTriangleFromPointCloud(pts, b0, b1, b2, n, initGuess, fittingSuccess);
    triPatchStore ~= pushObj!(BezierTrianglePatch, BezierTrianglePatchMT)(L, bezTri);
    lua_pushboolean(L, fittingSuccess);
    return 2;
}

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
    lua_pushcfunction(L, &areaOfSurface!(CoonsPatch, CoonsPatchMT));
    lua_setfield(L, -2, "area");

    lua_setglobal(L, CoonsPatchMT.toStringz);
    lua_getglobal(L, CoonsPatchMT.toStringz); lua_setglobal(L, "CoonsSurface"); // alias

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
    lua_pushcfunction(L, &areaOfSurface!(AOPatch, AOPatchMT));
    lua_setfield(L, -2, "area");

    lua_setglobal(L, AOPatchMT.toStringz);
    lua_getglobal(L, AOPatchMT.toStringz); lua_setglobal(L, "AOSurface"); // alias

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
    lua_pushcfunction(L, &areaOfSurface!(ChannelPatch, ChannelPatchMT));
    lua_setfield(L, -2, "area");
    lua_pushcfunction(L, &make_bridging_path_ChannelPatch);
    lua_setfield(L, -2, "make_bridging_path");

    lua_setglobal(L, ChannelPatchMT.toStringz);
    lua_getglobal(L, ChannelPatchMT.toStringz); lua_setglobal(L, "ChannelSurface"); // alias

    // Register the RuledSurface object
    luaL_newmetatable(L, RuledSurfaceMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newRuledSurface);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallSurface!(RuledSurface, RuledSurfaceMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallSurface!(RuledSurface, RuledSurfaceMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(RuledSurface, RuledSurfaceMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &areaOfSurface!(RuledSurface, RuledSurfaceMT));
    lua_setfield(L, -2, "area");

    lua_setglobal(L, RuledSurfaceMT.toStringz);
    lua_getglobal(L, RuledSurfaceMT.toStringz);

    // Register the NozzleExpansionpatch object
    luaL_newmetatable(L, NozzleExpansionPatchMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newNozzleExpansionPatch);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallSurface!(NozzleExpansionPatch, NozzleExpansionPatchMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallSurface!(NozzleExpansionPatch, NozzleExpansionPatchMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(NozzleExpansionPatch, NozzleExpansionPatchMT));
    lua_setfield(L, -2, "__tostring");

    lua_setglobal(L, NozzleExpansionPatchMT.toStringz);
    lua_getglobal(L, NozzleExpansionPatchMT.toStringz); lua_setglobal(L, "ExpandingChannelPatch"); // alias

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
    lua_getglobal(L, SweptPathPatchMT.toStringz); lua_setglobal(L, "SweptPathSurface"); // alias

    // Register the SpherePatch object
    luaL_newmetatable(L, SpherePatchMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newSpherePatch);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallSurface!(SpherePatch, SpherePatchMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallSurface!(SpherePatch, SpherePatchMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(SpherePatch, SpherePatchMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &areaOfSurface!(SpherePatch, SpherePatchMT));
    lua_setfield(L, -2, "area");

    lua_setglobal(L, SpherePatchMT.toStringz);
    lua_getglobal(L, SpherePatchMT.toStringz); lua_setglobal(L, "SphereSurface"); // alias

    // Register the CubePatch object
    luaL_newmetatable(L, CubePatchMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newCubePatch);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallSurface!(CubePatch, CubePatchMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallSurface!(CubePatch, CubePatchMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(CubePatch, CubePatchMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &areaOfSurface!(CubePatch, CubePatchMT));
    lua_setfield(L, -2, "area");

    lua_setglobal(L, CubePatchMT.toStringz);
    lua_getglobal(L, CubePatchMT.toStringz); lua_setglobal(L, "CubeSurface"); // alias

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
    lua_getglobal(L, MeshPatchMT.toStringz); lua_setglobal(L, "MeshSurface"); // alias

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
    lua_pushcfunction(L, &areaOfSurface!(LuaFnSurface, LuaFnSurfaceMT));
    lua_setfield(L, -2, "area");

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
    lua_pushcfunction(L, &areaOfSurface!(SubRangedSurface, SubRangedSurfaceMT));
    lua_setfield(L, -2, "area");

    lua_setglobal(L, SubRangedSurfaceMT.toStringz);

    // Register the BezierPatch object
    luaL_newmetatable(L, BezierPatchMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newBezierPatch);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallSurface!(BezierPatch, BezierPatchMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallSurface!(BezierPatch, BezierPatchMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(BezierPatch, BezierPatchMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &areaOfSurface!(BezierPatch, BezierPatchMT));
    lua_setfield(L, -2, "area");

    lua_setglobal(L, BezierPatchMT.toStringz);
    lua_getglobal(L, BezierPatchMT.toStringz); lua_setglobal(L, "BezierSurface"); // alias

    // Register function to write Bezier control points.
    lua_pushcfunction(L, &luaWriteCtrlPtsAsVtkXml); lua_setglobal(L, "writeCtrlPtsAsVtkXml");

    // Register the ControlPointPatch object
    luaL_newmetatable(L, ControlPointPatchMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newControlPointPatch);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallSurface!(ControlPointPatch, ControlPointPatchMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallSurface!(ControlPointPatch, ControlPointPatchMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(ControlPointPatch, ControlPointPatchMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &areaOfSurface!(ControlPointPatch, ControlPointPatchMT));
    lua_setfield(L, -2, "area");
    lua_pushcfunction(L, &luaGetCtrlPts);
    lua_setfield(L, -2, "getCtrlPts");
    lua_pushcfunction(L, &luaSetCtrlPtInPatch);
    lua_setfield(L, -2, "setCtrlPt");
    lua_pushcfunction(L, &luaWriteCtrlPtsInPatchAsVtkXml);
    lua_setfield(L, -2, "writeCtrlPtsAsVtkXml");

    lua_setglobal(L, ControlPointPatchMT.toStringz);
    lua_getglobal(L, ControlPointPatchMT.toStringz); lua_setglobal(L, "ControlPointSurface"); // alias


    // Register the NURBSPatch object
    luaL_newmetatable(L, NURBSPatchMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newNURBSPatch);
    lua_setfield(L, -2, "new");

    lua_pushcfunction(L, &opCallSurface!(NURBSSurface, NURBSPatchMT));
    lua_setfield(L, -2, "__call");

    lua_pushcfunction(L, &opCallSurface!(NURBSSurface, NURBSPatchMT));
    lua_setfield(L, -2, "eval");

    lua_pushcfunction(L, &toStringObj!(NURBSSurface, NURBSPatchMT));
    lua_setfield(L, -2, "__tostring");

    lua_pushcfunction(L, &areaOfSurface!(NURBSSurface, NURBSPatchMT));
    lua_setfield(L, -2, "area");

    lua_setglobal(L, NURBSPatchMT.toStringz);
    lua_getglobal(L, NURBSPatchMT.toStringz); lua_setglobal(L, "NURBSSurface"); // alias


    // Register utility functions.
    lua_pushcfunction(L, &isSurface); lua_setglobal(L, "isSurface");
    lua_pushcfunction(L, &makePatch); lua_setglobal(L, "makePatch");
    lua_pushcfunction(L, &makePatch); lua_setglobal(L, "makeSurface"); // alias
    lua_pushcfunction(L, &luaWriteSurfaceAsVtkXml); lua_setglobal(L, "writeSurfaceAsVtkXml");

    // Register the BezierTrianglePatch object
    luaL_newmetatable(L, BezierTrianglePatchMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newBezierTrianglePatch);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallBezierTrianglePatch);
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallBezierTrianglePatch);
    lua_setfield(L, -2, "eval");

    lua_setglobal(L, BezierTrianglePatchMT.toStringz);
    lua_getglobal(L, BezierTrianglePatchMT.toStringz); lua_setglobal(L, "BezierTriangleSurface"); // alias

    lua_pushcfunction(L, &writeBezierTriangleCtrlPtsAsText);
    lua_setglobal(L, "writeBezierTriangleCtrlPtsAsText");
    lua_pushcfunction(L, &writeBezierTriangleCtrlPtsAsVtkXml);
    lua_setglobal(L, "writeBezierTriangleCtrlPtsAsVtkXml");
    lua_pushcfunction(L, &writeBezierTriangleAsVtkXml);
    lua_setglobal(L, "writeBezierTriangleAsVtkXml");
    lua_pushcfunction(L, &writeBezierTriangleAsDat);
    lua_setglobal(L, "writeBezierTriangleAsDat");
    lua_pushcfunction(L, &bezierTriangleFromPointCloud);
    lua_setglobal(L, "bezierTriangleFromPointCloud");

} // end registerSurfaces()
