/**
 * A Lua interface for the D gpath module.
 *
 * Authors: Rowan G. and Peter J.
 * Date: 2015-02-22 First code
 *       2015-04-22 greatly expanded with Arc, Bezier, Polyline, LuaFnPath
 */

module geom.luawrap.luagpath;

import util.lua;
import std.stdio;
import std.string;
import std.conv;
import util.lua_service;
import geom;
import geom.luawrap.luageom;

immutable string LineMT = "Line"; // Name of Line metatable
immutable string ArcMT = "Arc";
immutable string Arc3MT = "Arc3";
immutable string HelixMT = "Helix";
immutable string BezierMT = "Bezier";
immutable string NURBSMT = "NURBS";
immutable string PolylineMT = "Polyline";
immutable string SplineMT = "Spline";
immutable string Spline2MT = "Spline2";
immutable string XSplineMT = "XSpline";
immutable string XSpline2MT = "XSpline2";
immutable string XSplineLsqMT = "XSplineLsq";
immutable string XSplineLsq2MT = "XSplineLsq2";
immutable string SVGPathMT = "SVGPath";
immutable string LuaFnPathMT = "LuaFnPath";
immutable string ArcLengthParameterizedPathMT = "ArcLengthParameterizedPath";
immutable string SubRangedPathMT = "SubRangedPath";
immutable string ReversedPathMT = "ReversedPath";
immutable string TranslatedPathMT = "TranslatedPath";
immutable string MirrorImagePathMT = "MirrorImagePath";
immutable string RotatedAboutZAxisPathMT = "RotatedAboutZAxisPath";

// A place to hang on to references to objects that are pushed into the Lua domain.
// We don't want the D garbage collector to prematurely dispose of said objects.
static const(Path)[] pathStore;

Path checkPath(lua_State* L, int index) {
    if ( isObjType(L, index, LineMT) ) {
        return checkObj!(Line, LineMT)(L, index);
    }
    if ( isObjType(L, index, ArcMT) ) {
        return checkObj!(Arc, ArcMT)(L, index);
    }
    if ( isObjType(L, index, Arc3MT) ) {
        return checkObj!(Arc3, Arc3MT)(L, index);
    }
    if ( isObjType(L, index, HelixMT) ) {
        return checkObj!(Helix, HelixMT)(L, index);
    }
    if ( isObjType(L, index, BezierMT) ) {
        return checkObj!(Bezier, BezierMT)(L, index);
    }
    if ( isObjType(L, index, NURBSMT) ) {
        return checkObj!(NURBS, NURBSMT)(L, index);
    }
    if ( isObjType(L, index, PolylineMT) ) {
        return checkObj!(Polyline, PolylineMT)(L, index);
    }
    if ( isObjType(L, index, XSplineMT) ) {
        return checkObj!(XSpline, XSplineMT)(L, index);
    }
    if ( isObjType(L, index, XSplineLsqMT) ) {
        return checkObj!(XSplineLsq, XSplineLsqMT)(L, index);
    }
    if ( isObjType(L, index, SVGPathMT) ) {
        return checkObj!(SVGPath, SVGPathMT)(L, index);
    }
    if ( isObjType(L, index, LuaFnPathMT) ) {
        return checkObj!(LuaFnPath, LuaFnPathMT)(L, index);
    }
    if ( isObjType(L, index, ArcLengthParameterizedPathMT) ) {
        return checkObj!(ArcLengthParameterizedPath,
                         ArcLengthParameterizedPathMT)(L, index);
    }
    if ( isObjType(L, index, SubRangedPathMT) ) {
        return checkObj!(SubRangedPath, SubRangedPathMT)(L, index);
    }
    if ( isObjType(L, index, ReversedPathMT) ) {
        return checkObj!(ReversedPath, ReversedPathMT)(L, index);
    }
    if ( isObjType(L, index, TranslatedPathMT) ) {
        return checkObj!(TranslatedPath, TranslatedPathMT)(L, index);
    }
    if ( isObjType(L, index, MirrorImagePathMT) ) {
        return checkObj!(MirrorImagePath, MirrorImagePathMT)(L, index);
    }
    if ( isObjType(L, index, RotatedAboutZAxisPathMT) ) {
        return checkObj!(RotatedAboutZAxisPath, RotatedAboutZAxisPathMT)(L, index);
    }
    // if all else fails
    return null;
}

extern(C) int opCallPath(T, string MTname)(lua_State* L)
{
    auto path = checkObj!(T, MTname)(L, 1);
    auto t = luaL_checknumber(L, 2);
    auto pt = path(t);
    return pushVector3(L, pt);
}

extern(C) int copyPath(T, string MTname)(lua_State* L)
{
    // Sometimes it's convenient to get a copy of a path.
    auto path = checkObj!(T, MTname)(L, 1);
    pathStore ~= pushObj!(T, MTname)(L, path.dup()); // new object
    return 1;
}

extern(C) int pathIntersect2D(T, string MTname)(lua_State* L)
// Example of use:
// found, t = path:intersect2D{ps=pstart, d=mydir, nseg=20}
{
    auto path = checkObj!(T, MTname)(L, 1);
    int narg = lua_gettop(L);
    if ( narg == 1 || !lua_istable(L, 2) ) {
        string errMsg = "Error in call to Path:intersect2D{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 2, ["ps", "d", "nseg"])) {
        string errMsg = "Error in call to Path:intersect2D{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect Vector3 for starting point.
    lua_getfield(L, 2, "ps");
    if (lua_isnil(L, -1)) {
        string errMsg = "Error in call to Path:intersect2D{}. No ps entry found." ~
            "Check that the keyword argument 'ps' is present,\n" ~
            "and that a valid object is passed as value.\n";
        luaL_error(L, errMsg.toStringz());
    }
    auto ps = toVector3(L, -1);
    lua_pop(L, 1);
    // Expect Vector3 for direction vector.
    lua_getfield(L, 2, "d");
    if (lua_isnil(L, -1)) {
        string errMsg = "Error in call to Path:intersect2D{}. No d entry found.\n" ~
            "Check that the keyword argument 'd' is present,\n" ~
            "and that a valid object is passed as value.";
        luaL_error(L, errMsg.toStringz());
    }
    auto d = toVector3(L, -1);
    lua_pop(L, 1);
    int nseg = 20; // default value
    lua_getfield(L, 2, "nseg");
    if (lua_isnumber(L, -1)) {
        nseg = to!int(lua_tonumber(L, -1));
    }
    lua_pop(L, 1);
    //
    double t = 0.0;
    bool found = path.intersect2D(ps, d, t, nseg);
    lua_settop(L, 0); // clear stack
    lua_pushboolean(L, found);
    lua_pushnumber(L, t);
    return 2;
} // end pathIntersect2D()()

extern(C) int pathClosestDistance(T, string MTname)(lua_State *L)
{
    auto path = checkObj!(T, MTname)(L, 1);

    int narg = lua_gettop(L);
    if (narg == 1) {
        string errMsg = "Error in call to Path:closestDistance().; " ~
            "A vector as argument is expected, but no arguments found.";
        luaL_error(L, errMsg.toStringz);
    }
    auto p = toVector3(L, 2);

    double t = 0.0;
    double dist = path.closestDistance(p, t);
    lua_settop(L, 0); //clear stack
    lua_pushnumber(L, dist);
    lua_pushnumber(L, t);
    return 2;
}

extern(C) int length(T, string MTname)(lua_State *L)
{
    auto path = checkObj!(T, MTname)(L, 1);

    auto length = path.length();
    lua_settop(L, 0); // clear stack
    lua_pushnumber(L, length.re);
    return 1;
}

/* ----------------- Specific constructors --------------- */

/**
 * The Lua constructor for a Line.
 *
 * Example construction in Lua:
 * ---------------------------------
 * a = Vector3:new{}
 * b = Vector3:new{x=1, y=1}
 * ab = Line:new{p0=a, p1=b}
 * ---------------------------------
 */
extern(C) int newLine(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected Line:new{p0=a, p1=b}; ";
        errMsg ~= "maybe you tried Line.new{p0=a, p1=b}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to Line:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["p0", "p1"])) {
        string errMsg = "Error in call to Line:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect Vector3 for starting point.
    lua_getfield(L, 1, "p0");
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to Line:new{}. No p0 entry found." ~
            "Check that the keyword argument 'p0' is present,\n" ~
            "and that a valid object is passed as value.\n";
        luaL_error(L, errMsg.toStringz());
    }
    auto p0 = toVector3(L, -1);
    lua_pop(L, 1);
    // Expect Vector3 for end point.
    lua_getfield(L, 1, "p1");
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to Line:new{}. No p1 entry found.\n" ~
            "Check that the keyword argument 'p1' is present,\n" ~
            "and that a valid object is passed as value.";
        luaL_error(L, errMsg.toStringz());
    }
    auto p1 = toVector3(L, -1);
    lua_pop(L, 1);
    auto my_line = new Line(p0, p1);
    pathStore ~= pushObj!(Line, LineMT)(L, my_line);
    return 1;
} // end newLine()


/**
 * The Lua constructor for an Arc.
 *
 * Example construction in Lua:
 * ---------------------------------
 * a = Vector3:new{x=1}
 * b = Vector3:new{y=1}
 * c = Vector3:new{}
 * abc = Arc:new{p0=a, p1=b, centre=c}
 * ---------------------------------
 */
extern(C) int newArc(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected Arc:new{p0=a, p1=b, centre=c}; ";
        errMsg ~= "maybe you tried Arc.new{p0=a, p1=b, centre=c}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to Arc:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["p0", "p1", "centre"])) {
        string errMsg = "Error in call to Arc:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect Vector3 for starting point.
    lua_getfield(L, 1, "p0");
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to Arc:new{}. No p0 entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto p0 = toVector3(L, -1);
    lua_pop(L, 1);
    // Expect Vector3 for end point.
    lua_getfield(L, 1, "p1");
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to Arc:new{}. No p1 entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto p1 = toVector3(L, -1);
    lua_pop(L, 1);
    // Expect Vector3 at centre.
    lua_getfield(L, 1, "centre");
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to Arc:new{}. No centre entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto centre = toVector3(L, -1);
    lua_pop(L, 1);
    auto my_arc = new Arc(p0, p1, centre);
    pathStore ~= pushObj!(Arc, ArcMT)(L, my_arc);
    return 1;
} // end newArc()


/**
 * The Lua constructor for an Arc3.
 *
 * Example construction in Lua:
 * ---------------------------------
 * a = Vector3:new{x=1}
 * m = Vector3:new{x=0.707107, y=0.707107}
 * b = Vector3:new{y=1}
 * amb = Arc3:new{p0=a, pmid=m, p1=b}
 * ---------------------------------
 */
extern(C) int newArc3(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected Arc3:new{p0=a, pmid=m, p1=b}; ";
        errMsg ~= "maybe you tried Arc3.new{p0=a, pmid=m, p1=b}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to Arc3:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["p0", "pmid", "p1"])) {
        string errMsg = "Error in call to Arc3:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect Vector3 for start point p0.
    lua_getfield(L, 1, "p0");
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to Arc3:new{}. No p0 entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto p0 = toVector3(L, -1);
    lua_pop(L, 1);
    // Expect Vector3 at mid-point pmid.
    lua_getfield(L, 1, "pmid");
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to Arc3:new{}. No pmid entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto pmid = toVector3(L, -1);
    lua_pop(L, 1);
    // Expect Vector3 at end point p1.
    lua_getfield(L, 1, "p1");
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to Arc3:new{}. No p1 entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto p1 = toVector3(L, -1);
    lua_pop(L, 1);
    auto my_arc = new Arc3(p0, pmid, p1);
    pathStore ~= pushObj!(Arc3, Arc3MT)(L, my_arc);
    return 1;
} // end newArc3()


/**
 * The Lua constructor for a Helix.
 *
 * Example construction in Lua:
 * ---------------------------------
 * axis0 = Vector3:new{x=0}
 * axis1 = Vector3:new{x=1}
 * pstart = Vector3:new{y=1}
 * pend = Vector3:new{x=1, z=1}
 * h1 = Helix:new{point_start=pstart, point_end=pend, axis0=axis0, axis1=axis1}
 * h2 = Helix:new{a0=Vector3:new{x=0.0}, a1=Vector3:new{x=1.0},
 *                xlocal=Vector3:new{y=1.0},
 *                r0=1.0, r1=1.0, dtheta=math.pi/2};
 * ---------------------------------
 */
extern(C) int newHelix(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected Helix:new{table-of-args}; ";
        errMsg ~= "maybe you tried Helix.new{table-of-args}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to Helix:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["a0", "a1", "xlocal", "r0", "r1", "dtheta",
                                  "point_start", "point_end", "axis0", "axis1"])) {
        string errMsg = "Error in call to Helix:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    // There are two ways to specify the Helix:
    // (a) with the fundamental parameters defining the local axis, radii and angles
    // (b) with start and end points about an axis.
    Helix my_helix;
    lua_getfield(L, 1, "a0");
    if (lua_isnil(L, -1)) {
        // Since we did not find the beginning point on the local axis,
        // we assume that the specification is with start and end points
        // about an axis.
        lua_pop(L, 1);
        // Expect Vector3 for start point.
        lua_getfield(L, 1, "point_start");
        auto pstart = toVector3(L, -1);
        lua_pop(L, 1);
        // Expect Vector3 at end point.
        lua_getfield(L, 1, "point_end");
        if (lua_isnil(L, -1)) {
            string errMsg = "Error in call to Helix:new{}. No point_end entry found.";
            luaL_error(L, errMsg.toStringz());
        }
        auto pend = toVector3(L, -1);
        lua_pop(L, 1);
        // Expect Vector3 specifying start of axis.
        lua_getfield(L, 1, "axis0");
        if (lua_isnil(L, -1)) {
            string errMsg = "Error in call to Helix:new{}. No axis0 entry found.";
            luaL_error(L, errMsg.toStringz());
        }
        auto axis0 = toVector3(L, -1);
        lua_pop(L, 1);
        // Expect Vector3 at end of axis, axis1.
        lua_getfield(L, 1, "axis1");
        if (lua_isnil(L, -1)) {
            string errMsg = "Error in call to Helix:new{}. No axis1 entry found.";
            luaL_error(L, errMsg.toStringz());
        }
        auto axis1 = toVector3(L, -1);
        lua_pop(L, 1);
        my_helix = new Helix(pstart, pend, axis0, axis1);
    } else {
        // Proceed with the specification using fundamental parameters.
        // Expect Vector3 for start point on local axis, a0.
        auto a0 = toVector3(L, -1);
        lua_pop(L, 1);
        // Expect Vector3 at end point on local axis, a1.
        lua_getfield(L, 1, "a1");
        if (lua_isnil(L, -1)) {
            string errMsg = "Error in call to Helix:new{}. No a1 entry found.";
            luaL_error(L, errMsg.toStringz());
        }
        auto a1 = toVector3(L, -1);
        lua_pop(L, 1);
        // Expect Vector3 specifying local x-direction.
        lua_getfield(L, 1, "xlocal");
        if (lua_isnil(L, -1)) {
            string errMsg = "Error in call to Helix:new{}. No xlocal entry found.";
            luaL_error(L, errMsg.toStringz());
        }
        auto xlocal = toVector3(L, -1);
        lua_pop(L, 1);
        string errMsgTmplt = "Error in call to Helix:new{}. " ~
            "A valid value for '%s' was not found in list of arguments. " ~
            "The value should be a number.";
        double r0 = getNumberFromTable(L, 1, "r0", true, 1.0, true, format(errMsgTmplt, "r0"));
        double r1 = getNumberFromTable(L, 1, "r1", true, 1.0, true, format(errMsgTmplt, "r1"));
        double dtheta = getNumberFromTable(L, 1, "dtheta", true, 0.0, true, format(errMsgTmplt, "dtheta"));
        my_helix = new Helix(a0, a1, xlocal, r0, r1, dtheta);
    }
    assert(my_helix !is null, "Did not successfully make a Helix object.");
    pathStore ~= pushObj!(Helix, HelixMT)(L, my_helix);
    return 1;
} // end newHelix()


/**
 * The Lua constructor for a Bezier.
 *
 * Example construction in Lua:
 * ---------------------------------
 * p0 = Vector3:new{x=1}
 * p1 = Vector3:new{x=0.7071, y=0.7071}
 * p2 = Vector3:new{y=1}
 * -- For an arbitrary number of points in the table.
 * bez = Bezier:new{points={p0, p1, p2}}
 * ---------------------------------
 */
extern(C) int newBezier(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected Bezier:new{points={p0, p1, p2}}; ";
        errMsg ~= "maybe you tried Bezier.new{points={p0, p1, p2}}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to Bezier:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["points"])) {
        string errMsg = "Error in call to Bezier:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_getfield(L, 1, "points".toStringz());
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to Bezier:new{}. No points entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    if ( !lua_istable(L, -1) ) {
        string errMsg = "Error in call to Bezier:new{}.; " ~
            "A table containing Vector3 points is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect Vector3 objects at array positions within that table.
    Vector3[] B;
    int position = 1;
    while ( true ) {
        lua_rawgeti(L, -1, position);
        if ( lua_isnil(L, -1) ) { lua_pop(L, 1); break; }
        auto a = toVector3(L, -1);
        lua_pop(L, 1);
        B ~= a;
        ++position;
    }
    lua_pop(L, 1); // dispose of points table
    if ( B.length == 0 ) {
        string errMsg = "Error in call to Bezier:new{}. No valid Vector3 objects found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto bez = new Bezier(B);
    pathStore ~= pushObj!(Bezier, BezierMT)(L, bez);
    return 1;
} // end newBezier()

/**
 * A getter for number of Bezier control points.
 */
extern(C) int numberBezierCtrlPoints(lua_State* L)
{
    auto bezier = checkObj!(Bezier, BezierMT)(L, 1);
    lua_pushinteger(L, bezier.B.length);
    return 1;
}

/**
 * A getter/setter for Bezier control points.
 */
extern(C) int bezierCtrlPoint(lua_State* L)
{
    auto bezier = checkObj!(Bezier, BezierMT)(L, 1);
    int narg = lua_gettop(L);
    if (narg < 2) {
        string errMsg = "Error in call to bez:ctrlPt(): not enough arguments.\n";
        luaL_error(L, errMsg.toStringz);
    }
    int i = to!int(luaL_checkint(L, 2));
    if (i < 0 || i >= bezier.B.length) {
        string errMsg = "Error in call to bez:ctrPt(): index out of range.\n";
        errMsg ~= format("Index is: %d  No. of control points is: %d\n", i, bezier.B.length);
    }
    if (narg == 2) {
        // Treat as getter
        pushVector3(L, bezier.B[i]);
        return 1;
    }
    // Treat as setter
    bezier.B[i] = toVector3(L, 3);
    return 0;
}

/**
 * Access to the elevateDegree Bezier method.
 */
extern(C) int elevateDegree(lua_State* L)
{
    auto bezier = checkObj!(Bezier, BezierMT)(L, 1);
    int narg = lua_gettop(L);
    if (narg < 2) {
        string errMsg = "Error in call to bez:elevateDegree(): not enough arguments.\n";
        errMsg ~= "An integer for the desired degree elevation is required as single argument.\n";
        luaL_error(L, errMsg.toStringz);
    }
    int newDegree = to!int(luaL_checkint(L, 2));
    bezier.elevateDegree(newDegree);
    return 0;
}

/**
 * The Lua constructor for a NURBS.
 *
 * Example construction in Lua:
 * ---------------------------------
 * p0 = Vector3:new{x=0}
 * p1 = Vector3:new{x=0.7071, y=0.7071}
 * p2 = Vector3:new{y=1}
 * w = {1.0, 1.0, 1.0}
 * U = {0.0, 0.0, 0.5, 1.0, 1.0}
 * p = 2
 * nrb = NURBS:new{points={p0, p1, p2}, weights=w, knots=U, degree=p}
 * ---------------------------------
 */
extern(C) int newNURBS(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected NURBS:new{points=..., weights=..., knots=..., degree=...};\n ";
        errMsg ~= "maybe you tried NURBS.new{}}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to NURBS:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["points", "weights", "knots", "degree"])) {
        string errMsg = "Error in call to NURBS:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    // Get "points"
    lua_getfield(L, 1, "points".toStringz());
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to NURBS:new{}. No points entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    if ( !lua_istable(L, -1) ) {
        string errMsg = "Error in call to NURBS:new{}.; " ~
            "A table containing Vector3 points is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect Vector3 objects at array positions within that table.
    Vector3[] P;
    int position = 1;
    while ( true ) {
        lua_rawgeti(L, -1, position);
        if ( lua_isnil(L, -1) ) { lua_pop(L, 1); break; }
        auto a = toVector3(L, -1);
        lua_pop(L, 1);
        P ~= a;
        ++position;
    }
    lua_pop(L, 1); // dispose of points table
    if (P.length == 0) {
        string errMsg = "Error in call to NURBS:new{}. No valid Vector3 objects found.";
        luaL_error(L, errMsg.toStringz());
    }
    // Get "weights"
    double[] w;
    getArrayOfDoubles(L, 1, "weights", w);
    // Get "knots"
    double[] U;
    getArrayOfDoubles(L, 1, "knots", U);
    // Get "degree"
    int p = getInt(L, 1, "degree");
    // Prepare Pw array and build NURBS
    if (P.length != w.length) {
        string errMsg = "Error in call to NURBS:new{}. The length of the points table and weights table do not agree.\n";
        errMsg ~= format("#points= %d, #weights= %d", P.length, w.length);
    }
    double[4][] Pw;
    Pw.length = P.length;
    foreach (i; 0 .. P.length) Pw[i] = [P[i].x.re*w[i], P[i].y.re*w[i], P[i].z.re*w[i], w[i]];
    auto nrb = new NURBS(Pw, U, p);
    pathStore ~= pushObj!(NURBS, NURBSMT)(L, nrb);
    return 1;
} // end newBezier()


/**
 * The Lua constructor for a Polyline.
 *
 * Example construction in Lua:
 * ---------------------------------
 * p0 = Vector3:new{x=1}
 * p1 = Vector3:new{x=0.7071, y=0.7071}
 * p2 = Vector3:new{y=1}
 * -- A couple of paths to combine.
 * line1 = Line:new{p0=p0, p1=p1}
 * line2 = Line:new{p0=p1, p1=p2}
 * -- An arbitrary number of Path objects in the table.
 * poly = Polyline:new{segments={line1, line2}}
 * ---------------------------------
 */
extern(C) int newPolyline(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected Polyline:new{segments={path1, path2}}; ";
        errMsg ~= "maybe you tried Polyline.new{segments={path1, path2}}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to Polyline:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["segments"])) {
        string errMsg = "Error in call to Polyline:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_getfield(L, 1, "segments".toStringz());
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to Polyline:new{}. No segments entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    if ( !lua_istable(L, -1) ) {
        string errMsg = "Error in call to Polyline:new{}.; " ~
            "A table containing Vector3 points is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect Path objects at array positions within the segments table.
    Path[] segments;
    int position = 1;
    while ( true ) {
        lua_rawgeti(L, -1, position);
        if ( lua_isnil(L, -1) ) { lua_pop(L, 1); break; }
        auto seg = checkPath(L, -1);
        lua_pop(L, 1);
        if ( seg is null ) break;
        segments ~= seg;
        ++position;
    }
    lua_pop(L, 1); // dispose of segments table
    if ( segments.length == 0 ) {
        string errMsg = "Error in call to Polyline:new{}. No valid Path objects found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto poly = new Polyline(segments);
    pathStore ~= pushObj!(Polyline, PolylineMT)(L, poly);
    return 1;
} // end newPolyline()

/**
 * Access to the Polyline::toGmshString method.
 *
 * pTag, cTag, lTag, str = polyline:toGmshString(pointsTag, curvesTag, loopsTag, ?label?, ?len?)
 */
extern(C) int polyline_toGmshString(lua_State* L)
{
    auto polyline = checkObj!(Polyline, PolylineMT)(L, 1);
    int narg = lua_gettop(L);
    if (narg < 4) {
        string errMsg = "Error in call to Polyline:toGmshString(): not enough arguments.\n";
        errMsg ~= "Integer values for the start tags for points, curves and loops are required.\n";
        luaL_error(L, errMsg.toStringz);
    }
    int pointsTag = to!int(luaL_checkint(L, 2));
    int curvesTag = to!int(luaL_checkint(L, 3));
    int loopsTag = to!int(luaL_checkint(L, 4));
    string label = "";
    if (narg >= 5) { label = to!string(luaL_checkstring(L, 5)); }
    double len = 1.0e-2;
    if (narg >= 6) { len = to!double(luaL_checknumber(L, 6)); }
    string str = polyline.toGmshString(pointsTag, curvesTag, loopsTag, label, len);
    lua_settop(L, 0); // clear stack
    lua_pushnumber(L, pointsTag);
    lua_pushnumber(L, curvesTag);
    lua_pushnumber(L, loopsTag);
    lua_pushstring(L, str.toStringz);
    return 4;
}

/**
 * The Lua constructor for a Spline (Polyline).
 *
 * Example construction in Lua:
 * ---------------------------------
 * p0 = Vector3:new{ x=0, y=-1, z=0}
 * p1 = Vector3:new{x=-1,  y=0, z=0}
 * p2 = Vector3:new{ x=0,  y=1, z=0}
 * p3 = Vector3:new{ x=1,  y=0, z=0}
 * p4 = Vector3:new{ x=0, y=-1, z=0}
 * -- For an arbitrary number of points in the table.
 * spl = Spline:new{points={p0, p1, p2, p3, p4}}
 * ---------------------------------
 */
extern(C) int newSpline(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected Spline:new{points={p0, p1, p2, p3, p4}}; ";
        errMsg ~= "maybe you tried Spline.new{points={p0, p1, p2, p3, p4}}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to Spline:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["points"])) {
        string errMsg = "Error in call to Spline:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_getfield(L, 1, "points".toStringz());
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to Spline:new{}. No points entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    if ( !lua_istable(L, -1) ) {
        string errMsg = "Error in call to Spline:new{}.; " ~
            "A table containing Vector3 points is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect Vector3 objects at array positions.
    Vector3[] B;
    int position = 1;
    while ( true ) {
        lua_rawgeti(L, -1, position);
        if ( lua_isnil(L, -1) ) { lua_pop(L, 1); break; }
        auto a = toVector3(L, -1);
        lua_pop(L, 1);
        B ~= a;
        ++position;
    }
    lua_pop(L, 1); // dispose of points table
    if ( B.length == 0 ) {
        string errMsg = "Error in call to Spline:new{}. No valid Vector3 objects found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto spline = new Polyline(B);
    pathStore ~= pushObj!(Polyline, PolylineMT)(L, spline);
    return 1;
} // end newSpline()


/**
 * The Lua constructor for a Spline2 (Polyline).
 *
 * Example construction in Lua:
 * ---------------------------------
 * spl = Spline2:new{filename="something.dat"}
 * -- Expecting 3 numbers per line, space-separated.
 * ---------------------------------
 */
extern(C) int newSpline2(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected Spline2:new{filename=\"something.dat\"}; ";
        errMsg ~= "maybe you tried Spline2.new{filename=\"something.dat\"}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to Spline2:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["filename"])) {
        string errMsg = "Error in call to Spline2:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_getfield(L, 1, "filename".toStringz());
    if ( !lua_isstring(L, -1) ) {
        string errMsg = "Error in call to Spline2:new{}.; " ~
            "A string containing the file name is expected, but no string was found.";
        luaL_error(L, errMsg.toStringz);
    }
    auto fileName = to!string(lua_tostring(L, -1));
    lua_pop(L, 1); // dispose of filename string
    auto spline = new Polyline(fileName);
    pathStore ~= pushObj!(Polyline, PolylineMT)(L, spline);
    return 1;
} // end newSpline2()


/**
 * The Lua constructor for XSpline (CubicSpline function y(x)).
 *
 * Example construction in Lua:
 * ---------------------------------
 * -- For an arbitrary number of points.
 * spl = XSpline:new{xs={1, 2, 3, 4}, ys={1.5, 2.5, 2.5, 1.5}}
 * ---------------------------------
 */
extern(C) int newXSpline(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected XSpline:new{xs={...},ys={...}}; ";
        errMsg ~= "maybe you tried XSpline.new{xs={...},ys={...}}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to XSpline:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["xs", "ys"])) {
        string errMsg = "Error in call to XSpline:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_getfield(L, 1, "xs".toStringz());
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to XSpline:new{}. No xs entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    if ( !lua_istable(L, -1) ) {
        string errMsg = "Error in call to XSpline:new{}.; " ~
            "A table containing numbers for xs is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect numbers at array positions.
    double[] xs;
    int position = 1;
    while ( true ) {
        lua_rawgeti(L, -1, position);
        if ( lua_isnil(L, -1) ) { lua_pop(L, 1); break; }
        double a = to!double(lua_tonumber(L, -1));
        lua_pop(L, 1);
        xs ~= a;
        ++position;
    }
    lua_pop(L, 1); // dispose of xs table
    if ( xs.length == 0 ) {
        string errMsg = "Error in call to XSpline:new{}. No valid numbers found for xs.";
        luaL_error(L, errMsg.toStringz());
    }
    //
    lua_getfield(L, 1, "ys".toStringz());
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to XSpline:new{}. No ys entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    if ( !lua_istable(L, -1) ) {
        string errMsg = "Error in call to XSpline:new{}.; " ~
            "A table containing numbers is expected for ys, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    double[] ys;
    position = 1;
    while ( true ) {
        lua_rawgeti(L, -1, position);
        if ( lua_isnil(L, -1) ) { lua_pop(L, 1); break; }
        double a = to!double(lua_tonumber(L, -1));
        lua_pop(L, 1);
        ys ~= a;
        ++position;
    }
    lua_pop(L, 1); // dispose of ys table
    if ( ys.length == 0 ) {
        string errMsg = "Error in call to XSpline:new{}. No valid numbers found for ys.";
        luaL_error(L, errMsg.toStringz());
    }
    auto spline = new XSpline(xs, ys);
    pathStore ~= pushObj!(XSpline, XSplineMT)(L, spline);
    return 1;
} // end newXSpline()


/**
 * The Lua constructor for XSpline2 (CubicSpline function y(x)).
 *
 * Example construction in Lua:
 * ---------------------------------
 * spl = XSpline2:new{filename="something.dat"}
 * -- Expecting 2 numbers per line, space-separated.
 * ---------------------------------
 */
extern(C) int newXSpline2(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected XSpline2:new{filename=\"something.dat\"}; ";
        errMsg ~= "maybe you tried XSpline2.new{filename=\"something.dat\"}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to XSpline2:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["filename"])) {
        string errMsg = "Error in call to XSpline2:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_getfield(L, 1, "filename".toStringz());
    if ( !lua_isstring(L, -1) ) {
        string errMsg = "Error in call to XSpline2:new{}.; " ~
            "A string containing the file name is expected, but no string was found.";
        luaL_error(L, errMsg.toStringz);
    }
    auto fileName = to!string(lua_tostring(L, -1));
    lua_pop(L, 1); // dispose of filename string
    auto spline = new XSpline(fileName);
    pathStore ~= pushObj!(XSpline, XSplineMT)(L, spline);
    return 1;
} // end newXSpline2()


/**
 * The Lua constructor for XSplineLsq (CubicSplineLsq function y(x)).
 *
 * Example construction in Lua:
 * ---------------------------------
 * -- For an arbitrary number of points.
 * spl = XSplineLsq:new{xd={...}, yd={...}, wd={...}, xs={...}, nseg}
 * ---------------------------------
 */
extern(C) int newXSplineLsq(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected XSplineLsq:new{...}; ";
        errMsg ~= "maybe you tried XSplineLsq.new{...}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to XSplineLsq:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["xd", "yd", "wd", "xs", "nseg"])) {
        string errMsg = "Error in call to XSplineLsq:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_getfield(L, 1, "xd".toStringz());
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to XSplineLsq:new{}. No xd entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    if ( !lua_istable(L, -1) ) {
        string errMsg = "Error in call to XSplineLsq:new{}.; " ~
            "A table containing numbers for xd is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect numbers at array positions.
    double[] xd;
    int position = 1;
    while ( true ) {
        lua_rawgeti(L, -1, position);
        if ( lua_isnil(L, -1) ) { lua_pop(L, 1); break; }
        double a = to!double(lua_tonumber(L, -1));
        lua_pop(L, 1);
        xd ~= a;
        ++position;
    }
    lua_pop(L, 1); // dispose of xd table
    if ( xd.length == 0 ) {
        string errMsg = "Error in call to XSplineLsq:new{}. No valid numbers found for xd.";
        luaL_error(L, errMsg.toStringz());
    }
    //
    lua_getfield(L, 1, "yd".toStringz());
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to XSplineLsq:new{}. No yd entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    if ( !lua_istable(L, -1) ) {
        string errMsg = "Error in call to XSplineLsq:new{}.; " ~
            "A table containing numbers is expected for yd, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    double[] yd;
    position = 1;
    while ( true ) {
        lua_rawgeti(L, -1, position);
        if ( lua_isnil(L, -1) ) { lua_pop(L, 1); break; }
        double a = to!double(lua_tonumber(L, -1));
        lua_pop(L, 1);
        yd ~= a;
        ++position;
    }
    lua_pop(L, 1); // dispose of yd table
    if ( yd.length == 0 ) {
        string errMsg = "Error in call to XSpline:new{}. No valid numbers found for yd.";
        luaL_error(L, errMsg.toStringz());
    }
    //
    double[] wd;
    lua_getfield(L, 1, "wd".toStringz());
    if ( lua_istable(L, -1) ) {
        position = 1;
        while ( true ) {
            lua_rawgeti(L, -1, position);
            if ( lua_isnil(L, -1) ) { lua_pop(L, 1); break; }
            double a = to!double(lua_tonumber(L, -1));
            lua_pop(L, 1);
            wd ~= a;
            ++position;
        }
    }
    lua_pop(L, 1); // dispose of wd item (maybe nil)
    // It is ok to have no elements in wd.
    //
    double[] xs;
    lua_getfield(L, 1, "xs".toStringz());
    if ( lua_istable(L, -1) ) {
        position = 1;
        while ( true ) {
            lua_rawgeti(L, -1, position);
            if ( lua_isnil(L, -1) ) { lua_pop(L, 1); break; }
            double a = to!double(lua_tonumber(L, -1));
            lua_pop(L, 1);
            xs ~= a;
            ++position;
        }
    }
    lua_pop(L, 1); // dispose of xs item (maybe nil)
    // It is ok to have no elements in xs.
    //
    lua_getfield(L, 1, "nseg".toStringz());
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to XSplineLsq:new{}. No nseg entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    int nseg = to!int(lua_tointeger(L, -1));
    lua_pop(L, 1); // dispose of nseg item
    //
    auto spline = new XSplineLsq(xd, yd, wd, xs, nseg);
    pathStore ~= pushObj!(XSplineLsq, XSplineLsqMT)(L, spline);
    return 1;
} // end newXSplineLsq()


/**
 * The Lua constructor for XSplineLsq2 (CubicSplineLsq function y(x)).
 *
 * Example construction in Lua:
 * ---------------------------------
 * spl = XSplineLsq2:new{filename="something.dat", xs={...}, nseg=3}
 * -- Expecting 3 numbers per line (x, y, weight), space-separated.
 * ---------------------------------
 */
extern(C) int newXSplineLsq2(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected XSpline2:new{filename=\"something.dat\"}; ";
        errMsg ~= "maybe you tried XSpline2.new{filename=\"something.dat\"}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to XSpline2:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["filename", "xs", "nseg"])) {
        string errMsg = "Error in call to XSpline2:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_getfield(L, 1, "filename".toStringz());
    if ( !lua_isstring(L, -1) ) {
        string errMsg = "Error in call to XSpline2:new{}.; " ~
            "A string containing the file name is expected, but no string was found.";
        luaL_error(L, errMsg.toStringz);
    }
    auto fileName = to!string(lua_tostring(L, -1));
    lua_pop(L, 1); // dispose of filename string
    //
    double[] xs;
    lua_getfield(L, 1, "xs".toStringz());
    if ( lua_istable(L, -1) ) {
        int position = 1;
        while ( true ) {
            lua_rawgeti(L, -1, position);
            if ( lua_isnil(L, -1) ) { lua_pop(L, 1); break; }
            double a = to!double(lua_tonumber(L, -1));
            lua_pop(L, 1);
            xs ~= a;
            ++position;
        }
    }
    lua_pop(L, 1); // dispose of xs item (maybe nil)
    // It is ok to have no elements in xs.
    //
    lua_getfield(L, 1, "nseg".toStringz());
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to XSplineLsq2:new{}. No nseg entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    int nseg = to!int(lua_tointeger(L, -1));
    lua_pop(L, 1); // dispose of nseg item
    //
    auto spline = new XSplineLsq(fileName, xs, nseg);
    pathStore ~= pushObj!(XSplineLsq, XSplineLsqMT)(L, spline);
    return 1;
} // end newXSplineLsq2()


/**
 * The Lua constructor for a SVGPath (subclass of Polyline).
 *
 * Example construction in Lua:
 * ---------------------------------
 * svgpth = SVGPath:new{path="M3.0,3.0;L4.0,3.0;v1.0;h-1.0;Z",
 *                      xtrans=0, ytrans=0, xscale=1, yscale=1,
 *                      tolerance=1.0e-10}
 * -- Only the path string is required.
 * -- In that string, we are expecting 3 numbers per line, space-separated.
 * ---------------------------------
 */
extern(C) int newSVGPath(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected SVGPath:new{path=\"some-path-commands\"}; ";
        errMsg ~= "maybe you tried SVGPath.new{path=\"some-path-commands\"}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to SVGPath:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["path", "xtrans", "ytrans", "xscale", "yscale", "tolerance"])) {
        string errMsg = "Error in call to SVGPath:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_getfield(L, 1, "path".toStringz());
    if ( !lua_isstring(L, -1) ) {
        string errMsg = "Error in call to SVGPath:new{}.; " ~
            "A string containing the path specification is expected, but no string was found.";
        luaL_error(L, errMsg.toStringz);
    }
    auto pathTxt = to!string(lua_tostring(L, -1));
    lua_pop(L, 1); // dispose of path string
    //
    double xtrans = 0.0;
    lua_getfield(L, 1, "xtrans".toStringz());
    if ( lua_isnumber(L, -1) ) { xtrans = to!double(lua_tonumber(L, -1)); }
    lua_pop(L, 1); // dispose of xtrans item
    double ytrans = 0.0;
    lua_getfield(L, 1, "ytrans".toStringz());
    if ( lua_isnumber(L, -1) ) { ytrans = to!double(lua_tonumber(L, -1)); }
    lua_pop(L, 1); // dispose of ytrans item
    double xscale = 1.0;
    lua_getfield(L, 1, "xscale".toStringz());
    if ( lua_isnumber(L, -1) ) { xscale = to!double(lua_tonumber(L, -1)); }
    lua_pop(L, 1); // dispose of xscale item
    double yscale = 1.0;
    lua_getfield(L, 1, "yscale".toStringz());
    if ( lua_isnumber(L, -1) ) { yscale = to!double(lua_tonumber(L, -1)); }
    lua_pop(L, 1); // dispose of yscale item
    double tolerance = 1.0e-10;
    lua_getfield(L, 1, "tolerance".toStringz());
    if ( lua_isnumber(L, -1) ) { tolerance = to!double(lua_tonumber(L, -1)); }
    lua_pop(L, 1); // dispose of tolerance item
    //
    auto svgpath = new SVGPath(pathTxt, xtrans, ytrans, xscale, yscale, tolerance);
    pathStore ~= pushObj!(SVGPath, SVGPathMT)(L, svgpath);
    return 1;
} // end newSVGPath()

/**
 * Access to the SVGPath::toGmshString method.
 *
 * pTag, cTag, lTag, str = svgpath:toGmshString(pointsTag, curvesTag, loopsTag, ?label?, ?len?)
 */
extern(C) int svgpath_toGmshString(lua_State* L)
{
    auto svgpath = checkObj!(SVGPath, SVGPathMT)(L, 1);
    int narg = lua_gettop(L);
    if (narg < 4) {
        string errMsg = "Error in call to SVGPath:toGmshString(): not enough arguments.\n";
        errMsg ~= "Integer values for the start tags for points, curves and loops are required.\n";
        luaL_error(L, errMsg.toStringz);
    }
    int pointsTag = to!int(luaL_checkint(L, 2));
    int curvesTag = to!int(luaL_checkint(L, 3));
    int loopsTag = to!int(luaL_checkint(L, 4));
    string label = "";
    if (narg >= 5) { label = to!string(luaL_checkstring(L, 5)); }
    double len = 1.0e-2;
    if (narg >= 6) { len = to!double(luaL_checknumber(L, 6)); }
    string str = svgpath.toGmshString(pointsTag, curvesTag, loopsTag, label, len);
    lua_settop(L, 0); // clear stack
    lua_pushnumber(L, pointsTag);
    lua_pushnumber(L, curvesTag);
    lua_pushnumber(L, loopsTag);
    lua_pushstring(L, str.toStringz);
    return 4;
}


/**
 * LuaFnPath class and it's Lua constructor.
 *
 * This is hangs onto a Lua call-back function that is invoked from the D domain.
 *
 * Example:
 * function myLuaFunction(t)
 *    -- Straight line from 0,0,0 to 1.0,2.0,3.0
 *    return {x=t, y=2*t, z=3*t}
 * end
 * myPath = LuaFnPath:new{luaFnName="myLuaFunction"}
 */

class LuaFnPath : Path {
public:
    lua_State* L; // a pointer to the Lua interpreter's state.
    // Even though some of the class methods claim that they don't change
    // the object state, we have to get the Lua interpreter to evaluate
    // things and that diddles with the Lua interpreter's internal state.
    // So the const on the lua_State pointer is more a statement that
    // "I'm not going to switch interpreters on you."
    // Hence the ugly but (hopefully safe) casts where ever we get
    // the Lua interpreter to do something.
    // This is the best I can do for the moment.  PJ, 2014-04-22
    string luaFnName;
    this(const lua_State* L, string luaFnName)
    {
        this.L = cast(lua_State*)L;
        this.luaFnName = luaFnName;
    }
    this(ref const(LuaFnPath) other)
    {
        L = cast(lua_State*)other.L;
        luaFnName = other.luaFnName;
    }
    override LuaFnPath dup() const
    {
        return new LuaFnPath(L, luaFnName);
    }
    override Vector3 opCall(double t) const
    {
        // Call back to the Lua function.
        lua_getglobal(cast(lua_State*)L, luaFnName.toStringz);
        lua_pushnumber(cast(lua_State*)L, t);
        if ( lua_pcall(cast(lua_State*)L, 1, 1, 0) != 0 ) {
            string errMsg = "Error in call to " ~ luaFnName ~
                " from LuaFnPath:opCall(): " ~
                to!string(lua_tostring(cast(lua_State*)L, -1));
            luaL_error(cast(lua_State*)L, errMsg.toStringz);
        }
        // We are expecting a table to be returned, containing three numbers.
        if ( !lua_istable(cast(lua_State*)L, -1) ) {
            string errMsg = "Error in call to LuaFnPath:opCall().; " ~
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
        return "LuaFnPath(luaFnName=\"" ~ luaFnName ~ "\")";
    }
    override string classString() const
    {
        return "LuaFnPath";
    }
} // end class LuaFnPath

extern(C) int newLuaFnPath(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected LuaFnPath:new{luaFnName=\"myLuaFunction\"}; ";
        errMsg ~= "maybe you tried LuaFnPath.new{luaFnName=\"myLuaFunction\"}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to LuaFnPath:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["luaFnName"])) {
        string errMsg = "Error in call to LuaFnPath:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect function name in table.
    string fnName = "";
    lua_getfield(L, 1, "luaFnName".toStringz());
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to LuaFnPath:new{}. No luaFnName entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    if ( lua_isstring(L, -1) ) {
        fnName ~= to!string(lua_tostring(L, -1));
    }
    lua_pop(L, 1);
    if ( fnName == "" ) {
        string errMsg = "Error in call to LuaFnPath:new{}. No function name found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto lfp = new LuaFnPath(L, fnName);
    pathStore ~= pushObj!(LuaFnPath, LuaFnPathMT)(L, lfp);
    return 1;
} // end newLuaFnPath()


/**
 * The Lua constructor for an ArcLengthParameterizedPath.
 *
 * Example construction in Lua:
 * ---------------------------------
 * p0 = Vector3:new{x=1}
 * p1 = Vector3:new{x=0.7071, y=0.7071}
 * p2 = Vector3:new{x=0, y=1}
 * original_path = Bezier:new{points={p0, p1, p2}}
 * alp_path = ArcLengthParameterizedPath:new{underlying_path=original_path}
 * ---------------------------------
 */
extern(C) int newArcLengthParameterizedPath(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected ArcLengthParameterizedPath:new{}; ";
        errMsg ~= "maybe you tried ArcLengthParameterizedPath.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to ArcLengthParameterizedPath:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if ( !checkAllowedNames(L, 1, ["underlying_path"]) ) {
        string errMsg = "Error in call to ArcLengthParameterizedPath:new{}. " ~
            "Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect the underlying_path object in the table.
    lua_getfield(L, 1, "underlying_path".toStringz());
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to ArcLengthParameterizedPath:new{}." ~
            " No underlying_path entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto underlying_path = checkPath(L, -1);
    lua_pop(L, 1);
    if ( underlying_path is null ) {
        string errMsg = "Error in call to ArcLengthParameterizedPath:new{};" ~
            " Not a valid Path object.";
        luaL_error(L, errMsg.toStringz());
    }
    auto alp_path = new ArcLengthParameterizedPath(underlying_path);
    pathStore ~= pushObj!(ArcLengthParameterizedPath,
                          ArcLengthParameterizedPathMT)(L, alp_path);
    return 1;
} // end newArcLengthParameterizedPath()


/**
 * The Lua constructor for an SubRangedPath.
 *
 * Example construction in Lua:
 * ---------------------------------
 * p0 = Vector3:new{x=1}
 * p1 = Vector3:new{x=0.7071, y=0.7071}
 * p2 = Vector3:new{y=1}
 * original_path = Bezier:new{points={p0, p1, p2}}
 * sr_path = SubRangedPath:new{underlying_path=original_path, t0=0.1, t1=0.9}
 * ---------------------------------
 */
extern(C) int newSubRangedPath(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected SubRangedPath:new{}; ";
        errMsg ~= "maybe you tried SubRangedPath.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to SubRangedPath:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["underlying_path", "t0", "t1"])) {
        string errMsg = "Error in call to SubRangedPath:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect the underlying_path object in the table.
    lua_getfield(L, 1, "underlying_path".toStringz());
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to SubRangedPath:new{}." ~
            " No underlying_path entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto underlying_path = checkPath(L, -1);
    lua_pop(L, 1);
    if ( underlying_path is null ) {
        string errMsg = "Error in call to SubRangedPath:new{};" ~
            " Not a valid Path object.";
        luaL_error(L, errMsg.toStringz());
    }
    string errMsgTmplt = "Error in call to SubRangedPath:new{}. " ~
        "A valid value for '%s' was not found in list of arguments. " ~
        "The value, if present, should be a number.";
    double t0 = getNumberFromTable(L, 1, "t0", false, 0.0, true, format(errMsgTmplt, "t0"));
    double t1 = getNumberFromTable(L, 1, "t1", false, 1.0, true, format(errMsgTmplt, "t1"));
    auto alp_path = new SubRangedPath(underlying_path, t0, t1);
    pathStore ~= pushObj!(SubRangedPath, SubRangedPathMT)(L, alp_path);
    return 1;
} // end newSubRangedPath()

extern(C) int t0Path(T, string MTname)(lua_State* L)
{
    // Not much error checking here because we assume
    // users are knowing what they are doing if
    // they are messing with the getter/setter functions.
    int narg = lua_gettop(L);
    auto path = checkObj!(T, MTname)(L, 1);
    if ( narg == 1 ) { // This is a getter
        lua_pushnumber(L, path.t0);
        return 1;
    }
    // else: treat as a setter
    path.t0 = luaL_checknumber(L, 2);
    return 0;
}

extern(C) int t1Path(T, string MTname)(lua_State* L)
{
    // Not much error checking here because we assume
    // users are knowing what they are doing if
    // they are messing with the getter/setter functions.
    int narg = lua_gettop(L);
    auto path = checkObj!(T, MTname)(L, 1);
    if ( narg == 1 ) { // This is a getter
        lua_pushnumber(L, path.t1);
        return 1;
    }
    // else: treat as a setter
    path.t1 = luaL_checknumber(L, 2);
    return 0;
}


/**
 * The Lua constructor for an ReversedPath.
 *
 * Example construction in Lua:
 * ---------------------------------
 * p0 = Vector3:new{x=1}
 * p1 = Vector3:new{x=0.7071, y=0.7071}
 * p2 = Vector3:new{y=1}
 * original_path = Bezier:new{points={p0, p1, p2}}
 * r_path = ReversedPath:new{underlying_path=original_path}
 * ---------------------------------
 */
extern(C) int newReversedPath(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected ReversedPath:new{}; ";
        errMsg ~= "maybe you tried ReversedPath.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to ReversedPath:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["underlying_path"])) {
        string errMsg = "Error in call to ReversedPath:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect the underlying_path object in the table.
    lua_getfield(L, 1, "underlying_path".toStringz());
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to ReversedPath:new{}." ~
            " No underlying_path entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto underlying_path = checkPath(L, -1);
    lua_pop(L, 1);
    if ( underlying_path is null ) {
        string errMsg = "Error in call to ReversedPath:new{};" ~
            " Not a valid Path object.";
        luaL_error(L, errMsg.toStringz());
    }
    auto alp_path = new ReversedPath(underlying_path);
    pathStore ~= pushObj!(ReversedPath, ReversedPathMT)(L, alp_path);
    return 1;
} // end newReversedPath()


/**
 * The Lua constructor for an TranslatedPath.
 *
 * Example construction in Lua:
 * ---------------------------------
 * p0 = Vector3:new{x=1}
 * p1 = Vector3:new{x=0.7071, y=0.7071}
 * p2 = Vector3:new{y=1}
 * opath = Bezier:new{points={p0, p1, p2}}
 * tr_path = TranslatedPath:new{original_path=opath, shift=Vector3:new{0.5, 0.5}}
 * ---------------------------------
 */
extern(C) int newTranslatedPath(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected TranslatedPath:new{}; ";
        errMsg ~= "maybe you tried TranslatedPath.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to TranslatedPath:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["original_path", "shift"])) {
        string errMsg = "Error in call to TranslatedPath:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect the original_path object in the table.
    lua_getfield(L, 1, "original_path".toStringz());
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to TranslatedPath:new{}." ~
            " No original_path entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto original_path = checkPath(L, -1);
    lua_pop(L, 1);
    if ( original_path is null ) {
        string errMsg = "Error in call to TranslatedPath:new{};" ~
            " Not a valid Path object.";
        luaL_error(L, errMsg.toStringz());
    }
    // Expect Vector3 for key "shift".
    lua_getfield(L, 1, "shift".toStringz());
    auto shift = toVector3(L, -1);
    lua_pop(L, 1);
    auto tr_path = new TranslatedPath(original_path, shift);
    pathStore ~= pushObj!(TranslatedPath, TranslatedPathMT)(L, tr_path);
    return 1;
} // end newTranslatedPath()

/**
 * The Lua constructor for an MirrorImagePath.
 *
 * Example construction in Lua:
 * ---------------------------------
 * p0 = Vector3:new{x=1}
 * p1 = Vector3:new{x=0.7071, y=0.7071}
 * p2 = Vector3:new{y=1}
 * opath = Bezier:new{points={p0, p1, p2}}
 * mi_path = MirrorImagePath:new{original_path=opath,
 *                               point=Vector3:new{1.0, 0.0},
 *                               normal=Vector3:new{1.0, 0.0}}
 * ---------------------------------
 */
extern(C) int newMirrorImagePath(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected MirrorImagePath:new{}; ";
        errMsg ~= "maybe you tried MirrorImagePath.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to MirrorImagePath:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["original_path", "point", "normal"])) {
        string errMsg = "Error in call to MirrorImagePath:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect the original_path object in the table.
    lua_getfield(L, 1, "original_path".toStringz());
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to MirrorImagePath:new{}." ~
            " No original_path entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto original_path = checkPath(L, -1);
    lua_pop(L, 1);
    if ( original_path is null ) {
        string errMsg = "Error in call to MirrorImagePath:new{}. Not a valid Path object.";
        luaL_error(L, errMsg.toStringz());
    }
    // Expect Vector3 for key "point".
    lua_getfield(L, 1, "point".toStringz());
    auto point = toVector3(L, -1);
    lua_pop(L, 1);
    // Expect Vector3 for key "normal".
    lua_getfield(L, 1, "normal".toStringz());
    auto normal = toVector3(L, -1);
    lua_pop(L, 1);
    auto mi_path = new MirrorImagePath(original_path, point, normal);
    pathStore ~= pushObj!(MirrorImagePath, MirrorImagePathMT)(L, mi_path);
    return 1;
} // end newMirrorImagePath()

/**
 * The Lua constructor for an RotatedAboutZAxisPath.
 *
 * Example construction in Lua:
 * ---------------------------------
 * p0 = Vector3:new{x=1}
 * p1 = Vector3:new{x=0.7071, y=0.7071}
 * p2 = Vector3:new{y=1}
 * opath = Bezier:new{points={p0, p1, p2}}
 * raza_path = RotatedAboutZAxisPath:new{original_path=opath, angle=math.pi/4}
 * ---------------------------------
 */
extern(C) int newRotatedAboutZAxisPath(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected RotatedAboutZAxisPath:new{}; ";
        errMsg ~= "maybe you tried RotatedAboutZAxisPath.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to RotatedAboutZAxisPath:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["original_path", "angle"])) {
        string errMsg = "Error in call to RotatedAboutZAxisPath:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect the original_path object in the table.
    lua_getfield(L, 1, "original_path".toStringz());
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to RotatedAboutZAxisPath:new{}." ~
            " No original_path entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto original_path = checkPath(L, -1);
    lua_pop(L, 1);
    if ( original_path is null ) {
        string errMsg = "Error in call to RotatedAboutZAxisPath:new{}. Not a valid Path object.";
        luaL_error(L, errMsg.toStringz());
    }
    string errMsgTmplt = "Error in call to RotatedAboutZAxisPath:new{}. " ~
        "A valid value for '%s' was not found in list of arguments. " ~
        "The value, if present, should be a number.";
    double angle = getNumberFromTable(L, 1, "angle", false, 0.0, true,
                                      format(errMsgTmplt, "angle"));
    auto raza_path = new RotatedAboutZAxisPath(original_path, angle);
    pathStore ~= pushObj!(RotatedAboutZAxisPath, RotatedAboutZAxisPathMT)(L, raza_path);
    return 1;
} // end newRotatedAboutZAxisPath()

//-------------------------------------------------------------------------------------

void registerPaths(lua_State* L)
{
    // Register the Line object
    luaL_newmetatable(L, LineMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newLine);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(Line, LineMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(Line, LineMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(Line, LineMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(Line, LineMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &pathIntersect2D!(Line, LineMT));
    lua_setfield(L, -2, "intersect2D");
    lua_pushcfunction(L, &pathClosestDistance!(Line, LineMT));
    lua_setfield(L, -2, "closestDistance");
    lua_pushcfunction(L, &length!(Line, LineMT));
    lua_setfield(L, -2, "length");

    lua_setglobal(L, LineMT.toStringz);

    // Register the Arc object
    luaL_newmetatable(L, ArcMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newArc);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(Arc, ArcMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(Arc, ArcMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(Arc, ArcMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(Arc, ArcMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &pathIntersect2D!(Arc, ArcMT));
    lua_setfield(L, -2, "intersect2D");
    lua_pushcfunction(L, &pathClosestDistance!(Arc, ArcMT));
    lua_setfield(L, -2, "closestDistance");
    lua_pushcfunction(L, &length!(Arc, ArcMT));
    lua_setfield(L, -2, "length");

    lua_setglobal(L, ArcMT.toStringz);

    // Register the Arc3 object
    luaL_newmetatable(L, Arc3MT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newArc3);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(Arc3, Arc3MT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(Arc3, Arc3MT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(Arc3, Arc3MT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(Arc3, Arc3MT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &pathIntersect2D!(Arc3, Arc3MT));
    lua_setfield(L, -2, "intersect2D");
    lua_pushcfunction(L, &pathClosestDistance!(Arc3, Arc3MT));
    lua_setfield(L, -2, "closestDistance");
    lua_pushcfunction(L, &length!(Arc3, Arc3MT));
    lua_setfield(L, -2, "length");

    lua_setglobal(L, Arc3MT.toStringz);

    // Register the Helix object
    luaL_newmetatable(L, HelixMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newHelix);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(Helix, HelixMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(Helix, HelixMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(Helix, HelixMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(Helix, HelixMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &pathClosestDistance!(Helix, HelixMT));
    lua_setfield(L, -2, "closestDistance");
    lua_pushcfunction(L, &length!(Helix, HelixMT));
    lua_setfield(L, -2, "length");

    lua_setglobal(L, HelixMT.toStringz);

    // Register the Bezier object
    luaL_newmetatable(L, BezierMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newBezier);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(Bezier, BezierMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(Bezier, BezierMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(Bezier, BezierMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(Bezier, BezierMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &pathIntersect2D!(Bezier, BezierMT));
    lua_setfield(L, -2, "intersect2D");
    lua_pushcfunction(L, &pathClosestDistance!(Bezier, BezierMT));
    lua_setfield(L, -2, "closestDistance");
    lua_pushcfunction(L, &length!(Bezier, BezierMT));
    lua_setfield(L, -2, "length");
    lua_pushcfunction(L, &numberBezierCtrlPoints);
    lua_setfield(L, -2, "numberCtrlPts");
    lua_pushcfunction(L, &bezierCtrlPoint);
    lua_setfield(L, -2, "ctrlPt");
    lua_pushcfunction(L, &elevateDegree);
    lua_setfield(L, -2, "elevateDegree");

    lua_setglobal(L, BezierMT.toStringz);

    // Register the NURBS object
    luaL_newmetatable(L, NURBSMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newNURBS);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(NURBS, NURBSMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(NURBS, NURBSMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(NURBS, NURBSMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(NURBS, NURBSMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &pathIntersect2D!(NURBS, NURBSMT));
    lua_setfield(L, -2, "intersect2D");
    lua_pushcfunction(L, &pathClosestDistance!(NURBS, NURBSMT));
    lua_setfield(L, -2, "closestDistance");
    lua_pushcfunction(L, &length!(NURBS, NURBSMT));
    lua_setfield(L, -2, "length");

    lua_setglobal(L, NURBSMT.toStringz);

    // Register the Polyline object
    luaL_newmetatable(L, PolylineMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newPolyline);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(Polyline, PolylineMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(Polyline, PolylineMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(Polyline, PolylineMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(Polyline, PolylineMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &pathIntersect2D!(Polyline, PolylineMT));
    lua_setfield(L, -2, "intersect2D");
    lua_pushcfunction(L, &pathClosestDistance!(Polyline, PolylineMT));
    lua_setfield(L, -2, "closestDistance");
    lua_pushcfunction(L, &length!(Polyline, PolylineMT));
    lua_setfield(L, -2, "length");
    lua_pushcfunction(L, &polyline_toGmshString);
    lua_setfield(L, -2, "toGmshString");

    lua_setglobal(L, PolylineMT.toStringz);

    // Register the Spline object which is actually a Polyline in Dlang.
    luaL_newmetatable(L, SplineMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newSpline);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(Polyline, SplineMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(Polyline, SplineMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(Polyline, SplineMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(Polyline, SplineMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &pathIntersect2D!(Polyline, SplineMT));
    lua_setfield(L, -2, "intersect2D");
    lua_pushcfunction(L, &pathClosestDistance!(Polyline, SplineMT));
    lua_setfield(L, -2, "closestDistance");
    lua_pushcfunction(L, &length!(Polyline, SplineMT));
    lua_setfield(L, -2, "length");

    lua_setglobal(L, SplineMT.toStringz);

    // Register the Spline2 object which is also a Polyline in Dlang.
    luaL_newmetatable(L, Spline2MT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newSpline2);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(Polyline, Spline2MT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(Polyline, Spline2MT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(Polyline, Spline2MT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(Polyline, Spline2MT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &pathIntersect2D!(Polyline, Spline2MT));
    lua_setfield(L, -2, "intersect2D");
    lua_pushcfunction(L, &pathClosestDistance!(Polyline, Spline2MT));
    lua_setfield(L, -2, "closestDistance");
    lua_pushcfunction(L, &length!(Polyline, Spline2MT));
    lua_setfield(L, -2, "length");

    lua_setglobal(L, Spline2MT.toStringz);

    // Register the XSpline object.
    luaL_newmetatable(L, XSplineMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newXSpline);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(XSpline, XSplineMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(XSpline, XSplineMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(XSpline, XSplineMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(XSpline, XSplineMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &pathIntersect2D!(XSpline, XSplineMT));
    lua_setfield(L, -2, "intersect2D");
    lua_pushcfunction(L, &pathClosestDistance!(XSpline, XSplineMT));
    lua_setfield(L, -2, "closestDistance");
    lua_pushcfunction(L, &length!(XSpline, XSplineMT));
    lua_setfield(L, -2, "length");

    lua_setglobal(L, XSplineMT.toStringz);

    // Register the XSpline2 object which is actually a XSpline object in D.
    luaL_newmetatable(L, XSpline2MT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newXSpline2);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(XSpline, XSpline2MT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(XSpline, XSpline2MT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(XSpline, XSpline2MT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(XSpline, XSpline2MT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &pathIntersect2D!(XSpline, XSpline2MT));
    lua_setfield(L, -2, "intersect2D");
    lua_pushcfunction(L, &pathClosestDistance!(XSpline, XSpline2MT));
    lua_setfield(L, -2, "closestDistance");
    lua_pushcfunction(L, &length!(XSpline, XSpline2MT));
    lua_setfield(L, -2, "length");

    lua_setglobal(L, XSpline2MT.toStringz);


    // Register the XSplineLsq object.
    luaL_newmetatable(L, XSplineLsqMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newXSplineLsq);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(XSplineLsq, XSplineLsqMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(XSplineLsq, XSplineLsqMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(XSplineLsq, XSplineLsqMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(XSplineLsq, XSplineLsqMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &pathIntersect2D!(XSplineLsq, XSplineLsqMT));
    lua_setfield(L, -2, "intersect2D");
    lua_pushcfunction(L, &pathClosestDistance!(XSplineLsq, XSplineLsqMT));
    lua_setfield(L, -2, "closestDistance");
    lua_pushcfunction(L, &length!(XSplineLsq, XSplineLsqMT));
    lua_setfield(L, -2, "length");

    lua_setglobal(L, XSplineLsqMT.toStringz);

    // Register the XSplineLsq2 object which is actually a XSplineLsq object in D.
    luaL_newmetatable(L, XSplineLsq2MT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newXSplineLsq2);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(XSplineLsq, XSplineLsq2MT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(XSplineLsq, XSplineLsq2MT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(XSplineLsq, XSplineLsq2MT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(XSplineLsq, XSplineLsq2MT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &pathIntersect2D!(XSplineLsq, XSplineLsq2MT));
    lua_setfield(L, -2, "intersect2D");
    lua_pushcfunction(L, &pathClosestDistance!(XSplineLsq, XSplineLsq2MT));
    lua_setfield(L, -2, "closestDistance");
    lua_pushcfunction(L, &length!(XSplineLsq, XSplineLsq2MT));
    lua_setfield(L, -2, "length");

    lua_setglobal(L, XSplineLsq2MT.toStringz);

    // Register the SVGPath object which is a subclass of Polyline in Dlang.
    luaL_newmetatable(L, SVGPathMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newSVGPath);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(SVGPath, SVGPathMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(SVGPath, SVGPathMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(SVGPath, SVGPathMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(SVGPath, SVGPathMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &pathIntersect2D!(SVGPath, SVGPathMT));
    lua_setfield(L, -2, "intersect2D");
    lua_pushcfunction(L, &pathClosestDistance!(SVGPath, SVGPathMT));
    lua_setfield(L, -2, "closestDistance");
    lua_pushcfunction(L, &length!(SVGPath, SVGPathMT));
    lua_setfield(L, -2, "length");
    lua_pushcfunction(L, &svgpath_toGmshString);
    lua_setfield(L, -2, "toGmshString");

    lua_setglobal(L, SVGPathMT.toStringz);

    // Register the LuaFnPath object
    luaL_newmetatable(L, LuaFnPathMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newLuaFnPath);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(LuaFnPath, LuaFnPathMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(LuaFnPath, LuaFnPathMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(LuaFnPath, LuaFnPathMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(LuaFnPath, LuaFnPathMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &pathIntersect2D!(LuaFnPath, LuaFnPathMT));
    lua_setfield(L, -2, "intersect2D");
    lua_pushcfunction(L, &pathClosestDistance!(LuaFnPath, LuaFnPathMT));
    lua_setfield(L, -2, "closestDistance");
    lua_pushcfunction(L, &length!(LuaFnPath, LuaFnPathMT));
    lua_setfield(L, -2, "length");

    lua_setglobal(L, LuaFnPathMT.toStringz);

    // Register the ArcLengthParameterized object
    luaL_newmetatable(L, ArcLengthParameterizedPathMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newArcLengthParameterizedPath);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(ArcLengthParameterizedPath,
                                      ArcLengthParameterizedPathMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(ArcLengthParameterizedPath,
                                      ArcLengthParameterizedPathMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(ArcLengthParameterizedPath,
                                       ArcLengthParameterizedPathMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(ArcLengthParameterizedPath,
                                    ArcLengthParameterizedPathMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &pathIntersect2D!(ArcLengthParameterizedPath,
                                           ArcLengthParameterizedPathMT));
    lua_setfield(L, -2, "intersect2D");
    lua_pushcfunction(L, &pathClosestDistance!(ArcLengthParameterizedPath,
                                               ArcLengthParameterizedPathMT));
    lua_setfield(L, -2, "closestDistance");
    lua_pushcfunction(L, &length!(ArcLengthParameterizedPath,
                                  ArcLengthParameterizedPathMT));
    lua_setfield(L, -2, "length");

    lua_setglobal(L, ArcLengthParameterizedPathMT.toStringz);

    // Register the SubRangedPath object
    luaL_newmetatable(L, SubRangedPathMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newSubRangedPath);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(SubRangedPath, SubRangedPathMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(SubRangedPath, SubRangedPathMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(SubRangedPath, SubRangedPathMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(SubRangedPath, SubRangedPathMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &t0Path!(SubRangedPath, SubRangedPathMT));
    lua_setfield(L, -2, "t0");
    lua_pushcfunction(L, &t1Path!(SubRangedPath, SubRangedPathMT));
    lua_setfield(L, -2, "t1");
    lua_pushcfunction(L, &pathIntersect2D!(SubRangedPath, SubRangedPathMT));
    lua_setfield(L, -2, "intersect2D");
    lua_pushcfunction(L, &pathClosestDistance!(SubRangedPath, SubRangedPathMT));
    lua_setfield(L, -2, "closestDistance");

    lua_setglobal(L, SubRangedPathMT.toStringz);

    // Register the ReversedPath object
    luaL_newmetatable(L, ReversedPathMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newReversedPath);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(ReversedPath, ReversedPathMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(ReversedPath, ReversedPathMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(ReversedPath, ReversedPathMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(ReversedPath, ReversedPathMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &t0Path!(ReversedPath, ReversedPathMT));
    lua_setfield(L, -2, "t0");
    lua_pushcfunction(L, &t1Path!(ReversedPath, ReversedPathMT));
    lua_setfield(L, -2, "t1");
    lua_pushcfunction(L, &pathIntersect2D!(ReversedPath, ReversedPathMT));
    lua_setfield(L, -2, "intersect2D");
    lua_pushcfunction(L, &pathClosestDistance!(ReversedPath, ReversedPathMT));
    lua_setfield(L, -2, "closestDistance");
    lua_pushcfunction(L, &length!(ReversedPath, ReversedPathMT));
    lua_setfield(L, -2, "length");

    lua_setglobal(L, ReversedPathMT.toStringz);

    // Register the TranslatedPath object
    luaL_newmetatable(L, TranslatedPathMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newTranslatedPath);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(TranslatedPath, TranslatedPathMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(TranslatedPath, TranslatedPathMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(TranslatedPath, TranslatedPathMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(TranslatedPath, TranslatedPathMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &pathIntersect2D!(TranslatedPath, TranslatedPathMT));
    lua_setfield(L, -2, "intersect2D");
    lua_pushcfunction(L, &pathClosestDistance!(TranslatedPath, TranslatedPathMT));
    lua_setfield(L, -2, "closestDistance");
    lua_pushcfunction(L, &length!(TranslatedPath, TranslatedPathMT));
    lua_setfield(L, -2, "length");

    lua_setglobal(L, TranslatedPathMT.toStringz);

    // Register the MirrorImagePath object
    luaL_newmetatable(L, MirrorImagePathMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newMirrorImagePath);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(MirrorImagePath, MirrorImagePathMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(MirrorImagePath, MirrorImagePathMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(MirrorImagePath, MirrorImagePathMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(MirrorImagePath, MirrorImagePathMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &pathIntersect2D!(MirrorImagePath, MirrorImagePathMT));
    lua_setfield(L, -2, "intersect2D");
    lua_pushcfunction(L, &pathClosestDistance!(MirrorImagePath, MirrorImagePathMT));
    lua_setfield(L, -2, "closestDistance");
    lua_pushcfunction(L, &length!(MirrorImagePath, MirrorImagePathMT));
    lua_setfield(L, -2, "length");

    lua_setglobal(L, MirrorImagePathMT.toStringz);

    // Register the RotatedAboutZAxisPath object
    luaL_newmetatable(L, RotatedAboutZAxisPathMT.toStringz);

    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    lua_pushcfunction(L, &newRotatedAboutZAxisPath);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallPath!(RotatedAboutZAxisPath, RotatedAboutZAxisPathMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallPath!(RotatedAboutZAxisPath, RotatedAboutZAxisPathMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(RotatedAboutZAxisPath, RotatedAboutZAxisPathMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyPath!(RotatedAboutZAxisPath, RotatedAboutZAxisPathMT));
    lua_setfield(L, -2, "copy");
    lua_pushcfunction(L, &length!(RotatedAboutZAxisPath, RotatedAboutZAxisPathMT));
    lua_setfield(L, -2, "length");

    lua_setglobal(L, RotatedAboutZAxisPathMT.toStringz);
} // end registerPaths()
