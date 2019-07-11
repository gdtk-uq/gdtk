/**
 * A Lua interface to the gpath_utils module.
 *
 * Authors: Rowan G. and Peter J.
 * Date: 2019-04-03
 */

module geom.luawrap.luagpath_utils;

import std.string;
import std.conv : to;
import util.lua;
import util.lua_service;
import geom;
import geom.luawrap.luageom;
import geom.luawrap.luagpath;

extern(C) int optimiseBezierPoints(lua_State* L)
{
    int narg = lua_gettop(L);
    if (narg < 2) {
        string errMsg = "Error in call to optimiseBezierPoints():\n" ~
            "A least two arguments are required: an array of points and\n" ~
            "a desired number of control points.";
        luaL_error(L, errMsg.toStringz);
    }
    // Grab points
    if (!lua_istable(L, 1)) {
        string errMsg = "Error in call to optimiseBezierPoints():\n" ~
            "An array of points (as Vector3s) is expected as argument 1,\n"~
            "but no array was found.";
        luaL_error(L, errMsg.toStringz);
    }
    int nPts = to!int(lua_objlen(L, 1));
    Vector3[] pts;
    foreach (i; 1 .. nPts+1) {
        lua_rawgeti(L, 1, i);
        pts ~= toVector3(L, -1);
        lua_pop(L, 1);
    }
    // Grab number of desired control points
    int nCtrlPts = luaL_checkint(L, 2);
    // Look for some optional arguments.
    Bezier initGuess;
    if (narg >= 3) {
        if (!lua_isnil(L, 3)) {
            initGuess = checkObj!(Bezier, BezierMT)(L, 3);
        }
    }
    double tol = 1.0e-6;
    if (narg >= 4) {
        tol = luaL_checknumber(L, 4);
    }
    int maxSteps = 10000;
    if (narg >= 5) {
        maxSteps = luaL_checkint(L, 5);
    }
    int dim = 2;
    if (narg >= 6) {
        dim = luaL_checkint(L, 6);
    }

    double[] ts; // We aren't interested in handing back t parameters
                 // through the Lua interface
    bool fittingSuccess;
    auto myBez = geom.optimiseBezierPoints(pts, nCtrlPts, initGuess, ts, fittingSuccess, tol, maxSteps, dim);
    pathStore ~= pushObj!(Bezier, BezierMT)(L, myBez);
    lua_pushboolean(L, fittingSuccess);
    return 2;
}

extern(C) int optimiseBezierPoints2(lua_State* L)
{
    int narg = lua_gettop(L);
    if (narg < 4) {
        string errMsg = "Error in call to optimiseBezierPoints2():\n" ~
            "A least four arguments are required:\n" ~
            "an array of points;\n" ~
            "a desired number of control points;\n" ~
            "a start slope; and\n" ~
            "an end slope.\n";
        luaL_error(L, errMsg.toStringz);
    }
    // Grab points
    if (!lua_istable(L, 1)) {
        string errMsg = "Error in call to optimiseBezierPoints():\n" ~
            "An array of points (as Vector3s) is expected as argument 1,\n"~
            "but no array was found.";
        luaL_error(L, errMsg.toStringz);
    }
    int nPts = to!int(lua_objlen(L, 1));
    Vector3[] pts;
    foreach (i; 1 .. nPts+1) {
        lua_rawgeti(L, 1, i);
        pts ~= toVector3(L, -1);
        lua_pop(L, 1);
    }
    // Grab number of desired control points
    int nCtrlPts = luaL_checkint(L, 2);
    // Grab start slope
    Vector3 startSlope = toVector3(L, 3);
    // Grab end slope
    Vector3 endSlope = toVector3(L, 4);
    // Look for some optional arguments.
    Bezier initGuess;
    if (narg >= 5) {
        if (!lua_isnil(L, 5)) {
            initGuess = checkObj!(Bezier, BezierMT)(L, 5);
        }
    }
    double tol = 1.0e-6;
    if (narg >= 6) {
        tol = luaL_checknumber(L, 6);
    }
    int maxSteps = 10000;
    if (narg >= 7) {
        maxSteps = luaL_checkint(L, 7);
    }
    int dim = 2;
    if (narg >= 8) {
        dim = luaL_checkint(L, 8);
    }

    double[] ts; // We aren't interested in handing back t parameters
                 // through the Lua interface
    bool fittingSuccess;
    auto myBez = geom.optimiseBezierPoints2(pts, nCtrlPts, initGuess, startSlope, endSlope, ts, fittingSuccess, tol, maxSteps, dim);
    pathStore ~= pushObj!(Bezier, BezierMT)(L, myBez);
    lua_pushboolean(L, fittingSuccess);
    return 2;
}


void registerGpathUtils(lua_State* L)
{
    lua_pushcfunction(L, &optimiseBezierPoints);
    lua_setglobal(L, "optimiseBezierPoints");
    lua_pushcfunction(L, &optimiseBezierPoints2);
    lua_setglobal(L, "optimiseBezierPoints2");
}
