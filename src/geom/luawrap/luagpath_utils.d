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
        pts ~= *(checkVector3(L, -1));
        lua_pop(L, 1);
    }
    // Grab number of desired control points
    int nCtrlPts = luaL_checkint(L, 2);
    // Look for some optional arguments.
    double tol = 1.0e-6;
    if (narg >= 3) {
        tol = luaL_checknumber(L, 3);
    }
    int maxSteps = 10000;
    if (narg >= 4) {
        maxSteps = luaL_checkint(L, 4);
    }
    int dim = 2;
    if (narg >= 5) {
        dim = luaL_checkint(L, 5);
    }

    double[] ts; // We aren't interested in handing back t parameters
                 // through the Lua interface
    
    auto myBez = geom.optimiseBezierPoints(pts, nCtrlPts, ts, tol, maxSteps, dim);
    pathStore ~= pushObj!(Bezier, BezierMT)(L, myBez);
    return 1;
}

void registerGpathUtils(lua_State* L)
{
    lua_pushcfunction(L, &optimiseBezierPoints);
    lua_setglobal(L, "optimiseBezierPoints");
}
