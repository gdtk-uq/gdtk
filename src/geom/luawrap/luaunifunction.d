/**
 * A Lua interface for the D univariatefunctions module.
 *
 * Authors: Rowan G. and Peter J.
 * Date: 2015-02-26
 */

module geom.luawrap.luaunifunction;

import std.stdio;
import std.string;
import std.conv;
import util.lua;
import util.lua_service;
import geom.misc.univariatefunctions;

// Name of metatables
immutable string LinearFunctionMT = "LinearFunction";
immutable string QuadraticFunctionMT = "QuadraticFunction";
immutable string RobertsFunctionMT = "RobertsFunction";
immutable string LuaFnClusteringMT = "LuaFnClustering";
immutable string GeometricFunctionMT = "GeometricFunction";
immutable string GaussianFunctionMT = "GaussianFunction";
immutable string GaussGeomHybridFunctionMT = "GaussGeomHybridFunction";

static const(UnivariateFunction)[] functionStore;

UnivariateFunction checkUnivariateFunction(lua_State* L, int index) {
    if ( isObjType(L, index, LinearFunctionMT) ) {
        return checkObj!(LinearFunction, LinearFunctionMT)(L, index);
    }
    if ( isObjType(L, index, QuadraticFunctionMT) ) {
        return checkObj!(QuadraticFunction, QuadraticFunctionMT)(L, index);
    }
    if ( isObjType(L, index, RobertsFunctionMT) ) {
        return checkObj!(RobertsFunction, RobertsFunctionMT)(L, index);
    }
    if ( isObjType(L, index, LuaFnClusteringMT) ) {
        return checkObj!(LuaFnUnivariateFunction, LuaFnClusteringMT)(L, index);
    }
    if ( isObjType(L, index, GeometricFunctionMT) ) {
        return checkObj!(GeometricFunction, GeometricFunctionMT)(L, index);
    }
    if ( isObjType(L, index, GaussianFunctionMT) ) {
        return checkObj!(GaussianFunction, GaussianFunctionMT)(L, index);
    }
    if ( isObjType(L, index, GaussGeomHybridFunctionMT) ) {
        return checkObj!(GaussGeomHybridFunction, GaussGeomHybridFunctionMT)(L, index);
    }
    // if all else fails
    return null;
}

extern(C) int opCallUnivariateFunction(T, string MTname)(lua_State* L)
{
    auto f = checkObj!(T, MTname)(L, 1);
    auto t = luaL_checknumber(L, 2);
    lua_pushnumber(L, f(t));
    return 1;
}

extern(C) int copyUnivariateFunction(T, string MTname)(lua_State* L)
{
    // Sometimes it's convenient to get a copy of a function.
    auto f = checkObj!(T, MTname)(L, 1);
    functionStore ~= pushObj!(T, MTname)(L, f);
    return 1;
}

/* ----------------- LinearFunction specific functions --------------- */

/**
 * The Lua constructor for a LinearFunction.
 *
 * Example construction in Lua:
 * --------------------------------------
 * f = LinearFunction:new{t0=0.0, t1=1.0}
 * --------------------------------------
 */
extern(C) int newLinearFunction(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected LinearFunction:new{}; ";
        errMsg ~= "maybe you tried LinearFunction.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to LinearFunction:new{}.; ";
        errMsg ~= "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["t0", "t1"])) {
        string errMsg = "Error in call to LinearFunction:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    //
    string errMsgTmplt = "Error in call to LinearFunction:new{}. ";
    errMsgTmplt ~= "A valid value for '%s' was not found in list of arguments. ";
    errMsgTmplt ~= "The value, if present, should be a number.";
    double t0 = getNumberFromTable(L, 1, "t0", false, 0.0, true, format(errMsgTmplt, "t0"));
    double t1 = getNumberFromTable(L, 1, "t1", false, 1.0, true, format(errMsgTmplt, "t1"));
    auto f = new LinearFunction(t0, t1);
    functionStore ~= pushObj!(LinearFunction, LinearFunctionMT)(L, f);
    return 1;
}

/* ----------------- QuadraticFunction specific functions --------------- */

/**
 * The Lua constructor for a QuadraticFunction.
 *
 * Example construction in Lua:
 * --------------------------------------
 * f = QuadraticFunction:new{ratio=2.0, reverse=false}
 * --------------------------------------
 */
extern(C) int newQuadraticFunction(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected QuadraticFunction:new{}; ";
        errMsg ~= "maybe you tried QuadraticFunction.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to QuadraticFunction:new{}.; ";
        errMsg ~= "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["ratio", "reverse"])) {
        string errMsg = "Error in call to QuadraticFunction:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    //
    string errMsgTmplt1 = "Error in call to QuadraticFunction:new{}. ";
    errMsgTmplt1 ~= "A valid value for '%s' was not found in list of arguments. ";
    errMsgTmplt1 ~= "The value, if present, should be a number.";
    double ratio = getNumberFromTable(L, 1, "ratio", false, 1.0, true, format(errMsgTmplt1, "ratio"));
    string errMsgTmplt2 = "Error in call to QuadraticFunction:new{}. ";
    errMsgTmplt2 ~= "A valid value for '%s' was not found in list of arguments. ";
    errMsgTmplt2 ~= "The value, if present, should be true or false.";
    bool reverse = getBooleanFromTable(L, 1, "reverse", false, false, true, format(errMsgTmplt2, "reverse"));
    auto f = new QuadraticFunction(ratio, reverse);
    functionStore ~= pushObj!(QuadraticFunction, QuadraticFunctionMT)(L, f);
    return 1;
}

/**
 * The Lua constructor for a RobertsFunction.
 *
 * Example construction in Lua:
 * -------------------------------------------------------
 * f = RobertsFunction:new{end0=true, end1=true, beta=1.0}
 * -------------------------------------------------------
 */
extern(C) int newRobertsFunction(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected RobertsFunction:new{}; ";
        errMsg ~= "maybe you tried RobertsFunction.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to RobertsFunction:new{}. ";
        errMsg ~= "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["end0", "end1", "beta"])) {
        string errMsg = "Error in call to RobertsFunction:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
     //
    string errMsgTmpltNumber = "Error in call to RobertsFunction:new{}. ";
    errMsgTmpltNumber ~= "A valid value for '%s' was not found in list of arguments. ";
    errMsgTmpltNumber ~= "The value, if present, should be a number.";
    string errMsgTmpltBool = "Error in call to RobertsFunction:new{}. ";
    errMsgTmpltBool ~= "A valid value for '%s' was not found in list of arguments. ";
    errMsgTmpltBool ~= "The value, if present, should be boolean (true or false).";
    bool end0 = getBooleanFromTable(L, 1, "end0", false, false, true, format(errMsgTmpltBool, "end0"));
    bool end1 = getBooleanFromTable(L, 1, "end1", false, false, true, format(errMsgTmpltBool, "end1"));
    double beta = getNumberFromTable(L, 1, "beta", false, 1.0, true, format(errMsgTmpltNumber, "beta"));
    auto f = new RobertsFunction(end0, end1, beta);
    functionStore ~= pushObj!(RobertsFunction, RobertsFunctionMT)(L, f);
    return 1;
}

/**
 * LuaFnClustering and its Lua constructor.
 */

class LuaFnUnivariateFunction : UnivariateFunction {
public:
    lua_State *L;
    string luaFnName;

    this(const lua_State *L, string luaFnName)
    {
        this.L = cast(lua_State*)L;
        this.luaFnName = luaFnName;
    }
    this(ref const(LuaFnUnivariateFunction) other)
    {
        L = cast(lua_State*)other.L;
        luaFnName = other.luaFnName;
    }
    LuaFnUnivariateFunction dup() const
    {
        return new LuaFnUnivariateFunction(this.L, this.luaFnName);
    }
    override double opCall(double t) const
    {
        // Call back to the Lua function.
        lua_getglobal(cast(lua_State*)L, luaFnName.toStringz);
        lua_pushnumber(cast(lua_State*)L, t);
        if ( lua_pcall(cast(lua_State*)L, 1, 1, 0) != 0 ) {
            string errMsg = "Error in call to " ~ luaFnName ~ 
                " from LuaFnClustering:opCall(): " ~ 
                to!string(lua_tostring(cast(lua_State*)L, -1));
            luaL_error(cast(lua_State*)L, errMsg.toStringz);
        }
        // We are expecting a double value to be returned.
        if ( !lua_isnumber(cast(lua_State*)L, -1) ) {
            string errMsg = "Error in call to LuaFnClustering:opCall().; " ~
                "A single floating point number is expected, but none found.";
            luaL_error(cast(lua_State*)L, errMsg.toStringz);
        }
        double u = luaL_checknumber(cast(lua_State*)L, -1);
        lua_settop(cast(lua_State*)L, 0); // clear the stack
        return u;
    } // end opCall()

    override string toString() const
    {
        return "LuaFnUnivariateFunction()";
    }
}

extern(C) int newLuaFnUnivariateFunction(lua_State *L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected RobertsFunction:new{}; ";
        errMsg ~= "maybe you tried RobertsFunction.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first arugment "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to LuaFnClustering:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["luaFnName"])) {
        string errMsg = "Error in call to LuaFnUnivariateFunction:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
    string fnName = "";
    lua_getfield(L, 1, toStringz("luaFnName"));
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to LuaFnUnivariateFunction.new{}. No luaFnName entry found.";
        luaL_error(L, errMsg.toStringz);
    }
    if ( lua_isstring(L, -1) ) {
        fnName ~= to!string(lua_tostring(L, -1));
    }
    lua_pop(L, 1);
    if ( fnName == "" ) {
        string errMsg = "Error in call to LuaFnUnivariateFunction:new{}. No function name is empty.";
        luaL_error(L, errMsg.toStringz);
    }
    auto lfc = new LuaFnUnivariateFunction(L, fnName);
    functionStore ~= pushObj!(LuaFnUnivariateFunction, LuaFnClusteringMT)(L, lfc);
    return 1;
}

/**
 * The Lua constructor for a GeometricFunction.
 *
 * Example construction in Lua:
 * --------------------------------------------------------------
 * f = GeometricFunction:new{a=0.005, r=1.1, N=40, reverse=false}
 * --------------------------------------------------------------
 */
extern(C) int newGeometricFunction(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected GeometricFunction:new{}; ";
        errMsg ~= "maybe you tried GeometricFunction.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to GeometricFunction:new{}. ";
        errMsg ~= "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["a", "r", "N", "reverse"])) {
        string errMsg = "Error in call to GeometricFunction:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
     //
    string errMsgTmpltNumber = "Error in call to GeometricFunction:new{}. ";
    errMsgTmpltNumber ~= "A valid value for '%s' was not found in list of arguments. ";
    errMsgTmpltNumber ~= "The value, if present, should be a number.";
    string errMsgTmpltBool = "Error in call to GeometricFunction:new{}. ";
    errMsgTmpltBool ~= "A valid value for '%s' was not found in list of arguments. ";
    errMsgTmpltBool ~= "The value, if present, should be boolean (true or false).";
    string errMsgTmpltInt = "Error in call to GeometricFunction:new{}. ";
    errMsgTmpltInt ~= "A valid value for '%s' was not found in list of arguments. ";
    errMsgTmpltInt ~= "The value, if present, should be an integer.";
    double a = getNumberFromTable(L, 1, "a", true, 1.0, true, format(errMsgTmpltNumber, "a"));
    double r = getNumberFromTable(L, 1, "r", true, 1.0, true, format(errMsgTmpltNumber, "r"));
    int N = getIntegerFromTable(L, 1, "N", true, 1, true, format(errMsgTmpltInt, "N"));
    bool reverse = getBooleanFromTable(L, 1, "reverse", false, false, true, format(errMsgTmpltBool, "reverse"));
    auto f = new GeometricFunction(a, r, N, reverse);
    functionStore ~= pushObj!(GeometricFunction, GeometricFunctionMT)(L, f);
    return 1;
}

/**
 * The Lua constructor for a GaussianFunction.
 *
 * Example construction in Lua:
 * --------------------------------------------------------------
 * f = GeometricFunction:new{m=0.5, s=0.1, ratio=0.2}
 * --------------------------------------------------------------
 */
extern(C) int newGaussianFunction(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected GaussianFunction:new{}; ";
        errMsg ~= "maybe you tried GaussianFunction.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to GaussianFunction:new{}. ";
        errMsg ~= "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["m", "s", "ratio"])) {
        string errMsg = "Error in call to GaussianFunction:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
     //
    string errMsgTmpltNumber = "Error in call to GaussianFunction:new{}. ";
    errMsgTmpltNumber ~= "A valid value for '%s' was not found in list of arguments. ";
    errMsgTmpltNumber ~= "The value, if present, should be a number.";
    double m = getNumberFromTable(L, 1, "m", true, 0.5, true, format(errMsgTmpltNumber, "m"));
    double s = getNumberFromTable(L, 1, "s", true, 0.1, true, format(errMsgTmpltNumber, "s"));
    double ratio = getNumberFromTable(L, 1, "ratio", true, 0.2, true, format(errMsgTmpltNumber, "ratio"));
    auto f = new GaussianFunction(m, s, ratio);
    functionStore ~= pushObj!(GaussianFunction, GaussianFunctionMT)(L, f);
    return 1;
}

/**
 * The Lua constructor for a GaussGeomHybridFunction.
 *
 * Example construction in Lua:
 * --------------------------------------------------------------
 * f = GaussGeomHybridFunction:new{A=0.01, R=1.2, N=40, m=0.8, s=0.1, ratio=0.2, reverse=false}
 * --------------------------------------------------------------
 */
extern(C) int newGaussGeomHybridFunction(lua_State* L)
{
    int narg = lua_gettop(L);
    if ( !(narg == 2 && lua_istable(L, 1)) ) {
        // We did not get what we expected as arguments.
        string errMsg = "Expected GaussGeomHybridFunction:new{}; ";
        errMsg ~= "maybe you tried GaussGeomHybridFunction.new{}.";
        luaL_error(L, errMsg.toStringz);
    }
    lua_remove(L, 1); // remove first argument "this"
    if ( !lua_istable(L, 1) ) {
        string errMsg = "Error in call to GaussGeomHybridFunction:new{}. ";
        errMsg ~= "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    if (!checkAllowedNames(L, 1, ["A","R","N","m", "s","ratio","reverse"])) {
        string errMsg = "Error in call to GaussGeomHybridFunction:new{}. Invalid name in table.";
        luaL_error(L, errMsg.toStringz);
    }
     //
    string errMsgTmpltNumber = "Error in call to GaussGeomHybridFunction:new{}. ";
    errMsgTmpltNumber ~= "A valid value for '%s' was not found in list of arguments. ";
    errMsgTmpltNumber ~= "The value, if present, should be a number.";
    string errMsgTmpltBool = "Error in call to GaussGeomHybridFunction:new{}. ";
    errMsgTmpltBool ~= "A valid value for '%s' was not found in list of arguments. ";
    errMsgTmpltBool ~= "The value, if present, should be boolean (true or false).";
    string errMsgTmpltInt = "Error in call to GaussGeomHybridFunction:new{}. ";
    errMsgTmpltInt ~= "A valid value for '%s' was not found in list of arguments. ";
    errMsgTmpltInt ~= "The value, if present, should be an integer.";

    double A = getNumberFromTable(L, 1, "A", true, 1.0, true, format(errMsgTmpltNumber, "A"));
    double R = getNumberFromTable(L, 1, "R", true, 1.0, true, format(errMsgTmpltNumber, "R"));
    int N = getIntegerFromTable(L, 1, "N", true, 1, true, format(errMsgTmpltInt, "N"));
    double m = getNumberFromTable(L, 1, "m", true, 0.5, true, format(errMsgTmpltNumber, "m"));
    double s = getNumberFromTable(L, 1, "s", true, 0.1, true, format(errMsgTmpltNumber, "s"));
    double ratio = getNumberFromTable(L, 1, "ratio", true, 0.2, true, format(errMsgTmpltNumber, "ratio"));
    bool reverse = getBooleanFromTable(L, 1, "reverse", false, false, true, format(errMsgTmpltBool, "reverse"));
    auto f = new GaussGeomHybridFunction(A, R, N, m, s, ratio, reverse);
    functionStore ~= pushObj!(GaussGeomHybridFunction, GaussGeomHybridFunctionMT)(L, f);
    return 1;
}

void registerUnivariateFunctions(lua_State* L)
{
    // Register the LinearFunction object
    luaL_newmetatable(L, LinearFunctionMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newLinearFunction);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallUnivariateFunction!(LinearFunction, LinearFunctionMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallUnivariateFunction!(LinearFunction, LinearFunctionMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(LinearFunction, LinearFunctionMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyUnivariateFunction!(LinearFunction, LinearFunctionMT));
    lua_setfield(L, -2, "copy");

    lua_setglobal(L, LinearFunctionMT.toStringz);


    // Register the QuadraticFunction object
    luaL_newmetatable(L, QuadraticFunctionMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newQuadraticFunction);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallUnivariateFunction!(QuadraticFunction, QuadraticFunctionMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallUnivariateFunction!(QuadraticFunction, QuadraticFunctionMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(QuadraticFunction, QuadraticFunctionMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyUnivariateFunction!(QuadraticFunction, QuadraticFunctionMT));
    lua_setfield(L, -2, "copy");

    lua_setglobal(L, QuadraticFunctionMT.toStringz);


    // Register the RobertsFunction object
    luaL_newmetatable(L, RobertsFunctionMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newRobertsFunction);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallUnivariateFunction!(RobertsFunction, RobertsFunctionMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallUnivariateFunction!(RobertsFunction, RobertsFunctionMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(RobertsFunction, RobertsFunctionMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyUnivariateFunction!(RobertsFunction, RobertsFunctionMT));
    lua_setfield(L, -2, "copy");

    lua_setglobal(L, RobertsFunctionMT.toStringz);

    // Register the LuaFnClustering object
    luaL_newmetatable(L, LuaFnClusteringMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newLuaFnUnivariateFunction);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallUnivariateFunction!(LuaFnUnivariateFunction, LuaFnClusteringMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallUnivariateFunction!(LuaFnUnivariateFunction, LuaFnClusteringMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(LuaFnUnivariateFunction, LuaFnClusteringMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyUnivariateFunction!(LuaFnUnivariateFunction, LuaFnClusteringMT));
    lua_setfield(L, -2, "copy");

    lua_setglobal(L, LuaFnClusteringMT.toStringz);

    // Register the GeometricFunction object
    luaL_newmetatable(L, GeometricFunctionMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newGeometricFunction);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallUnivariateFunction!(GeometricFunction, GeometricFunctionMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallUnivariateFunction!(GeometricFunction, GeometricFunctionMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(GeometricFunction, GeometricFunctionMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyUnivariateFunction!(GeometricFunction, GeometricFunctionMT));
    lua_setfield(L, -2, "copy");

    lua_setglobal(L, GeometricFunctionMT.toStringz);

    // Register the GaussianFunction object
    luaL_newmetatable(L, GaussianFunctionMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newGaussianFunction);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallUnivariateFunction!(GaussianFunction, GaussianFunctionMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallUnivariateFunction!(GaussianFunction, GaussianFunctionMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(GaussianFunction, GaussianFunctionMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyUnivariateFunction!(GaussianFunction, GaussianFunctionMT));
    lua_setfield(L, -2, "copy");

    lua_setglobal(L, GaussianFunctionMT.toStringz);

    // Register the GaussGeomHybridFunction object
    luaL_newmetatable(L, GaussGeomHybridFunctionMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newGaussGeomHybridFunction);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &opCallUnivariateFunction!(GaussGeomHybridFunction, GaussGeomHybridFunctionMT));
    lua_setfield(L, -2, "__call");
    lua_pushcfunction(L, &opCallUnivariateFunction!(GaussGeomHybridFunction, GaussGeomHybridFunctionMT));
    lua_setfield(L, -2, "eval");
    lua_pushcfunction(L, &toStringObj!(GaussGeomHybridFunction, GaussGeomHybridFunctionMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &copyUnivariateFunction!(GaussGeomHybridFunction, GaussGeomHybridFunctionMT));
    lua_setfield(L, -2, "copy");

    lua_setglobal(L, GaussGeomHybridFunctionMT.toStringz);
}
