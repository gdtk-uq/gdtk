/**
 * A Lua interface for the D univariatefunctions module.
 *
 * Authors: Rowan G. and Peter J.
 * Date: 2015-02-26
 */

module luaunifunction;

// We cheat to get the C Lua headers by using LuaD.
import std.stdio;
import std.string;
import util.lua;
import util.lua_service;
import univariatefunctions;

// Name of metatables
immutable string LinearFunctionMT = "LinearFunction";
immutable string RobertsFunctionMT = "RobertsFunction";

static const(UnivariateFunction)[] functionStore;

UnivariateFunction checkUnivariateFunction(lua_State* L, int index) {
    if ( isObjType(L, index, LinearFunctionMT) ) {
	return checkObj!(LinearFunction, LinearFunctionMT)(L, index);
    }
    if ( isObjType(L, index, RobertsFunctionMT) ) {
	return checkObj!(RobertsFunction, RobertsFunctionMT)(L, index);
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
    lua_remove(L, 1); // remove first argument "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = `Error in call to LinearFunction:new{}.;
A table containing arguments is expected, but no table was found.`;
	luaL_error(L, errMsg.toStringz);
    }
    string errMsgTmplt = `Error in call to LinearFunction:new{}.
A valid value for '%s' was not found in list of arguments.
The value, if present, should be a number.`;
    double t0 = getNumberFromTable(L, 1, "t0", false, 0.0, true, format(errMsgTmplt, "t0"));
    double t1 = getNumberFromTable(L, 1, "t1", false, 1.0, true, format(errMsgTmplt, "t1"));
    auto f = new LinearFunction(t0, t1);
    functionStore ~= pushObj!(LinearFunction, LinearFunctionMT)(L, f);
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
    lua_remove(L, 1); // remove first argument "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = `Error in call to RobertsFunction:new{}.;
A table containing arguments is expected, but no table was found.`;
	luaL_error(L, errMsg.toStringz);
    }
    string errMsgTmpltNumber = `Error in call to RobertsFunction:new{}.
A valid value for '%s' was not found in list of arguments.
The value, if present, should be a number.`;
    string errMsgTmpltBool = `Error in call to RobertsFunction:new{}.
A valid value for '%s' was not found in list of arguments.
The value, if present, should be boolean (true or false).`;
    bool end0 = getBooleanFromTable(L, 1, "end0", false, false, true, format(errMsgTmpltBool, "end0"));
    bool end1 = getBooleanFromTable(L, 1, "end1", false, false, true, format(errMsgTmpltBool, "end1"));
    double beta = getNumberFromTable(L, 1, "beta", false, 1.0, true, format(errMsgTmpltNumber, "beta"));
    auto f = new RobertsFunction(end0, end1, beta);
    functionStore ~= pushObj!(RobertsFunction, RobertsFunctionMT)(L, f);
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
}
    






