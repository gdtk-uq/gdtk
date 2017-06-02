/**
 * A Lua interface for the D univariatefunctions module.
 *
 * Authors: Rowan G. and Peter J.
 * Date: 2015-02-26
 */

module luaunifunction;

import std.stdio;
import std.string;
import std.conv;
import util.lua;
import util.lua_service;
import univariatefunctions;

// Name of metatables
immutable string LinearFunctionMT = "LinearFunction";
immutable string RobertsFunctionMT = "RobertsFunction";
immutable string LuaFnClusteringMT = "LuaFnClustering";

static const(UnivariateFunction)[] functionStore;

UnivariateFunction checkUnivariateFunction(lua_State* L, int index) {
    if ( isObjType(L, index, LinearFunctionMT) ) {
	return checkObj!(LinearFunction, LinearFunctionMT)(L, index);
    }
    if ( isObjType(L, index, RobertsFunctionMT) ) {
	return checkObj!(RobertsFunction, RobertsFunctionMT)(L, index);
    }
    if ( isObjType(L, index, LuaFnClusteringMT) ) {
	return checkObj!(LuaFnUnivariateFunction, LuaFnClusteringMT)(L, index);
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
    if (!checkAllowedNames(L, 1, ["t0", "t1"])) {
	string errMsg = "Error in call to LinearFunction:new{}. Invalid name in table.";
	luaL_error(L, errMsg.toStringz);
    }
    //
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
    if (!checkAllowedNames(L, 1, ["end0", "end1", "beta"])) {
	string errMsg = "Error in call to RobertsFunction:new{}. Invalid name in table.";
	luaL_error(L, errMsg.toStringz);
    }
     //
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
    lua_remove(L, 1); // remove first arugment "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
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

}
    






