/**
 * lua_service.d
 * Home to some commonly used idioms (collected as functions)
 * when using interfacing with Lua.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2015-01-14
 */

module util.lua_service;

import core.stdc.stdlib : exit;
import std.stdio;
import std.string;
import std.conv;
import std.algorithm;

import util.lua;

lua_State* init_lua_State()
{
    lua_State* L = luaL_newstate();
    luaL_openlibs(L);
    return L;
}

void doLuaFile(lua_State* L, string fname)
{
    if ( luaL_dofile(L, fname.toStringz) != 0 ) {
        string errMsg = to!string(lua_tostring(L, -1));
        throw new LuaInputException(errMsg);
    }
}

void warnOnStackChange(lua_State* L, int old_top, string file=__FILE__, size_t line=__LINE__)
{
    int new_top = lua_gettop(L);
    if (new_top - old_top != 0) {
        writefln("WARNING unexpected stack size change old_top=%d new_top=%d at line %d in file %s",
                 old_top, new_top, line, file);
    }
}

bool tableEmpty(lua_State* L, int tblIdx)
{
    lua_pushnil(L); // first key
    // If we are looking from the top of the stack,
    // the table is now one further down.
    if (tblIdx < 0) { tblIdx -= 1; }
    if (lua_next(L, tblIdx) != 0) {
        lua_pop(L, 2); // discard key and value
        return false;
    } else {
        return true;
    }
}

string getString(lua_State* L, int tblIdx, string key)
{
    int old_top = lua_gettop(L);
    if (!lua_istable(L, tblIdx)) {
        string errMsg = format("getString was expecting a table at stack index: %d", tblIdx);
        throw new Error(errMsg);
    }
    lua_getfield(L, tblIdx, key.toStringz);
    if ( !lua_isstring(L, -1) ) {
        string errMsg = format("A string was expected in field: %s", key);
        lua_pop(L, 1);
        throw new LuaInputException(errMsg);
    }
    string val = to!string(lua_tostring(L, -1));
    lua_pop(L, 1);
    warnOnStackChange(L, old_top);
    return val;
}

string getString(lua_State* L, string key)
{
    int old_top = lua_gettop(L);
    lua_getglobal(L, key.toStringz);
    if ( !lua_isstring(L, -1) ) {
        string errMsg = format("A string was expected in field: %s", key);
        lua_pop(L, 1);
        throw new LuaInputException(errMsg);
    }
    string val = to!string(lua_tostring(L, -1));
    lua_pop(L, 1);
    warnOnStackChange(L, old_top);
    return val;
}

string getStringWithDefault(lua_State* L, int tblIdx, string key, string defaultValue)
{
    int old_top = lua_gettop(L);
    string val = defaultValue;
    if (!lua_istable(L, tblIdx)) {
        string errMsg = format("getStringWithDefault was expecting a table at stack index: %d", tblIdx);
        throw new Error(errMsg);
    }
    lua_getfield(L, tblIdx, key.toStringz);
    if (lua_isstring(L, -1)) {
        val = to!string(lua_tostring(L, -1));
    }
    lua_pop(L, 1);
    warnOnStackChange(L, old_top);
    return val;
}

double getDouble(lua_State* L, int tblIdx, string key)
{
    int old_top = lua_gettop(L);
    if (!lua_istable(L, tblIdx)) {
        string errMsg = format("getDouble was expecting a table at stack index: %d", tblIdx);
        throw new Error(errMsg);
    }
    lua_getfield(L, tblIdx, key.toStringz);
    if ( !lua_isnumber(L, -1) ) {
        string errMsg = format("A double was expected in field: %s", key);
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    double val = lua_tonumber(L, -1);
    lua_pop(L, 1);
    warnOnStackChange(L, old_top);
    return val;
}

double getDouble(lua_State* L, string key)
{
    int old_top = lua_gettop(L);
    lua_getglobal(L, key.toStringz);
    if ( !lua_isnumber(L, -1) ) {
        string errMsg = format("A double was expected in field: %s", key);
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    double val = lua_tonumber(L, -1);
    lua_pop(L, 1);
    warnOnStackChange(L, old_top);
    return val;
}

double getDoubleWithDefault(lua_State* L, int tblIdx, string key, double defaultValue)
{
    int old_top = lua_gettop(L);
    if (!lua_istable(L, tblIdx)) {
        string errMsg = format("getDoubleWithDefault was expecting a table at stack index: %d", tblIdx);
        throw new Error(errMsg);
    }
    double val = defaultValue;
    lua_getfield(L, tblIdx, key.toStringz);
    if (lua_isnumber(L, -1)) {
        val = lua_tonumber(L, -1);
    }
    lua_pop(L, 1);
    warnOnStackChange(L, old_top);
    return val;
}

int getInt(lua_State* L, int tblIdx, string key)
{
    int old_top = lua_gettop(L);
    if (!lua_istable(L, tblIdx)) {
        string errMsg = format("getInt was expecting a table at stack index: %d", tblIdx);
        throw new Error(errMsg);
    }
    lua_getfield(L, tblIdx, key.toStringz);
    if ( !lua_isinteger(L, -1) ) {
        string errMsg = format("An integer was expected in field: %s", key);
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    int val = to!int(lua_tointeger(L, -1));
    lua_pop(L, 1);
    warnOnStackChange(L, old_top);
    return val;
}

int getInt(lua_State* L, string key)
{
    int old_top = lua_gettop(L);
    lua_getglobal(L, key.toStringz);
    if ( !lua_isinteger(L, -1) ) {
        string errMsg = format("An integer was expected in field: %s", key);
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    int val = to!int(lua_tointeger(L, -1));
    lua_pop(L, 1);
    warnOnStackChange(L, old_top);
    return val;
}

int getIntWithDefault(lua_State* L, int tblIdx, string key, int defaultValue)
{
    int old_top = lua_gettop(L);
    if (!lua_istable(L, tblIdx)) {
        string errMsg = format("getIntWithDefault was expecting a table at stack index: %d", tblIdx);
        throw new Error(errMsg);
    }
    int val = defaultValue;
    lua_getfield(L, tblIdx, key.toStringz);
    if (lua_isinteger(L, -1)) {
        val = to!int(lua_tointeger(L, -1));
    }
    lua_pop(L, 1);
    warnOnStackChange(L, old_top);
    return val;
}

bool getBool(lua_State* L, int tblIdx, string key)
{
    int old_top = lua_gettop(L);
    if (!lua_istable(L, tblIdx)) {
        string errMsg = format("getBool was expecting a table at stack index: %d", tblIdx);
        throw new Error(errMsg);
    }
    lua_getfield(L, tblIdx, key.toStringz);
    if ( !lua_isboolean(L, -1) ) {
        string errMsg = format("A boolean value was expected in field: %s", key);
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    bool val = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);
    warnOnStackChange(L, old_top);
    return val;
}

bool getBool(lua_State* L, string key)
{
    int old_top = lua_gettop(L);
    lua_getglobal(L, key.toStringz);
    if ( !lua_isboolean(L, -1) ) {
        string errMsg = format("A boolean value was expected in field: %s", key);
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    bool val = to!bool(lua_toboolean(L, -1));
    lua_pop(L, 1);
    warnOnStackChange(L, old_top);
    return val;
}

bool getBoolWithDefault(lua_State* L, int tblIdx, string key, bool defaultValue)
{
    int old_top = lua_gettop(L);
    if (!lua_istable(L, tblIdx)) {
        string errMsg = format("getBoolWithDefault was expecting a table at stack index: %d", tblIdx);
        throw new Error(errMsg);
    }
    bool val = defaultValue;
    lua_getfield(L, tblIdx, key.toStringz);
    if (lua_isboolean(L, -1)) {
        val = to!bool(lua_toboolean(L, -1));
    }
    lua_pop(L, 1);
    warnOnStackChange(L, old_top);
    return val;
}

void getArrayOfStrings(lua_State* L, int tblIdx, string key, out string[] values)
{
    int old_top = lua_gettop(L);
    if (!lua_istable(L, tblIdx)) {
        string errMsg = format("getArrayOfStrings was expecting a table at stack index: %d", tblIdx);
        throw new Error(errMsg);
    }
    values.length = 0;
    lua_getfield(L, tblIdx, key.toStringz);
    if ( !lua_istable(L, -1) ) {
        string errMsg = format("A table of strings was expected in field: %s", key);
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    auto n = to!int(lua_objlen(L, -1));
    foreach ( i; 1..n+1 ) {
        lua_rawgeti(L, -1, i);
        if ( lua_isstring(L, -1) ) values ~= to!string(lua_tostring(L, -1));
        // Silently ignore anything that isn't a string value.
        lua_pop(L, 1);
    }
    lua_pop(L, 1);
    warnOnStackChange(L, old_top);
}

void getArrayOfStrings(lua_State* L, string key, out string[] values)
{
    int old_top = lua_gettop(L);
    values.length = 0;
    lua_getglobal(L, key.toStringz);
    if ( !lua_istable(L, -1) ) {
        string errMsg = format("A table of strings was expected in field: %s", key);
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    auto n = to!int(lua_objlen(L, -1));
    foreach ( i; 1..n+1 ) {
        lua_rawgeti(L, -1, i);
        if ( lua_isstring(L, -1) ) values ~= to!string(lua_tostring(L, -1));
        // Silently ignore anything that isn't a string value.
        lua_pop(L, 1);
    }
    lua_pop(L, 1);
    warnOnStackChange(L, old_top);
}

/**
 * Get named array of numbers from index in Lua stack.
 */

void getArrayOfDoubles(T)(lua_State* L, int tblIdx, string key, out T[] values)
    // NOTE: This checks if we can cast T -> double, which isn't exactly what
    //       we want, but it's a reasonable assumption for now. 
    //       A better method will check if `to!` is valid, but `to!` is *very* powerful
    if (canCastTo!(T, double))
{
    int old_top = lua_gettop(L);
    if (!lua_istable(L, tblIdx)) {
        string errMsg = format("getArrayOfDoubles was expecting a table at stack index: %d", tblIdx);
        throw new Error(errMsg);
    }
    values.length = 0;
    lua_getfield(L, tblIdx, key.toStringz);
    if ( !lua_istable(L, -1) ) {
        string errMsg = format("A table of numbers was expected in field: %s", key);
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    auto n = to!int(lua_objlen(L, -1));
    foreach ( i; 1..n+1 ) {
        lua_rawgeti(L, -1, i);
        if ( lua_isnumber(L, -1) ) values ~= to!(T)(lua_tonumber(L, -1));
        // Silently ignore anything that isn't a value.
        lua_pop(L, 1);
    }
    lua_pop(L, 1);
    warnOnStackChange(L, old_top);
}

/**
 * Get named array of numbers from index in Lua stack.
 */

void getArrayOfInts(lua_State* L, int tblIdx, string key, out int[] values)
{
    int old_top = lua_gettop(L);
    if (!lua_istable(L, tblIdx)) {
        string errMsg = format("getArrayOfInts was expecting a table at stack index: %d", tblIdx);
        throw new Error(errMsg);
    }
    values.length = 0;
    lua_getfield(L, tblIdx, key.toStringz);
    if ( !lua_istable(L, -1) ) {
        string errMsg = format("A table of integers was expected in field: %s", key);
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    auto n = to!int(lua_objlen(L, -1));
    foreach ( i; 1..n+1 ) {
        lua_rawgeti(L, -1, i);
        if ( lua_isinteger(L, -1) ) values ~= to!int(lua_tointeger(L, -1));
        // Silently ignore anything that isn't a value.
        lua_pop(L, 1);
    }
    lua_pop(L, 1);
    warnOnStackChange(L, old_top);
}

void getAssocArrayOfDoubles(T)(lua_State* L, string key, string[] pList, out T[string] params)
{
    int old_top = lua_gettop(L);
    if (!lua_istable(L, -1)) {
        string errMsg = format("getAssocArrayOfDoubles was expecting a table at stack index: %d", -1);
        throw new Error(errMsg);
    }
    lua_getfield(L, -1, key.toStringz);
    if ( !lua_istable(L, -1) ) {
        string errMsg = format("A table with key and values (doubles) was expected in field: %s", key);
        lua_pop(L, 1);
        throw new Error(errMsg);
    }
    foreach (p; pList) {
        params[p] = to!(T)(getDouble(L, -1, p));
    }
    lua_pop(L, 1);
    warnOnStackChange(L, old_top);
}

// Get matrix of numbers from table at a particular stack index.

void getMatrixOfDoubles(size_t N)(lua_State* L, int tblIdx, out double[N][N] values)
{
    int old_top = lua_gettop(L);
    if (!lua_istable(L, tblIdx)) {
        string errMsg = format("getArrayOfDoubles was expecting a table at stack index: %d", tblIdx);
        throw new Error(errMsg);
    }
    foreach (i; 0 .. N) {
        lua_rawgeti(L, -1, to!int(i+1));
        foreach (j; 0 .. N) {
            lua_rawgeti(L, -1, to!int(j+1));
            values[i][j] = (lua_isnumber(L, -1)) ? luaL_checknumber(L, -1) : 0.0;
            lua_pop(L, 1);
        }
        lua_pop(L, 1);
    }
    // Table is left on stack but no other items should be added to stack.
    warnOnStackChange(L, old_top);
}

// Push a matrix of numbers onto the stack as a table of tables.

void pushMatrixOfDoubles(size_t N)(lua_State* L, double[N][N] values)
{
    lua_newtable(L); // for whole matrix
    int old_top = lua_gettop(L);
    foreach (i; 0 .. N) {
        // Build a new row.
        lua_newtable(L);
        foreach (j; 0 .. N) {
            lua_pushnumber(L, values[i][j]);
            lua_rawseti(L, -2, to!int(j+1));
        }
        // Completed row sitting at top of stack.
        lua_rawseti(L, -2, to!int(i+1));
    }
    // Table is left on stack but we have cleaned up otherwise.
    warnOnStackChange(L, old_top);
}


/**
 * This creates a new userdata spot on the Lua stack and
 * populates it with an object of type T.
 *
 * Notes:
 * (1) We explicitly set the metatable name as
 * a string and do NOT try to get the type name
 * using T.stringof (or a more elaborate __traits function).
 * The reason is that we want control of the name that
 * appears in the Lua script, no matter what name the various D
 * compilers decide to give your class type.
 * (2) We also want to keep a reference to the object in the D domain
 * so that the D garbage collector will not try to remove it and
 * any of its internally-referenced objects from that domain.
 * We cannot tell how long the Lua interpreter will want
 * to have access to the object.
 */

const(T) pushObj(T, string metatableName)(lua_State* L, T obj)
{
    auto ptr = cast(T*) lua_newuserdata(L, obj.sizeof);
    *ptr = obj; // Make a copy of the object into the Lua domain.
    luaL_getmetatable(L, metatableName.toStringz);
    lua_setmetatable(L, -2);
    return obj; // Back in the D domain, we should keep a reference.
    // It is the responsibility of the calling code to do so.
}

/**
 * Retrieve reference to an object from lua_State.
 *
 * This function looks for an object of type T
 * at the specified index in lua_State. A error is
 * raised by luaL_checkudata if the object is not of type
 * metaTableName.
 */

T checkObj(T, string metatableName)(lua_State* L, int index)
{
    if (lua_isnil(L, index)) { luaL_error(L, "nil value as pointer to type %s", metatableName.toStringz); }
    auto ptr = cast(T*) luaL_checkudata(L, index, metatableName.toStringz);
    if (ptr is null) { luaL_error(L, "null pointer to type %s", metatableName.toStringz); }
    return *ptr;
}

/**
 * Test if the Lua-object type from index position in stack matches tname.
 *
 * Assuming the object is stored as Lua userdata, this returns
 * the metatable name, which is essentially the object's type.
 * If there is no valid object, the string "nil" is returned.
 *
 * This code basically duplicates the function luaL_testudata
 * in the Lua source file: lauxlib.c. IN LUA VERSION 5.2
 * (You won't find it in our code collection.)
 * The difference is that it only looks to see if the
 * metatable matches.
 */

bool isObjType(lua_State* L, int index, string tname)
{
    bool result;
    if (tname == "number") {
        if (lua_isnumber(L, index)) {
            return true;
        } else {
            return false;
        }
    }
    if (tname == "string") {
        if (lua_isstring(L, index)) {
            return true;
        } else {
            return false;
        }
    }
    void *p = lua_touserdata(L, index);
    if ( p ) {  // value is a userdata?
        if (lua_getmetatable(L, index)) {  // does it have a metatable?
            luaL_getmetatable(L, tname.toStringz);  // get correct metatable
            if ( lua_rawequal(L, -1, -2) )  // the same?
                result = true;
            else
                result = false;
            lua_pop(L, 2);  // remove both metatables
            return result;
        }
    }
    return false;  // value is not a userdata with a metatable
}

/**
 * Call an object's toString method and push result on Lua stack.
 */

extern(C) int toStringObj(T, string metatableName)(lua_State* L)
{
    if (lua_isnil(L, 1)) { luaL_error(L, "nil value as pointer to type %s", metatableName.toStringz); }
    auto path = checkObj!(T, metatableName)(L, 1);
    if (path is null) { luaL_error(L, "null pointer to type %s", metatableName.toStringz); }
    lua_pushstring(L, toStringz(path.toString));
    return 1;
}

/**
 * Attempt to retrieve a double from a field in a table.
 *
 * We can configure what happens in the event the value is missing or invalid.
 * The boolean errorIfNotFound controls what happens if the field is missing
 * (or nil -- we can't actually tell the difference in Lua). If this is set
 * to true, then we raise a Lua error and print the supplied error message.
 * If we don't care if it's missing, then we simply return the valIfError value.
 *
 * Another common case is that we don't care if the value is missing (so set
 * errorIfNotFound to false) BUT, it is set, then we want to make sure it
 * is valid. In other words, we don't want to ignore an incorrectly set value.
 * In this case, set errorIfNotValid to true. When an invalid value is encountered,
 * we will raise a Lua error using the supplied error message.
 *
 * If we really don't care if the attempt fails, then both bools should be set
 * to false. In this case, the value valIfError is returned and no errors are raised.
 */

double getNumberFromTable(lua_State* L, int index, string field,
                          bool errorIfNotFound=false, double valIfError=double.init,
                          bool errorIfNotValid=false, string errMsg="")
{
    int old_top = lua_gettop(L);
    lua_getfield(L, index, field.toStringz);
    if ( lua_isnil(L, -1) ) {
        if ( errorIfNotFound ) {
            luaL_error(L, errMsg.toStringz);
        }
        else { // We didn't really care
            lua_pop(L, 1);
            warnOnStackChange(L, old_top);
            return valIfError;
        }
    }
    // Presumably then we have something to look at.
    if ( lua_isnumber(L, -1) ) {
        auto val = lua_tonumber(L, -1);
        lua_pop(L, 1);
        warnOnStackChange(L, old_top);
        return val;
    }
    // else, failed to find a number value.
    if ( errorIfNotValid ) {
        luaL_error(L, errMsg.toStringz);
    }
    // We didn't want to fail, so give back double.init
    lua_pop(L, 1);
    warnOnStackChange(L, old_top);
    return valIfError;
} // end getNumberFromTable()



int getIntegerFromTable(lua_State* L, int index, string field,
                        bool errorIfNotFound=false, int valIfError=int.init,
                        bool errorIfNotValid=false, string errMsg="")
{
    int old_top = lua_gettop(L);
    lua_getfield(L, index, field.toStringz);
    if ( lua_isnil(L, -1) ) {
        if ( errorIfNotFound ) {
            luaL_error(L, errMsg.toStringz);
        }
        else { // We didn't really care
            lua_pop(L, 1);
            warnOnStackChange(L, old_top);
            return valIfError;
        }
    }
    // Presumably then we have something to look at.
    if ( lua_isnumber(L, -1) ) {
        auto val = to!int(lua_tonumber(L, -1));
        lua_pop(L, 1);
        warnOnStackChange(L, old_top);
        return val;
    }
    // else, failed to find an integer value.
    if ( errorIfNotValid ) {
        luaL_error(L, errMsg.toStringz);
    }
    // We didn't want to fail, so give back int.init
    warnOnStackChange(L, old_top);
    return valIfError;
} // end getIntegerFromTable()

bool getBooleanFromTable(lua_State* L, int index, string field,
                         bool errorIfNotFound=false, bool valIfError=bool.init,
                         bool errorIfNotValid=false, string errMsg="")
{
    int old_top = lua_gettop(L);
    lua_getfield(L, index, field.toStringz);
    if ( lua_isnil(L, -1) ) {
        if ( errorIfNotFound ) {
            luaL_error(L, errMsg.toStringz);
        }
        else { // We didn't really care
            lua_pop(L, 1);
            warnOnStackChange(L, old_top);
            return valIfError;
        }
    }
    // Presumably then we have something to look at.
    if ( lua_isboolean(L, -1) ) {
        auto val = lua_toboolean(L, -1);
        lua_pop(L, 1);
        warnOnStackChange(L, old_top);
        return val;
    }
    // else, failed to find a boolean value.
    if ( errorIfNotValid ) {
        luaL_error(L, errMsg.toStringz);
    }
    // We didn't want to fail, so give back bool.init
    warnOnStackChange(L, old_top);
    return valIfError;
} // end getBooleanFromTable()

bool checkAllowedNames(lua_State* L, int tblIndx, string[] allowedNames)
{
    int old_top = lua_gettop(L);
    bool namesOk = true;
    // Iterate through the table and check each key (name).
    lua_pushnil(L); // first key
    while (lua_next(L, tblIndx) != 0) {
        // next key is left on stack at -2, value at -1
        string key = to!string(lua_tostring(L, -2));
        if (!canFind(allowedNames, key)) {
            writeln("checkAllowedNames found invalid key: ", key);
            namesOk = false;
            // We stop on the first invalid key because continuing on may lead
            // to a failure message for next, and loss of the location information.
            // The resulting error message is really confusing to the user.
            break;
        }
        lua_pop(L, 1); // discard value but keep key for next
    } // end while
    warnOnStackChange(L, old_top);
    return namesOk;
} // end checkAllowedNames()

/**
 * Custom exception type for signalling Lua input errors.
 *
 * Recipe for custom execptions used from D Cookbook, p. 21.
 *
 * Reference:
 * Ruppe, A.D. (2014)
 * D Cookbook
 * Packt Publishing, Birmingham, UK
 */
class LuaInputException : Exception {
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}


/*
unittest
{
    auto lua = new LuaState;
    auto t = lua.newTable();
    t.set!(string,double)("A", 9.0);
    t.set!(string,double)("B", -15.8);
    t.set!(string,int)("C", 2);

    string[2] keys = ["A", "B"];
    string[1] keys2 = ["C"];
    string[3] keys3 = ["A", "B", "C"];

    double[string] vals;

    /// Test 1. Grab A and B as doubles.
    getValues(t, keys, vals, "test");
    assert(vals["A"] == 9.0);
    assert(vals["B"] == -15.8);

    /// Test 2. Grab C as int.
    int[string] vals2;
    getValues(t, keys2, vals2, "test2");
    assert(vals2["C"] == 2);

    /// Test 3. Grab A and B as doubles using slice from keys3
    getValues(t, keys3[0..2], vals, "test3");
    assert(vals["A"] == 9.0);
    assert(vals["B"] == -15.8);

    /// Test 4. Grab all values as doubles
    getValues(t, keys3, vals, "test4");
    assert(vals["A"] == 9.0);
    assert(vals["B"] == -15.8);
    assert(vals["C"] == 2.0);

    /// Test 5. Expect an exit exception when we go for an invalid key.
    keys3[0] = "AA";
    try {
        getValues(t, keys3, vals, "test6");
    }
    catch (Exception e) {
        assert(e);
    }
}
*/
