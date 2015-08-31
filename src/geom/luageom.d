/**
 * An Lua interface for the D geom module.
 *
 * This follows and adapts the examples given
 * in PIL in Chapter 28, specifically Section 28.3.
 *
 * Reference:
 * Ierusalimschy, R. (2006)
 * Programming in Lua, 2nd Edition
 * Lua.org, Rio de Janeiro 
 *
 * Authors: Rowan G. and Peter J.
 * Date: 2015-02-21
 */

module luageom;

// We cheat to get the C Lua headers by using LuaD.
import util.lua;
import std.stdio;
import std.string;
import geom;

immutable string Vector3MT = "Vector3"; // Name of Vector3 metatable

/**
 * This creates a new userdata spot on the stack
 * and populates it with a Vector3 struct.
 */
int pushVector3(lua_State *L, Vector3 vec)
{
    auto vPtr = cast(Vector3*) lua_newuserdata(L, vec.sizeof);
    *vPtr = vec;
    luaL_getmetatable(L, Vector3MT.toStringz);
    lua_setmetatable(L, -2);
    return 1;
}

/**
 * This function will serve as our "constructor"
 * in the Lua script.
 *
 * Construction from Lua can be any of:
 * ----------------------
 * a = Vector3:new(0.0, 1.0, 2.0)
 * b = Vector3:new{0.0, 1.0, 2.0}
 * c = Vector3:new{x=0.0, y=1.0, z=2.0}
 * d = Vector3:new{1.0, 3.0, x=5.0, z=8.0}
 * assert(d:x() == 1.0); assert(d:y() == 3.0); assert(d:z() == 0.0)
 * ----------------------
 * For any of the lists of arguments, missing values
 * are set to 0.0.
 * Note that if you try to mix-n-match in the table, then
 * the array-style of setting wins.
 * This constructor is fairly robust to bad parameters.
 * What will happen is that they are ignored and you get a 0.0.
 */
extern(C) int newVector3(lua_State *L)
{
    auto vec = Vector3(0.0, 0.0, 0.0);
    /* This is where we decide how the user will instantiate
     * an object in Lua-land.
     */
    lua_remove(L, 1); // remove first argument "this".

    int narg = lua_gettop(L);
    if ( narg == 1 ) {	// Could be a table or a single double value
	if ( lua_isnumber(L, 1) )  vec.refx = luaL_checknumber(L, 1);
	else if ( lua_istable(L, 1) ) {
	    // If it has a length > 0, then it's been populated array style.
	    // This style of setting beats any fields that are present.
	    size_t n = lua_objlen(L, 1);
	    if ( n >= 1 ) {
		lua_rawgeti(L, 1, 1);
		if ( lua_isnumber(L, -1) ) vec.refx = lua_tonumber(L, -1);
		lua_pop(L, 1);
	    }
	    if ( n >= 2 ) {
		lua_rawgeti(L, 1, 2);
		if ( lua_isnumber(L, -1) ) vec.refy = lua_tonumber(L, -1);
		lua_pop(L, 1);
	    }
	    if ( n >= 3 ) {
		lua_rawgeti(L, 1, 3);
		if ( lua_isnumber(L, -1) ) vec.refz = lua_tonumber(L, -1);
		lua_pop(L, 1);
	    }
	    if ( n == 0 ) { // then field based table.
		lua_getfield(L, 1, "x");
		if ( lua_isnumber(L, -1) ) vec.refx = lua_tonumber(L, -1);
		lua_pop(L, 1);
		lua_getfield(L, 1, "y");
		if ( lua_isnumber(L, -1) ) vec.refy = lua_tonumber(L, -1);
		lua_pop(L, 1);
		lua_getfield(L, 1, "z");
		if ( lua_isnumber(L, -1) ) vec.refz = lua_tonumber(L, -1);
		lua_pop(L, 1);
	    }
	}
	// else: You've given us something funny, so you're going to get
	// a Vector3(0.0, 0.0, 0.0)
    }
    else if ( narg == 2 ) {
	if ( lua_isnumber(L, 1) )  vec.refx = luaL_checknumber(L, 1);
	if ( lua_isnumber(L, 2) )  vec.refy = luaL_checknumber(L, 2);
    }
    else if ( narg >= 3 ) {
	if ( lua_isnumber(L, 1) )  vec.refx = luaL_checknumber(L, 1);
	if ( lua_isnumber(L, 2) )  vec.refy = luaL_checknumber(L, 2);
	if ( lua_isnumber(L, 3) )  vec.refz = luaL_checknumber(L, 3);
    }

    /* Regardless of how we filled in vec. We are now
     * ready to grab a piece of the lua stack and
     * place our new Vector3 there as userdata.
     */
    return pushVector3(L, vec);
}

/**
 * Provides a sanity check that the raw userdata
 * is in fact what we think it is.
 */
Vector3* checkVector3(lua_State *L, int index)
{
    auto vPtr = cast(Vector3*) luaL_checkudata(L, index, Vector3MT.toStringz);
    return vPtr;
}

/*-------- exposed Vector3 methods ------------ */

// The x(), y(), z() methods are a little funny
// because they act as both getters and setters.
// We are faking data access in a sense.
/**
 * Acts as both getter and setter for x component of Vector3.
 *
 * Example:
 * -------------------------------
 * a = Vector3()
 * a:x(0.8) -- used as a setter
 * b = a:x() -- used as a getter
 * -------------------------------
 */
extern(C) int xVector3(lua_State* L)
{
    int narg = lua_gettop(L);
    auto a = checkVector3(L, 1);
    if ( narg == 1 ) { // This is a getter
	lua_pushnumber(L, a.x);
	return 1;
    }
    // else: treat as a setter.
    a.refx = luaL_checknumber(L, 2);
    return 0;
}

/**
 * Acts as both a getter and setter for y component of Vector3.
 *
 * See example for xVector3()
 */
extern(C) int yVector3(lua_State* L)
{
    int narg = lua_gettop(L);
    auto a = checkVector3(L, 1);
    if( narg == 1 ) { // This is a getter
	lua_pushnumber(L, a.y);
	return 1;
    }
    // else: treat as a setter.
    a.refy = luaL_checknumber(L, 2);
    return 0;
}

/**
 * Acts as both a getter and setter for y component of Vector3.
 *
 * See example for xVector3()
 */
extern(C) int zVector3(lua_State* L)
{
    int narg = lua_gettop(L);
    auto a = checkVector3(L, 1);
    if( narg == 1 ) { // This is a getter
	lua_pushnumber(L, a.z);
	return 1;
    }
    // else: treat as a setter.
    a.refz = luaL_checknumber(L, 2);
    return 0;
}

/**
 * This provied the unary minus operator for Lua.
 */
extern(C) int opUnaryMinVector3(lua_State* L)
{
    auto a = checkVector3(L, 1);
    auto b = -(*a);
    return pushVector3(L, b);
}

/**
 * Adds two Vector3 objects. Exposes geom.Vector3.opBinary("+")
 */
extern(C) int addVector3(lua_State* L)
{
    auto a = checkVector3(L, 1);
    auto b = checkVector3(L, 2);
    auto c = (*a) + (*b);
    return pushVector3(L, c);
}

/**
 * Subtracts two Vector3 objects. Exposes geom.Vector3.opBinary("-")
 */
extern(C) int subVector3(lua_State *L)
{
    auto a = checkVector3(L, 1);
    auto b = checkVector3(L, 2);
    auto c = (*a) - (*b);
    return pushVector3(L, c);
}

/**
 * Multiplies a Vector3 object by scalar.
 */
extern(C) int mulVector3(lua_State *L)
{
    // Need to test order of arguments.
    // Could be:
    //   Vector3 * val  <or>
    //   val * Vector3
    Vector3* a;
    double b;
    if ( lua_isuserdata(L, 1) ) {
	if ( !lua_isnumber(L, 2) ) {
	    string errMsg = "can't multiply Vector3 by non-number";
	    luaL_error(L, errMsg.toStringz);
	}
	a = checkVector3(L, 1);
	b = luaL_checknumber(L, 2);
    }
    else {
	if ( !lua_isnumber(L, 1) ) {
	    string errMsg = "can't multiply Vector3 by non-number";
	    luaL_error(L, errMsg.toStringz);
	}   
	a = checkVector3(L, 2);
	b = luaL_checknumber(L, 1);
    }
    auto c = (*a) * b;
    return pushVector3(L, c);
}

/**
 * Divides a Vector3 by a scalar.
 */
extern(C) int divVector3(lua_State* L)
{
    /* Lua could pass us:
     *     Vector3 / scalar <or>
     *     scalar / Vector3
     * We can't do anything with the
     * second form, so signal an error.
     */
    Vector3* a;
    double b;
    if ( lua_isuserdata(L, 1) ) {
	if ( !lua_isnumber(L, 2) ) {
	    string errMsg = "can't divide Vector3 by non-number";
	    luaL_error(L, errMsg.toStringz);
	}
	a = checkVector3(L, 1);
	b = luaL_checknumber(L, 2);
    }
    else {
	string errMsg = "can't divide by a Vector3 object";
	luaL_error(L, errMsg.toStringz);
    }
    auto c = (*a) / b;
    return pushVector3(L, c);
}

/**
 * Normalizes a Vector3 object. Exposes geom.Vector3.normalize()
 */
extern(C) int normalizeVector3(lua_State* L)
{
    auto a = checkVector3(L, 1);
    a.normalize();
    return 0;
}

/**
 * Computes the dot product of two Vector3s.
 *
 * Note that this Lua function can service
 * as both the method form and function form
 * of the D dot function.
 */
extern(C) int dotVector3(lua_State* L)
{
    auto a = checkVector3(L, 1);
    auto b = checkVector3(L, 2);
    lua_pushnumber(L, dot(*a, *b));
    return 1;
}

/**
 * Computes the magnitude of a Vector3.
 */
extern(C) int absVector3(lua_State* L)
{
    auto a = checkVector3(L, 1);
    lua_pushnumber(L, abs(*a));
    return 1;
}

/**
 * Returns the unit vector in the same direction.
 */
extern(C) int unitVector3(lua_State *L)
{
    auto a = checkVector3(L, 1);
    // Copy before normalizing so that
    // we don't change 'a'.
    auto b = *a;
    b.normalize();
    return pushVector3(L, b); 
}

/**
 * Returns the cross product of two Vector3s.
 */
extern(C) int crossVector3(lua_State* L)
{
    auto a = checkVector3(L, 1);
    auto b = checkVector3(L, 2);
    return pushVector3(L, cross(*a, *b));
}

extern(C) int toStringVector3(lua_State* L)
{
    auto a = checkVector3(L, 1);
    lua_pushstring(L, a.toString.toStringz);
    return 1;
}

void registerVector3(lua_State* L)
{
    luaL_newmetatable(L, Vector3MT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newVector3);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &xVector3);
    lua_setfield(L, -2, "x");
    lua_pushcfunction(L, &yVector3);
    lua_setfield(L, -2, "y");
    lua_pushcfunction(L, &zVector3);
    lua_setfield(L, -2, "z");
    lua_pushcfunction(L, &toStringVector3);
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &opUnaryMinVector3);
    lua_setfield(L, -2, "__unm");
    lua_pushcfunction(L, &addVector3);
    lua_setfield(L, -2, "__add");
    lua_pushcfunction(L, &subVector3);
    lua_setfield(L, -2, "__sub");
    lua_pushcfunction(L, &mulVector3);
    lua_setfield(L, -2, "__mul");
    lua_pushcfunction(L, &divVector3);
    lua_setfield(L, -2, "__div");
    lua_pushcfunction(L, &normalizeVector3);
    lua_setfield(L, -2, "normalize");
    lua_pushcfunction(L, &dotVector3);
    lua_setfield(L, -2, "dot");
    lua_pushcfunction(L, &absVector3);
    lua_setfield(L, -2, "abs");
    lua_pushcfunction(L, &unitVector3);
    lua_setfield(L, -2, "unit");

    lua_setglobal(L, Vector3MT.toStringz);

    /* Also attempt to put "add" in the global namespace. */
    lua_pushcfunction(L, &addVector3);
    lua_setglobal(L, "add");
    lua_pushcfunction(L, &dotVector3);
    lua_setglobal(L, "dot");
    lua_pushcfunction(L, &absVector3);
    lua_setglobal(L, "vabs"); // to avoid name clash with math library
    lua_pushcfunction(L, &unitVector3);
    lua_setglobal(L, "unit");
    lua_pushcfunction(L, &crossVector3);
    lua_setglobal(L, "cross");
}
    
