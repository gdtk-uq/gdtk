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

import util.lua;
import std.stdio;
import std.string;
import util.lua_service;
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
 * Provides a sanity check that the raw userdata
 * is in fact what we think it is.
 */
Vector3* checkVector3(lua_State *L, int index)
{
    auto vPtr = cast(Vector3*) luaL_checkudata(L, index, Vector3MT.toStringz);
    return vPtr;
}

/**
 * This function will serve as our "constructor"
 * in the Lua script.
 *
 * Construction from Lua can be one of:
 * ----------------------
 * a = Vector3:new{otherVector3}
 * b = Vector3:new(otherVector3)
 * c = Vector3:new{x=1.0, y=1.0, z=2.0}
 * assert(c:x() == 1.0); assert(c:y() == 3.0); assert(c:z() == 2.0)
 * ----------------------
 * When a single argument is given, it may be another Vector3 object.
 * For the table of coordinates, missing values are assumed to be 0.0.
 */
extern(C) int newVector3(lua_State *L)
{
    auto vec = Vector3(0.0, 0.0, 0.0);
    lua_remove(L, 1); // remove first argument "this".
    int narg = lua_gettop(L);
    if ( narg == 1 ) {	// Could be a table or a single Vector3 object
	if ( lua_istable(L, 1) ) {
	    // If it has a length >= 1, then it's been populated array style.
	    // This style of setting beats any fields that are present.
	    if ( lua_objlen(L, 1) >= 1 ) {
		// A single item may be a Vector3 object already.
		lua_rawgeti(L, 1, 1);
		vec = *checkVector3(L, -1);
		lua_pop(L, 1);
	    } else {
		// Look for named fields containing coordinate values.
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
	} else {
	    vec = *checkVector3(L, -1);
	}
    } else {
	// Just leave the zero-filled values.
    } 
    // Regardless of how we filled in vec, we are now ready 
    // to place our new Vector3 on the stack as userdata.
    return pushVector3(L, vec);
} // end newVector3()

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
 * Moves a Vector3 object about the z-axis. 
 * Exposes geom.Vector3.rotate_about_z_axis(angle)
 */
extern(C) int rotateAboutZAxisVector3(lua_State* L)
{
    auto a = checkVector3(L, 1);
    double angle = luaL_checknumber(L, 2);
    a.rotate_about_zaxis(angle);
    return 0;
}

/**
 * Mirror-image a Vector3 object through a given plane that is
 * defined by a point and a normal. 
 * Exposes geom.Vector3.mirror_image(point, normal)
 */
extern(C) int mirrorImageVector3(lua_State* L)
{
    auto a = checkVector3(L, 1);
    auto point = checkVector3(L, 2);
    auto normal = checkVector3(L, 3);
    a.mirror_image(*point, *normal);
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

/**
 * Returns a table containing the named properties.
 * Expects a table of names corners.
 */
extern(C) int quadProperties(lua_State* L)
{
    if ( lua_istable(L, 1) ) {
	lua_getfield(L, 1, "p0");
	Vector3 p0 = *checkVector3(L, -1);
	lua_pop(L, 1);
	lua_getfield(L, 1, "p1");
	Vector3 p1 = *checkVector3(L, -1);
	lua_pop(L, 1);
	lua_getfield(L, 1, "p2");
	Vector3 p2 = *checkVector3(L, -1);
	lua_pop(L, 1);
	lua_getfield(L, 1, "p3");
        Vector3 p3 = *checkVector3(L, -1);
	lua_pop(L, 1);
	lua_settop(L, 0); // clear stack

	Vector3 centroid;
	Vector3 n;
	Vector3 t1;
	Vector3 t2;
	double area;
	quad_properties(p0, p1, p2, p3, centroid, n, t1, t2, area);

	lua_newtable(L); // anonymous table { }
	auto tblIndx = lua_gettop(L);
	pushVector3(L, centroid);
	lua_setfield(L, tblIndx, "centroid");
	pushVector3(L, n);
	lua_setfield(L, tblIndx, "n");
	pushVector3(L, t1);
	lua_setfield(L, tblIndx, "t1");
	pushVector3(L, t2);
	lua_setfield(L, tblIndx, "t2");
	lua_pushnumber(L, area);
	lua_setfield(L, tblIndx, "area");
	return 1;
    } else {
	string errMsg = "quadProperties enpected a table with names corners.";
	luaL_error(L, errMsg.toStringz);
	return 0;
    }
} // end quadProperties()

/**
 * Returns a table containing the named properties.
 * Expects a table of names corners.
 */
extern(C) int hexCellProperties(lua_State* L)
{
    if ( lua_istable(L, 1) ) {
	lua_getfield(L, 1, "p0");
	Vector3 p0 = *checkVector3(L, -1);
	lua_pop(L, 1);
	lua_getfield(L, 1, "p1");
	Vector3 p1 = *checkVector3(L, -1);
	lua_pop(L, 1);
	lua_getfield(L, 1, "p2");
	Vector3 p2 = *checkVector3(L, -1);
	lua_pop(L, 1);
	lua_getfield(L, 1, "p3");
        Vector3 p3 = *checkVector3(L, -1);
	lua_pop(L, 1);
	lua_getfield(L, 1, "p4");
	Vector3 p4 = *checkVector3(L, -1);
	lua_pop(L, 1);
	lua_getfield(L, 1, "p5");
	Vector3 p5 = *checkVector3(L, -1);
	lua_pop(L, 1);
	lua_getfield(L, 1, "p6");
	Vector3 p6 = *checkVector3(L, -1);
	lua_pop(L, 1);
	lua_getfield(L, 1, "p7");
        Vector3 p7 = *checkVector3(L, -1);
	lua_pop(L, 1);
	lua_settop(L, 0); // clear stack

	Vector3 centroid;
	double volume;
	double iLen, jLen, kLen;
	hex_cell_properties(p0, p1, p2, p3, p4, p5, p6, p7, 
			    centroid, volume, iLen, jLen, kLen);

	lua_newtable(L); // anonymous table { }
	auto tblIndx = lua_gettop(L);
	pushVector3(L, centroid);
	lua_setfield(L, tblIndx, "centroid");
	lua_pushnumber(L, volume);
	lua_setfield(L, tblIndx, "volume");
	lua_pushnumber(L, iLen);
	lua_setfield(L, tblIndx, "iLen");
	lua_pushnumber(L, jLen);
	lua_setfield(L, tblIndx, "jLen");
	lua_pushnumber(L, kLen);
	lua_setfield(L, tblIndx, "kLen");
	return 1;
    } else {
	string errMsg = "hexCellProperties enpected a table with names corners.";
	luaL_error(L, errMsg.toStringz);
	return 0;
    }
} // end hexCellProperties()

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
    lua_pushcfunction(L, &rotateAboutZAxisVector3);
    lua_setfield(L, -2, "rotateAboutZAxis");
    lua_pushcfunction(L, &mirrorImageVector3);
    lua_setfield(L, -2, "mirrorImage");

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
    lua_pushcfunction(L, &quadProperties);
    lua_setglobal(L, "quadProperties");
    lua_pushcfunction(L, &hexCellProperties);
    lua_setglobal(L, "hexCellProperties");
} // end registerVector3()
