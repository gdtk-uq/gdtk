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

module geom.luawrap.luageom;

import std.stdio;
import std.string;
import std.conv;
import ntypes.complex;
import nm.number;
import util.lua;
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
 * Returns a Vector3 constructed from the item on the Lua stack.
 * This item may be a true Vector3 pointer (to userdata) or it
 * may be a table with suitable fields, x, y, z.
 */
Vector3 toVector3(lua_State *L, int index)
{
    // First, see if it is a udata blob of Vector3 type.
    if (lua_isuserdata(L, index)) {
        return *checkVector3(L, index);
    }
    // If we arrive here, see if we have a table with x,y,z fields.
    Vector3 vec = Vector3(0.0, 0.0, 0.0);
    if (!lua_istable(L, index)) {
        luaL_error(L, "Did not get a Vector3 udata object nor a table.");
        return vec;
    }
    // Have table, now look for named fields containing coordinate values.
    lua_getfield(L, index, "x");
    if (lua_isnumber(L, -1)) { vec.x = lua_tonumber(L, -1); }
    lua_pop(L, 1);
    lua_getfield(L, index, "y");
    if (lua_isnumber(L, -1)) { vec.y = lua_tonumber(L, -1); }
    lua_pop(L, 1);
    lua_getfield(L, index, "z");
    if (lua_isnumber(L, -1)) { vec.z = lua_tonumber(L, -1); }
    lua_pop(L, 1);
    return vec;
} // end toVector3

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
 * assert(c.x == 1.0); assert(c.y == 3.0); assert(c.z == 2.0)
 * ----------------------
 * When a single argument is given, it may be another Vector3 object.
 * For the table of coordinates, missing values are assumed to be 0.0.
 */
extern(C) int newVector3(lua_State *L)
{
    auto vec = Vector3(0.0, 0.0, 0.0);
    // On being called, we are expecting at least two arguments.
    // The first will be a table being the prototype object "this"
    // and the next will may a table with the arguments with which
    // to set the data for our newly constructed point
    // or it may be another Vector3 object.
    int narg = lua_gettop(L);
    if (narg == 2 && lua_istable(L, 1)) {
        lua_remove(L, 1); // remove first argument "this".
        // Could be a table or a single Vector3 object
        if (lua_istable(L, 1)) {
            // If it has a length >= 1, then it's been populated array style.
            // This style of setting beats any fields that are present.
            if (lua_objlen(L, 1) >= 1) {
                // A single item may be a Vector3 object already.
                lua_rawgeti(L, 1, 1);
                vec = *checkVector3(L, -1);
                lua_pop(L, 1);
            } else {
                // Look for named fields containing coordinate values.
                lua_getfield(L, 1, "x");
                if ( lua_isnumber(L, -1) ) vec.x = lua_tonumber(L, -1);
                lua_pop(L, 1);
                lua_getfield(L, 1, "y");
                if ( lua_isnumber(L, -1) ) vec.y = lua_tonumber(L, -1);
                lua_pop(L, 1);
                lua_getfield(L, 1, "z");
                if ( lua_isnumber(L, -1) ) vec.z = lua_tonumber(L, -1);
                lua_pop(L, 1);
            }
        } else {
            vec = *checkVector3(L, -1);
        }
    } else {
        // We did not get what we expected as arguments.
        string errMsg = "Expected Vector3:new{x=number, y=number}; ";
        errMsg ~= "maybe you tried Vector3.new{x=number, y=number}.";
        luaL_error(L, errMsg.toStringz);
    }
    // Regardless of how we filled in vec, we are now ready
    // to place our new Vector3 on the stack as userdata.
    return pushVector3(L, vec);
} // end newVector3()


extern(C) int indexVector3(lua_State* L)
{
    auto a = checkVector3(L, 1);
    string key = to!string(luaL_checkstring(L, 2));
    // If we find an "x", "y" or "z", get value and return
    if ( key == "x" ) {
        lua_pushnumber(L, a.x);
        return 1;
    }
    if ( key == "y" ) {
        lua_pushnumber(L, a.y);
        return 1;
    }
    if ( key == "z") {
        lua_pushnumber(L, a.z);
        return 1;
    }
    // else forward through to the metatable
    lua_getmetatable(L, 1);
    lua_getfield(L, -1, toStringz(key));
    return 1;
}

extern(C) int newindexVector3(lua_State* L)
{
    auto a = checkVector3(L, 1);
    string key = to!string(luaL_checkstring(L, 2));
    double val = luaL_checknumber(L, 3);
    // If we find an "x", "y" or "z", get value and return
    if ( key == "x" ) {
        a.x = val;
        return 0;
    }
    if ( key == "y" ) {
        a.y = val;
        return 0;
    }
    if ( key == "z") {
        a.z = val;
        return 0;
    }
    // else just ignore the setter
    return 0;
}

/*-------- exposed Vector3 methods ------------ */

/**
 * This provides the unary minus operator for Lua.
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
 * Transforms a Vector3 object into a local frame defined by three unit vectors n,t1,t2.
 * Exposes geom.Vector3.transform_to_local_frame(n, t1, t2)
 */
extern(C) int transformToLocalFrameVector3(lua_State* L)
{
    auto a = checkVector3(L, 1);
    auto n = checkVector3(L, 2);
    auto t1 = checkVector3(L, 3);
    auto t2 = checkVector3(L, 4);
    a.transform_to_local_frame(*n, *t1, *t2);
    return 0;
}

/**
 * Transforms a Vector3 object from local frame n,t1,t2 to global frame.
 * Exposes geom.Vector3.transform_to_global_frame(n, t1, t2)
 */
extern(C) int transformToGlobalFrameVector3(lua_State* L)
{
    auto a = checkVector3(L, 1);
    auto n = checkVector3(L, 2);
    auto t1 = checkVector3(L, 3);
    auto t2 = checkVector3(L, 4);
    a.transform_to_global_frame(*n, *t1, *t2);
    return 0;
}

/**
 * Returns a table containing the named properties.
 * Expects a table of names corners.
 */
extern(C) int quadProperties(lua_State* L)
{
    if ( lua_istable(L, 1) ) {
        if (!checkAllowedNames(L, 1, ["p0", "p1", "p2", "p3"])) {
            string errMsg = "Error in call to quadProperties{}. Invalid name in table.";
            luaL_error(L, errMsg.toStringz);
        }
        lua_getfield(L, 1, "p0");
        Vector3 p0 = toVector3(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, 1, "p1");
        Vector3 p1 = toVector3(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, 1, "p2");
        Vector3 p2 = toVector3(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, 1, "p3");
        Vector3 p3 = toVector3(L, -1);
        lua_pop(L, 1);
        lua_settop(L, 0); // clear stack

        Vector3 centroid;
        Vector3 n;
        Vector3 t1;
        Vector3 t2;
        number area;
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
        string errMsg = "quadProperties expected a table with named corner points.";
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
        if (!checkAllowedNames(L, 1, ["p0", "p1", "p2", "p3", "p4", "p5", "p6", "p7"])) {
            string errMsg = "Error in call to hexCellProperties{}. Invalid name in table.";
            luaL_error(L, errMsg.toStringz);
        }
        lua_getfield(L, 1, "p0");
        Vector3 p0 = toVector3(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, 1, "p1");
        Vector3 p1 = toVector3(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, 1, "p2");
        Vector3 p2 = toVector3(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, 1, "p3");
        Vector3 p3 = toVector3(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, 1, "p4");
        Vector3 p4 = toVector3(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, 1, "p5");
        Vector3 p5 = toVector3(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, 1, "p6");
        Vector3 p6 = toVector3(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, 1, "p7");
        Vector3 p7 = toVector3(L, -1);
        lua_pop(L, 1);
        lua_settop(L, 0); // clear stack

        Vector3 centroid;
        number volume;
        number iLen, jLen, kLen;
        hex_cell_properties(p0, p1, p2, p3, p4, p5, p6, p7,
                            true, centroid, volume, iLen, jLen, kLen);

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
        string errMsg = "hexCellProperties expected a table with named corner points.";
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
    //    lua_pushvalue(L, -1); // duplicates the current metatable
    //    lua_setfield(L, -2, "__index");
    lua_pushcfunction(L, &indexVector3);
    lua_setfield(L, -2, "__index");
    lua_pushcfunction(L, &newindexVector3);
    lua_setfield(L, -2, "__newindex");

    /* Register methods for use. */
    lua_pushcfunction(L, &newVector3);
    lua_setfield(L, -2, "new");
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
    lua_pushcfunction(L, &transformToLocalFrameVector3);
    lua_setfield(L, -2, "transformToLocalFrame");
    lua_pushcfunction(L, &transformToGlobalFrameVector3);
    lua_setfield(L, -2, "transformToGlobalFrame");

    lua_setglobal(L, Vector3MT.toStringz);

    /* Also attempt to put "add" and friends in the global namespace. */
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
