/**
 * A Lua interface for the D gpath module.
 *
 * Authors: Rowan G. and Peter J.
 * Date: 2015-02-22 First code
 *       2015-04-22 greatly expanded with Arc, Bezier, Polyline, LuaFnPath
 */

module luagpath;

// We cheat to get the C Lua headers by using LuaD.
import util.lua;
import std.stdio;
import std.string;
import std.conv;
import util.lua_service;
import geom;
import gpath;
import luageom;

immutable string LineMT = "Line"; // Name of Line metatable
immutable string ArcMT = "Arc";
immutable string Arc3MT = "Arc3";
immutable string BezierMT = "Bezier";
immutable string PolylineMT = "Polyline";
immutable string SplineMT = "Spline";
immutable string LuaFnPathMT = "LuaFnPath";
immutable string ArcLengthParameterizedPathMT = "ArcLengthParameterizedPath";
immutable string SubRangedPathMT = "SubRangedPath";
immutable string ReversedPathMT = "ReversedPath";

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
    if ( isObjType(L, index, BezierMT) ) {
	return checkObj!(Bezier, BezierMT)(L, index);
    }
    if ( isObjType(L, index, PolylineMT) ) {
	return checkObj!(Polyline, PolylineMT)(L, index);
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

/* ----------------- Path-specific functions --------------- */

/**
 * The Lua constructor for a Line.
 *
 * Example construction in Lua:
 * ---------------------------------
 * a = Vector3:new{}
 * b = Vector3:new{1, 1}
 * ab = Line:new{a, b}
 * ---------------------------------
 */
extern(C) int newLine(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = `Error in call to Line:new{}.;
A table containing arguments is expected, but no table was found.`;
	luaL_error(L, errMsg.toStringz);
    }
    // Expect Vector3 at position 1.
    lua_rawgeti(L, 1, 1);
    auto a = checkVector3(L, -1);
    if ( a is null ) {
	string errMsg = `Error in call to Line:new{}.
A Vector3 object is expected as the first argument. No valid Vector3 was found.`;
	luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1);
    // Expect Vector3 at position 2.
    lua_rawgeti(L, 1, 2);
    auto b = checkVector3(L, -1);
    if ( b is null ) {
	string errMsg = `Error in call to Line:new{}.
A Vector3 object is expected as the second argument. No valid Vector3 was found.`;
	luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1);
    auto ab = new Line(*a, *b);
    pathStore ~= pushObj!(Line, LineMT)(L, ab);
    return 1;
} // end newLine()


/**
 * The Lua constructor for an Arc.
 *
 * Example construction in Lua:
 * ---------------------------------
 * a = Vector3:new{1, 0}
 * b = Vector3:new{0, 1}
 * c = Vector3:new{0, 0}
 * abc = Arc:new{a, b, c}
 * ---------------------------------
 */
extern(C) int newArc(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = `Error in call to Arc:new{}.;
A table containing arguments is expected, but no table was found.`;
	luaL_error(L, errMsg.toStringz);
    }
    // Expect Vector3 at position 1.
    lua_rawgeti(L, 1, 1);
    auto a = checkVector3(L, -1);
    if ( a is null ) {
	string errMsg = `Error in call to Arc:new{}.
A Vector3 object is expected as the first argument. No valid Vector3 was found.`;
	luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1);
    // Expect Vector3 at position 2.
    lua_rawgeti(L, 1, 2);
    auto b = checkVector3(L, -1);
    if ( b is null ) {
	string errMsg = `Error in call to Arc:new{}.
A Vector3 object is expected as the second argument. No valid Vector3 was found.`;
	luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1);
    // Expect Vector3 at position 3.
    lua_rawgeti(L, 1, 3);
    auto c = checkVector3(L, -1);
    if ( c is null ) {
	string errMsg = `Error in call to Arc:new{}.
A Vector3 object is expected as the third argument. No valid Vector3 was found.`;
	luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1);
    auto abc = new Arc(*a, *b, *c);
    pathStore ~= pushObj!(Arc, ArcMT)(L, abc);
    return 1;
} // end newArc()


/**
 * The Lua constructor for an Arc3.
 *
 * Example construction in Lua:
 * ---------------------------------
 * a = Vector3:new{1, 0}
 * m = Vector3:new{0.707107, 0.707107}
 * b = Vector3:new{0, 1}
 * amb = Arc:new{a, m, b}
 * ---------------------------------
 */
extern(C) int newArc3(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = `Error in call to Arc3:new{}.;
A table containing arguments is expected, but no table was found.`;
	luaL_error(L, errMsg.toStringz);
    }
    // Expect Vector3 at position 1.
    lua_rawgeti(L, 1, 1);
    auto a = checkVector3(L, -1);
    if ( a is null ) {
	string errMsg = `Error in call to Arc3:new{}.
A Vector3 object is expected as the first argument. No valid Vector3 was found.`;
	luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1);
    // Expect Vector3 at position 2.
    lua_rawgeti(L, 1, 2);
    auto m = checkVector3(L, -1);
    if ( m is null ) {
	string errMsg = `Error in call to Arc3:new{}.
A Vector3 object is expected as the second argument. No valid Vector3 was found.`;
	luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1);
    // Expect Vector3 at position 3.
    lua_rawgeti(L, 1, 3);
    auto b = checkVector3(L, -1);
    if ( b is null ) {
	string errMsg = `Error in call to Arc3:new{}.
A Vector3 object is expected as the third argument. No valid Vector3 was found.`;
	luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1);
    auto amb = new Arc3(*a, *m, *b);
    pathStore ~= pushObj!(Arc3, Arc3MT)(L, amb);
    return 1;
} // end newArc3()


/**
 * The Lua constructor for a Bezier.
 *
 * Example construction in Lua:
 * ---------------------------------
 * p0 = Vector3:new{1, 0}
 * p1 = Vector3:new{0.7071, 0.7071}
 * p2 = Vector3:new{0, 1}
 * -- For an arbitrary number of points in the table.
 * bez = Bezier:new{points={p0, p1, p2}}
 * ---------------------------------
 */
extern(C) int newBezier(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = `Error in call to Bezier:new{}.;
A table containing arguments is expected, but no table was found.`;
	luaL_error(L, errMsg.toStringz);
    }
    lua_getfield(L, 1, "points".toStringz());
    if ( !lua_istable(L, -1) ) {
	string errMsg = `Error in call to Bezier:new{}.;
A table containing Vector3 points is expected, but no table was found.`;
	luaL_error(L, errMsg.toStringz);
    }
    // Expect Vector3 objects at array positions within that table.
    Vector3[] B;
    int position = 1;
    while ( true ) {
	lua_rawgeti(L, -1, position);
	if ( lua_isnil(L, -1) ) { lua_pop(L, 1); break; }
	auto a = checkVector3(L, -1);
	lua_pop(L, 1);
	if ( a is null ) break;
	B ~= *a;
	++position;
    }
    lua_pop(L, 1); // dispose of points table
    if ( B.length == 0 ) {
	string errMsg = `Error in call to Bezier:new{}. No valid Vector3 objects found.`;
	luaL_error(L, errMsg.toStringz());
    }
    auto bez = new Bezier(B);
    pathStore ~= pushObj!(Bezier, BezierMT)(L, bez);
    return 1;
} // end newBezier()


/**
 * The Lua constructor for a Polyline.
 *
 * Example construction in Lua:
 * ---------------------------------
 * p0 = Vector3:new{1, 0}
 * p1 = Vector3:new{0.7071, 0.7071}
 * p2 = Vector3:new{0, 1}
 * -- A couple of paths to combine.
 * line1 = Line:new{p0, p1}
 * line2 = Line:new{p1, p2}
 * -- An arbitrary number of Path objects in the table.
 * poly = Polyline:new{line1, line2}
 * ---------------------------------
 */
extern(C) int newPolyline(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = `Error in call to Polyline:new{}.;
A table containing arguments is expected, but no table was found.`;
	luaL_error(L, errMsg.toStringz);
    }
    // Expect Path objects at array positions.
    Path[] segments;
    int position = 1;
    while ( true ) {
	lua_rawgeti(L, 1, position);
	if ( lua_isnil(L, -1) ) { lua_pop(L, 1); break; }
	auto seg = checkPath(L, -1);
	lua_pop(L, 1);
	if ( seg is null ) break;
	segments ~= seg;
	++position;
    }
    if ( segments.length == 0 ) {
	string errMsg = `Error in call to Polyline:new{}. No valid Path objects found.`;
	luaL_error(L, errMsg.toStringz());
    }
    auto poly = new Polyline(segments);
    pathStore ~= pushObj!(Polyline, PolylineMT)(L, poly);
    return 1;
} // end newPolyline()


/**
 * The Lua constructor for a Spline (Polyline).
 *
 * Example construction in Lua:
 * ---------------------------------
 * p0 = Vector3:new{ 0, -1, 0}
 * p1 = Vector3:new{-1,  0, 0}
 * p2 = Vector3:new{ 0,  1, 0}
 * p3 = Vector3:new{ 1,  0, 0}
 * p4 = Vector3:new{ 0, -1, 0}
 * -- For an arbitrary number of points in the table.
 * bez = Spline:new{points={p0, p1, p2, p3, p4}}
 * ---------------------------------
 */
extern(C) int newSpline(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = `Error in call to Spline:new{}.;
A table containing arguments is expected, but no table was found.`;
	luaL_error(L, errMsg.toStringz);
    }
    lua_getfield(L, 1, "points".toStringz());
    if ( !lua_istable(L, -1) ) {
	string errMsg = `Error in call to Spline:new{}.;
A table containing Vector3 points is expected, but no table was found.`;
	luaL_error(L, errMsg.toStringz);
    }
    // Expect Vector3 objects at array positions.
    Vector3[] B;
    int position = 1;
    while ( true ) {
	lua_rawgeti(L, -1, position);
	if ( lua_isnil(L, -1) ) { lua_pop(L, 1); break; }
	auto a = checkVector3(L, -1);
	lua_pop(L, 1);
	if ( a is null ) break;
	B ~= *a;
	++position;
    }
    lua_pop(L, 1); // dispose of points table
    if ( B.length == 0 ) {
	string errMsg = `Error in call to Spline:new{}. No valid Vector3 objects found.`;
	luaL_error(L, errMsg.toStringz());
    }
    auto spline = new Polyline(B);
    pathStore ~= pushObj!(Polyline, PolylineMT)(L, spline);
    return 1;
} // end newSpline()


/**
 * LuaFnPath class and it's Lua constructor.
 *
 * This is hangs onto a Lua call-back function that is invoked from the D domain.
 *
 * Example:
 * function myLuaFunction(t)
 *    -- Straight line from 0,0,0 to 1.0,2.0,3.0
 *    return {t, 2*t, 3*t}
 * end
 * myPath = LuaFnPath:new{"myLuaFunction"}
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
	    string errMsg = `Error in call to LuaFnPath:opCall().;
A table containing arguments is expected, but no table was found.`;
	    luaL_error(cast(lua_State*)L, errMsg.toStringz);
	}
	// Expect a number at position 1 for x.
	lua_rawgeti(cast(lua_State*)L, -1, 1);
	if ( !lua_isnumber(cast(lua_State*)L, -1) ) {
	    string errMsg = `Error in call to LuaFnPath:opCall().;
A number was expected in position 1, for x.`;
	    luaL_error(cast(lua_State*)L, errMsg.toStringz);
	}
	double x = to!double(lua_tonumber(cast(lua_State*)L, -1));
	lua_pop(cast(lua_State*)L, 1);
	// Expect a number at position 2 for y.
	lua_rawgeti(cast(lua_State*)L, -1, 2);
	if ( !lua_isnumber(cast(lua_State*)L, -1) ) {
	    string errMsg = `Error in call to LuaFnPath:opCall().;
A number was expected in position 2, for y.`;
	    luaL_error(cast(lua_State*)L, errMsg.toStringz);
	}
	double y = to!double(lua_tonumber(cast(lua_State*)L, -1));
	lua_pop(cast(lua_State*)L, 1);
	// Expect a number at position 3 for z.
	double z = 0.0; // default value
	lua_rawgeti(cast(lua_State*)L, -1, 3);
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
} // end class LuaFnPath

extern(C) int newLuaFnPath(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = `Error in call to LuaFnPath:new{}.;
A table containing arguments is expected, but no table was found.`;
	luaL_error(L, errMsg.toStringz);
    }
    // Expect function name at array position 1.
    string fnName = "";
    lua_rawgeti(L, 1, 1);
    if ( lua_isstring(L, -1) ) {
	fnName ~= to!string(lua_tostring(L, -1));
    }
    lua_pop(L, 1);
    if ( fnName == "" ) {
	string errMsg = `Error in call to LuaFnPath:new{}. No function name found.`;
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
 * p0 = Vector3:new{1, 0}
 * p1 = Vector3:new{0.7071, 0.7071}
 * p2 = Vector3:new{0, 1}
 * original_path = Bezier:new{points={p0, p1, p2}}
 * alp_path = ArcLengthParameterizedPath:new{original_path}
 * ---------------------------------
 */
extern(C) int newArcLengthParameterizedPath(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = `Error in call to ArcLengthParameterizedPath:new{}.;
A table containing arguments is expected, but no table was found.`;
	luaL_error(L, errMsg.toStringz);
    }
    // Expect a Path object at the first array position in the table.
    lua_rawgeti(L, 1, 1);
    if ( lua_isnil(L, -1) ) {
	string errMsg = `Error in call to ArcLengthParameterizedPath:new{}. No table entry found.`;
	luaL_error(L, errMsg.toStringz());
    }
    auto original_path = checkPath(L, -1);
    lua_pop(L, 1);
    if ( original_path is null ) {
	string errMsg = `Error in call to ArcLengthParameterizedPath:new{}. No valid Path object found.`;
	luaL_error(L, errMsg.toStringz());
    }
    auto alp_path = new ArcLengthParameterizedPath(original_path);
    pathStore ~= pushObj!(ArcLengthParameterizedPath,
			  ArcLengthParameterizedPathMT)(L, alp_path);
    return 1;
} // end newArcLengthParameterizedPath()


/**
 * The Lua constructor for an SubRangedPath.
 *
 * Example construction in Lua:
 * ---------------------------------
 * p0 = Vector3:new{1, 0}
 * p1 = Vector3:new{0.7071, 0.7071}
 * p2 = Vector3:new{0, 1}
 * original_path = Bezier:new{points={p0, p1, p2}}
 * sr_path = SubRangedPath:new{original_path, t0=0.1, t1=0.9}
 * ---------------------------------
 */
extern(C) int newSubRangedPath(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = `Error in call to SubRangedPath:new{}.;
A table containing arguments is expected, but no table was found.`;
	luaL_error(L, errMsg.toStringz);
    }
    // Expect a Path object at the first array position in the table.
    lua_rawgeti(L, 1, 1);
    if ( lua_isnil(L, -1) ) {
	string errMsg = `Error in call to SubRangedPath:new{}. No table entry found.`;
	luaL_error(L, errMsg.toStringz());
    }
    auto original_path = checkPath(L, -1);
    lua_pop(L, 1);
    if ( original_path is null ) {
	string errMsg = `Error in call to SubRangedPath:new{}. No valid Path object found.`;
	luaL_error(L, errMsg.toStringz());
    }
    string errMsgTmplt = `Error in call to SubRangedPath:new{}.
A valid value for '%s' was not found in list of arguments.
The value, if present, should be a number.`;
    double t0 = getNumberFromTable(L, 1, "t0", false, 0.0, true, format(errMsgTmplt, "t0"));
    double t1 = getNumberFromTable(L, 1, "t1", false, 1.0, true, format(errMsgTmplt, "t1"));
    auto alp_path = new SubRangedPath(original_path, t0, t1);
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
 * p0 = Vector3:new{1, 0}
 * p1 = Vector3:new{0.7071, 0.7071}
 * p2 = Vector3:new{0, 1}
 * original_path = Bezier:new{points={p0, p1, p2}}
 * r_path = ReversedPath:new{original_path}
 * ---------------------------------
 */
extern(C) int newReversedPath(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = `Error in call to ReversedPath:new{}.;
A table containing arguments is expected, but no table was found.`;
	luaL_error(L, errMsg.toStringz);
    }
    // Expect a Path object at the first array position in the table.
    lua_rawgeti(L, 1, 1);
    if ( lua_isnil(L, -1) ) {
	string errMsg = `Error in call to ReversedPath:new{}. No table entry found.`;
	luaL_error(L, errMsg.toStringz());
    }
    auto original_path = checkPath(L, -1);
    lua_pop(L, 1);
    if ( original_path is null ) {
	string errMsg = `Error in call to ReversedPath:new{}. No valid Path object found.`;
	luaL_error(L, errMsg.toStringz());
    }
    auto alp_path = new ReversedPath(original_path);
    pathStore ~= pushObj!(ReversedPath, ReversedPathMT)(L, alp_path);
    return 1;
} // end newReversedPath()


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

    lua_setglobal(L, Arc3MT.toStringz);

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

    lua_setglobal(L, BezierMT.toStringz);

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

    lua_setglobal(L, SplineMT.toStringz);

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

    lua_setglobal(L, ReversedPathMT.toStringz);
} // end registerPaths()
    






