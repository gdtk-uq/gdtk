/**
 * luasketch.d  A Lua interface to the Sketch module.
 *
 * Authors: Peter J. and Rowan G.
 * Version: 2016-08-11 First code (adapted from our other Lua interface modules)
 * 
 */

module luasketch;

import util.lua;
import std.stdio;
import std.string;
import std.conv;
import util.lua_service;
import geom;
import gpath;
import surface;
import sketch;
import luageom;
import luasurface;

immutable string SketchMT = "Sketch"; // Name of Sketch metatable

// A place to hang on to references to objects that are pushed into the Lua domain.
// We don't want the D garbage collector to prematurely dispose of said objects.
static const(Sketch)[] sketchStore;

Sketch checkSketch(lua_State* L, int index) {
    if ( isObjType(L, index, SketchMT) ) {
	return checkObj!(Sketch, SketchMT)(L, index);
    }
    // if check else fails
    return null;
}

/**
 * Lua constructor for a Sketch.
 *
 * Examples:
 * s = Sketch:new{renderer="svg", projection="xyortho"}
 * s = Sketch:new{} -- as above
 *
 */
extern(C) int newSketch(lua_State* L)
{
    lua_remove(L, 1); // remove first argument "this"
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = "Error in call to Sketch:new{}.; " ~
	    "A table containing arguments is expected, but no table was found.";
	luaL_error(L, errMsg.toStringz);
    }
    //
    string renderer_name = "";
    lua_getfield(L, 1, "renderer");
    if ( lua_isstring(L, -1) ) { renderer_name = to!string(lua_tostring(L, -1)); }
    lua_pop(L, 1);
    if (renderer_name == "") { renderer_name = "svg"; }
    //
    string projection_name = "";
    lua_getfield(L, 1, "projection");
    if ( lua_isstring(L, -1) ) { projection_name = to!string(lua_tostring(L, -1)); }
    lua_pop(L, 1);
    if (projection_name == "") { projection_name = "xyortho"; }
    //
    auto my_sketch = new Sketch(renderer_name, projection_name);
    sketchStore ~= pushObj!(Sketch, SketchMT)(L, my_sketch);
    return 1;
} // end newSketch()

extern(C) int startSketch(lua_State* L)
{
    auto my_sketch = checkObj!(Sketch, SketchMT)(L, 1); // first argument "this"
    lua_remove(L, 1); // remove first argument "this"
    // Now, I'm expecting a table of arguments as the only item on the stack.
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = "Error in call to Sketch:start{}.; " ~
	    "A table containing arguments is expected, but no table was found.";
	luaL_error(L, errMsg.toStringz);
    }
    //
    string file_name = "";
    lua_getfield(L, 1, "file_name");
    if ( lua_isstring(L, -1) ) {
	file_name = to!string(lua_tostring(L, -1));
    } else {
	string errMsg = "Error in call to Sketch:start{}. Expected a string for file_name.";
	luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1); // dispose of file_name item
    //
    my_sketch.start(file_name);
    return 0;
} // end startSketch()

extern(C) int finishSketch(lua_State* L)
{
    auto my_sketch = checkObj!(Sketch, SketchMT)(L, 1); // first argument "this"
    lua_remove(L, 1); // remove first argument "this"
    // Now, I could be expecting a table of arguments, but I don't care.
    my_sketch.finish();
    return 0;
} // end endSketch()

extern(C) int setSketch(lua_State* L)
{
    auto my_sketch = checkObj!(Sketch, SketchMT)(L, 1); // first argument "this"
    lua_remove(L, 1); // remove first argument "this"
    // Now, I'm expecting a table of arguments as the only item on the stack.
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = "Error in call to Sketch:set{}.; " ~
	    "A table containing arguments is expected, but no table was found.";
	luaL_error(L, errMsg.toStringz);
    }
    //
    lua_getfield(L, 1, "line_colour");
    if ( !lua_isnil(L, -1) ) {
	if ( lua_isstring(L, -1) ) {
	    string lineColour = to!string(lua_tostring(L, -1));
	    my_sketch.setLineColour(lineColour);
	} else {
	    string errMsg = "Error in call to Sketch:set{}. Expected a string for line_colour.";
	    luaL_error(L, errMsg.toStringz());
	}
    }
    lua_pop(L, 1); // dispose of line_colour item
    //
    lua_getfield(L, 1, "fill_colour");
    if ( !lua_isnil(L, -1) ) {
	if ( lua_isstring(L, -1) ) {
	    string fillColour = to!string(lua_tostring(L, -1));
	    my_sketch.setFillColour(fillColour);
	} else {
	    string errMsg = "Error in call to Sketch:set{}. Expected a string for fill_colour.";
	    luaL_error(L, errMsg.toStringz());
	}
    }
    lua_pop(L, 1); // dispose of fill_colour item
    //
    lua_getfield(L, 1, "line_width");
    if ( !lua_isnil(L, -1) ) {
	if ( lua_isnumber(L, -1) ) {
	    auto w = to!double(lua_tonumber(L, -1));
	    my_sketch.setLineWidth(w);
	} else {
	    string errMsg = "Error in call to Sketch:set{}. Expected a number for line_width.";
	    luaL_error(L, errMsg.toStringz());
	}
    }
    lua_pop(L, 1); // dispose of line_width item
    //
    lua_getfield(L, 1, "dash_array");
    if ( !lua_isnil(L, -1) ) {
	double dashLength = 2.0;
	double gapLength = 2.0;
	if ( lua_istable(L, -1) ) {
	    lua_rawgeti(L, -1, 1);
	    if (lua_isnumber(L, -1)) { dashLength = to!double(lua_tonumber(L, -1)); }
	    lua_pop(L, 1);
	    lua_rawgeti(L, -1, 2);
	    if (lua_isnumber(L, -1)) { gapLength = to!double(lua_tonumber(L, -1)); }
	    lua_pop(L, 1);
	    my_sketch.setDashArray(dashLength, gapLength);
	} else {
	    string errMsg = "Error in call to Sketch:set{}. Expected a table for dash_array.";
	    luaL_error(L, errMsg.toStringz());
	}
    }
    lua_pop(L, 1); // dispose of dash_array item
    //
    lua_getfield(L, 1, "canvas");
    if ( !lua_isnil(L, -1) ) {
	double x0 = 0.0; double y0 = 120.0;
	double x1 = 0.0; double y1 = 120.0;
	if ( lua_istable(L, -1) ) {
	    lua_rawgeti(L, -1, 1);
	    if (lua_isnumber(L, -1)) { x0 = to!double(lua_tonumber(L, -1)); }
	    lua_pop(L, 1);
	    lua_rawgeti(L, -1, 2);
	    if (lua_isnumber(L, -1)) { y0 = to!double(lua_tonumber(L, -1)); }
	    lua_pop(L, 1);
	    lua_rawgeti(L, -1, 3);
	    if (lua_isnumber(L, -1)) { x1 = to!double(lua_tonumber(L, -1)); }
	    lua_pop(L, 1);
	    lua_rawgeti(L, -1, 4);
	    if (lua_isnumber(L, -1)) { y1 = to!double(lua_tonumber(L, -1)); }
	    lua_pop(L, 1);
	    my_sketch.canvas.set(x0, y0, x1, y1);
	} else {
	    string errMsg = "Error in call to Sketch:set{}. Expected a table for canvas.";
	    luaL_error(L, errMsg.toStringz());
	}
    }
    lua_pop(L, 1); // dispose of viewport item
    //
    lua_getfield(L, 1, "viewport");
    if ( !lua_isnil(L, -1) ) {
	double x0 = 0.0; double y0 = 120.0;
	double x1 = 0.0; double y1 = 120.0;
	if ( lua_istable(L, -1) ) {
	    lua_rawgeti(L, -1, 1);
	    if (lua_isnumber(L, -1)) { x0 = to!double(lua_tonumber(L, -1)); }
	    lua_pop(L, 1);
	    lua_rawgeti(L, -1, 2);
	    if (lua_isnumber(L, -1)) { y0 = to!double(lua_tonumber(L, -1)); }
	    lua_pop(L, 1);
	    lua_rawgeti(L, -1, 3);
	    if (lua_isnumber(L, -1)) { x1 = to!double(lua_tonumber(L, -1)); }
	    lua_pop(L, 1);
	    lua_rawgeti(L, -1, 4);
	    if (lua_isnumber(L, -1)) { y1 = to!double(lua_tonumber(L, -1)); }
	    lua_pop(L, 1);
	    my_sketch.viewport.set(x0, y0, x1, y1);
	} else {
	    string errMsg = "Error in call to Sketch:set{}. Expected a table for viewport.";
	    luaL_error(L, errMsg.toStringz());
	}
    }
    lua_pop(L, 1); // dispose of viewport item
    return 0;
} // end setSketch()

extern(C) int lineSketch(lua_State* L)
{
    auto my_sketch = checkObj!(Sketch, SketchMT)(L, 1); // first argument "this"
    lua_remove(L, 1); // remove first argument "this"
    // Now, I'm expecting a table of arguments as the only item on the stack.
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = "Error in call to Sketch:line{}.; " ~
	    "A table containing arguments is expected, but no table was found.";
	luaL_error(L, errMsg.toStringz);
    }
    //
    lua_getfield(L, 1, "p0");
    if ( lua_isnil(L, -1) ) {
	string errMsg = "Error in call to Sketch:line{}. No p0 entry found.";
	luaL_error(L, errMsg.toStringz());
    }
    auto p0 = checkVector3(L, -1);
    if ( p0 is null ) {
	string errMsg = "Error in call to Sketch:line{}. " ~
	    "A Vector3 object is expected as the p0 argument. No valid Vector3 was found.";
	luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1); // dispose of p0 item
    //
    lua_getfield(L, 1, "p1");
    if ( lua_isnil(L, -1) ) {
	string errMsg = "Error in call to Sketch:line{}. No p1 entry found.";
	luaL_error(L, errMsg.toStringz());
    }
    auto p1 = checkVector3(L, -1);
    if ( p1 is null ) {
	string errMsg = "Error in call to Sketch:line{}. " ~
	    "A Vector3 object is expected as the p1 argument. No valid Vector3 was found.";
	luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1); // dispose of p1 item
    //
    bool dashed = false;
    lua_getfield(L, 1, "dashed");
    if ( !lua_isnil(L, -1) ) {
	if ( lua_isboolean(L, -1) ) {
	    dashed = to!bool(lua_toboolean(L, -1));
	} else {
	    string errMsg = "Error in call to Sketch:line{}. Expected a bool for dashed.";
	    luaL_error(L, errMsg.toStringz());
	}
    }
    lua_pop(L, 1); // dispose of dashed item
    //
    my_sketch.line(*p0, *p1, dashed);
    return 0;
} // end lineSketch()

extern(C) int dotlabelSketch(lua_State* L)
{
    auto my_sketch = checkObj!(Sketch, SketchMT)(L, 1); // first argument "this"
    lua_remove(L, 1); // remove first argument "this"
    // Now, I'm expecting a table of arguments as the only item on the stack.
    int narg = lua_gettop(L);
    if ( narg == 0 || !lua_istable(L, 1) ) {
	string errMsg = "Error in call to Sketch:dotlabel{}.; " ~
	    "A table containing arguments is expected, but no table was found.";
	luaL_error(L, errMsg.toStringz);
    }
    //
    lua_getfield(L, 1, "point");
    if ( lua_isnil(L, -1) ) {
	string errMsg = "Error in call to Sketch:dotlabel{}. No point entry found.";
	luaL_error(L, errMsg.toStringz());
    }
    auto point = checkVector3(L, -1);
    if ( point is null ) {
	string errMsg = "Error in call to Sketch:dotlabel{}. " ~
	    "A Vector3 object is expected as the point argument. No valid Vector3 was found.";
	luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1); // dispose of p0 item
    //
    string label = "";
    lua_getfield(L, 1, "label");
    if ( lua_isstring(L, -1) ) {
	label = to!string(lua_tostring(L, -1));
    } else {
	string errMsg = "Error in call to Sketch:dotlabel{}. Expected a string for label.";
	luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1); // dispose of label item
    //
    string anchor = "middle";
    lua_getfield(L, 1, "anchor");
    if ( !lua_isnil(L, -1) ) {
	if ( lua_isstring(L, -1) ) {
	    anchor = to!string(lua_tostring(L, -1));
	} else {
	    string errMsg = "Error in call to Sketch:dotlabel{}. Expected a string for anchor.";
	    luaL_error(L, errMsg.toStringz());
	}
    }
    lua_pop(L, 1); // dispose of anchor item
    //
    double dotSize = 2.0;
    lua_getfield(L, 1, "dot_size");
    if ( !lua_isnil(L, -1) ) {
	if ( lua_isnumber(L, -1) ) {
	    dotSize = to!double(lua_tonumber(L, -1));
	} else {
	    string errMsg = "Error in call to Sketch:dotlabel{}. Expected a number for dot_size.";
	    luaL_error(L, errMsg.toStringz());
	}
    }
    lua_pop(L, 1); // dispose of dot_size item
    //
    my_sketch.dotlabel(*point, label, anchor, dotSize);
    return 0;
} // end dotlabelSketch()


void registerSketch(lua_State* L)
{
    // Register the Sketch object
    luaL_newmetatable(L, SketchMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    /* Register methods for use. */
    lua_pushcfunction(L, &newSketch);
    lua_setfield(L, -2, "new");
    lua_pushcfunction(L, &toStringObj!(Sketch, SketchMT));
    lua_setfield(L, -2, "__tostring");
    lua_pushcfunction(L, &startSketch);
    lua_setfield(L, -2, "start");
    lua_pushcfunction(L, &finishSketch);
    lua_setfield(L, -2, "finish");
    lua_pushcfunction(L, &setSketch);
    lua_setfield(L, -2, "set");
    lua_pushcfunction(L, &lineSketch);
    lua_setfield(L, -2, "line");
    lua_pushcfunction(L, &dotlabelSketch);
    lua_setfield(L, -2, "dotlabel");

    lua_setglobal(L, SketchMT.toStringz);
} // end registerSurfaces()
