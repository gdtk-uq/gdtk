/**
 * luasketch.d  A Lua interface to the Sketch module.
 *
 * Authors: Peter J. and Rowan G.
 * Version: 2016-08-11 First code (adapted from our other Lua interface modules)
 * 
 */

module geom.luawrap.luasketch;

import util.lua;
import std.stdio;
import std.string;
import std.conv;
import util.lua_service;
import geom;
import geom.luawrap.luageom;
import geom.luawrap.luagpath;
import geom.luawrap.luasurface;
import geom.luawrap.luavolume;

immutable string SketchMT = "Sketch"; // Name of Sketch metatable

// A place to hang on to references to objects that are pushed into the Lua domain.
// We don't want the D garbage collector to prematurely dispose of said objects.
static const(Sketch)[] sketchStore;

Sketch checkSketch(lua_State* L, int index) {
    if (isObjType(L, index, SketchMT)) {
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
    if (narg == 0 || !lua_istable(L, 1)) {
        string errMsg = "Error in call to Sketch:new{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    //
    string renderer_name = "";
    lua_getfield(L, 1, "renderer");
    if (lua_isstring(L, -1)) { renderer_name = to!string(lua_tostring(L, -1)); }
    lua_pop(L, 1);
    if (renderer_name == "") { renderer_name = "svg"; }
    //
    string projection_name = "";
    lua_getfield(L, 1, "projection");
    if (lua_isstring(L, -1)) { projection_name = to!string(lua_tostring(L, -1)); }
    lua_pop(L, 1);
    if (projection_name == "") { projection_name = "xyortho"; }
    //
    double x0 = 0.0; double y0 = 0.0;
    double x1 = 120.0; double y1 = 120.0;
    lua_getfield(L, 1, "canvas_mm"); // optional field
    if (!lua_isnil(L, -1)) {
        if (lua_istable(L, -1)) {
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
        }
    }
    lua_pop(L, 1); // dispose of canvas item
    //
    auto my_sketch = new Sketch(renderer_name, projection_name, [x0, y0, x1, y1]);
    sketchStore ~= pushObj!(Sketch, SketchMT)(L, my_sketch);
    return 1;
} // end newSketch()

extern(C) int startSketch(lua_State* L)
{
    auto my_sketch = checkObj!(Sketch, SketchMT)(L, 1); // first argument "this"
    lua_remove(L, 1); // remove first argument "this"
    // Now, I'm expecting a table of arguments as the only item on the stack.
    int narg = lua_gettop(L);
    if (narg == 0 || !lua_istable(L, 1)) {
        string errMsg = "Error in call to Sketch:start{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    //
    string file_name = "";
    lua_getfield(L, 1, "file_name");
    if (lua_isstring(L, -1)) {
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
    if (narg == 0 || !lua_istable(L, 1)) {
        string errMsg = "Error in call to Sketch:set{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    //
    lua_getfield(L, 1, "line_colour");
    if (!lua_isnil(L, -1)) {
        if (lua_isstring(L, -1)) {
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
    if (!lua_isnil(L, -1)) {
        if (lua_isstring(L, -1)) {
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
    if (!lua_isnil(L, -1)) {
        if (lua_isnumber(L, -1)) {
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
    if (!lua_isnil(L, -1)) {
        double dashLength = 2.0;
        double gapLength = 2.0;
        if (lua_istable(L, -1)) {
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
    /+
     We set the canvas corners above, in the Sketch constructor.
     They should remain fixed.

    lua_getfield(L, 1, "canvas");
    if (!lua_isnil(L, -1)) {
        double x0 = 0.0; double y0 = 0.0;
        double x1 = 120.0; double y1 = 120.0;
        if (lua_istable(L, -1)) {
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
    lua_pop(L, 1); // dispose of canvas item
    +/
    //
    lua_getfield(L, 1, "viewport");
    if (!lua_isnil(L, -1)) {
        double x0 = 0.0; double y0 = 120.0;
        double x1 = 0.0; double y1 = 120.0;
        if (lua_istable(L, -1)) {
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
    if (narg == 0 || !lua_istable(L, 1)) {
        string errMsg = "Error in call to Sketch:line{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    //
    lua_getfield(L, 1, "p0");
    if (lua_isnil(L, -1)) {
        string errMsg = "Error in call to Sketch:line{}. No p0 entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto p0 = toVector3(L, -1);
    lua_pop(L, 1); // dispose of p0 item
    //
    lua_getfield(L, 1, "p1");
    if (lua_isnil(L, -1)) {
        string errMsg = "Error in call to Sketch:line{}. No p1 entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto p1 = toVector3(L, -1);
    lua_pop(L, 1); // dispose of p1 item
    //
    bool dashed = false;
    lua_getfield(L, 1, "dashed");
    if (!lua_isnil(L, -1)) {
        if (lua_isboolean(L, -1)) {
            dashed = to!bool(lua_toboolean(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:line{}. Expected a bool for dashed.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of dashed item
    //
    my_sketch.line(p0, p1, dashed);
    return 0;
} // end lineSketch()

extern(C) int polylineSketch(lua_State* L)
{
    auto my_sketch = checkObj!(Sketch, SketchMT)(L, 1); // first argument "this"
    lua_remove(L, 1); // remove first argument "this"
    // Now, I'm expecting a table of arguments as the only item on the stack.
    int narg = lua_gettop(L);
    if (narg == 0 || !lua_istable(L, 1)) {
        string errMsg = "Error in call to Sketch:polyline{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    //
    lua_getfield(L, 1, "points".toStringz());
    if (lua_isnil(L, -1)) {
        string errMsg = "Error in call to Sketch:polyline{}. No points entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    if (!lua_istable(L, -1)) {
        string errMsg = "Error in call to Sketch:polyline{}.; " ~
            "A table containing Vector3 points is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect Vector3 objects at array positions within that table.
    Vector3[] points;
    int position = 1;
    while (true) {
        lua_rawgeti(L, -1, position);
        if (lua_isnil(L, -1)) { lua_pop(L, 1); break; }
        auto p = toVector3(L, -1);
        lua_pop(L, 1);
        points ~= p;
        ++position;
    }
    lua_pop(L, 1); // dispose of points table
    if (points.length == 0) {
        string errMsg = "Error in call to Sketch:polyline{}. No valid Vector3 objects found.";
        luaL_error(L, errMsg.toStringz());
    }
    //
    bool dashed = false;
    lua_getfield(L, 1, "dashed");
    if (!lua_isnil(L, -1)) {
        if (lua_isboolean(L, -1) ) {
            dashed = to!bool(lua_toboolean(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:polyline{}. Expected a bool for dashed.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of dashed item
    //
    my_sketch.polyline(points, dashed);
    return 0;
} // end polylineSketch()

extern(C) int arcSketch(lua_State *L)
{
    auto my_sketch = checkObj!(Sketch, SketchMT)(L, 1); // first argument "this"
    lua_remove(L, 1); // remove first argument "this"
    // Now, I'm expecting a table of arguments as the only item on the stack.
    int narg = lua_gettop(L);
    if (narg == 0 || !lua_istable(L, 1)) {
        string errMsg = "Error in call to Sketch:arc{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    //
    lua_getfield(L, 1, "p0");
    if (lua_isnil(L, -1)) {
        string errMsg = "Error in call to Sketch:arc{}. No p0 entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto p0 = toVector3(L, -1);
    lua_pop(L, 1); // dispose of p0 item
    //
    lua_getfield(L, 1, "p1");
    if (lua_isnil(L, -1)) {
        string errMsg = "Error in call to Sketch:arc{}. No p1 entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto p1 = toVector3(L, -1);
    lua_pop(L, 1); // dispose of p1 item
    //
    lua_getfield(L, 1, "centre");
    if (lua_isnil(L, -1)) {
        string errMsg = "Error in call to Sketch:arc{}. No centre entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto pc = toVector3(L, -1);
    lua_pop(L, 1); // dispose of centre item
    //
    bool dashed = false;
    lua_getfield(L, 1, "dashed");
    if (!lua_isnil(L, -1)) {
        if (lua_isboolean(L, -1)) {
            dashed = to!bool(lua_toboolean(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:arc{}. Expected a bool for dashed.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of dashed item
    //
    my_sketch.arc(p0, p1, pc, dashed);
    return 0;
} // end arcSketch()

extern(C) int polygonSketch(lua_State* L)
{
    auto my_sketch = checkObj!(Sketch, SketchMT)(L, 1); // first argument "this"
    lua_remove(L, 1); // remove first argument "this"
    // Now, I'm expecting a table of arguments as the only item on the stack.
    int narg = lua_gettop(L);
    if (narg == 0 || !lua_istable(L, 1)) {
        string errMsg = "Error in call to Sketch:polygon{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    //
    lua_getfield(L, 1, "points".toStringz());
    if ( lua_isnil(L, -1) ) {
        string errMsg = "Error in call to Sketch:polygon{}. No points entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    if (!lua_istable(L, -1)) {
        string errMsg = "Error in call to Sketch:polygon{}.; " ~
            "A table containing Vector3 points is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    // Expect Vector3 objects at array positions within that table.
    Vector3[] points;
    int position = 1;
    while (true) {
        lua_rawgeti(L, -1, position);
        if (lua_isnil(L, -1)) { lua_pop(L, 1); break; }
        auto p = toVector3(L, -1);
        lua_pop(L, 1);
        points ~= p;
        ++position;
    }
    lua_pop(L, 1); // dispose of points table
    if (points.length == 0) {
        string errMsg = "Error in call to Sketch:polygon{}. No valid Vector3 objects found.";
        luaL_error(L, errMsg.toStringz());
    }
    //
    bool fill = true;
    lua_getfield(L, 1, "fill");
    if (!lua_isnil(L, -1)) {
        if (lua_isboolean(L, -1)) {
            fill = to!bool(lua_toboolean(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:polygon{}. Expected a bool for fill.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of fill item
    //
    bool stroke = true;
    lua_getfield(L, 1, "stroke");
    if (!lua_isnil(L, -1)) {
        if (lua_isboolean(L, -1)) {
            stroke = to!bool(lua_toboolean(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:polygon{}. Expected a bool for stroke.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of stroke item
    //
    bool dashed = false;
    lua_getfield(L, 1, "dashed");
    if (!lua_isnil(L, -1)) {
        if (lua_isboolean(L, -1)) {
            dashed = to!bool(lua_toboolean(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:polygon{}. Expected a bool for dashed.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of dashed item
    //
    my_sketch.polygon(points, fill, stroke, dashed);
    return 0;
} // end polygonSketch()

extern(C) int textSketch(lua_State* L)
{
    auto my_sketch = checkObj!(Sketch, SketchMT)(L, 1); // first argument "this"
    lua_remove(L, 1); // remove first argument "this"
    // Now, I'm expecting a table of arguments as the only item on the stack.
    int narg = lua_gettop(L);
    if (narg == 0 || !lua_istable(L, 1)) {
        string errMsg = "Error in call to Sketch:text{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    //
    lua_getfield(L, 1, "point");
    if (lua_isnil(L, -1)) {
        string errMsg = "Error in call to Sketch:text{}. No point entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    Vector3 point = toVector3(L, -1);
    lua_pop(L, 1); // dispose of point item
    //
    string text = "";
    lua_getfield(L, 1, "text");
    if (lua_isstring(L, -1)) {
        text = to!string(lua_tostring(L, -1));
    } else {
        string errMsg = "Error in call to Sketch:text{}. Expected a string for text.";
        luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1); // dispose of text item
    //
    double angle = 0.0;
    lua_getfield(L, 1, "angle");
    if (!lua_isnil(L, -1)) {
        if (lua_isnumber(L, -1)) {
            angle = to!double(lua_tonumber(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:text{}. Expected a number for angle.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of angle item
    //
    string anchor = "middle";
    lua_getfield(L, 1, "anchor");
    if (!lua_isnil(L, -1)) {
        if (lua_isstring(L, -1)) {
            anchor = to!string(lua_tostring(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:text{}. Expected a string for anchor.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of anchor item
    //
    int fontSize = 10;
    lua_getfield(L, 1, "font_size");
    if (!lua_isnil(L, -1)) {
        if (lua_isnumber(L, -1)) {
            fontSize = to!int(lua_tointeger(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:text{}. Expected a number for font_size.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of font_size item
    //
    string colour = "black";
    lua_getfield(L, 1, "colour");
    if (!lua_isnil(L, -1)) {
        if (lua_isstring(L, -1)) {
            colour = to!string(lua_tostring(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:text{}. Expected a string for colour.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of colour item
    //
    string fontFamily = "sanserif";
    lua_getfield(L, 1, "font_family");
    if (!lua_isnil(L, -1)) {
        if (lua_isstring(L, -1)) {
            colour = to!string(lua_tostring(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:text{}. Expected a string for font_family.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of font_family item
    //
    my_sketch.text(point, text, angle, anchor, fontSize, colour, fontFamily);
    return 0;
} // end textSketch()

extern(C) int dotlabelSketch(lua_State* L)
{
    auto my_sketch = checkObj!(Sketch, SketchMT)(L, 1); // first argument "this"
    lua_remove(L, 1); // remove first argument "this"
    // Now, I'm expecting a table of arguments as the only item on the stack.
    int narg = lua_gettop(L);
    if (narg == 0 || !lua_istable(L, 1)) {
        string errMsg = "Error in call to Sketch:dotlabel{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    //
    lua_getfield(L, 1, "point");
    if (lua_isnil(L, -1)) {
        string errMsg = "Error in call to Sketch:dotlabel{}. No point entry found.";
        luaL_error(L, errMsg.toStringz());
    }
    auto point = toVector3(L, -1);
    lua_pop(L, 1); // dispose of point item
    //
    string label = "";
    lua_getfield(L, 1, "label");
    if (lua_isstring(L, -1)) {
        label = to!string(lua_tostring(L, -1));
    } else {
        string errMsg = "Error in call to Sketch:dotlabel{}. Expected a string for label.";
        luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1); // dispose of label item
    //
    string anchor = "middle";
    lua_getfield(L, 1, "anchor");
    if (!lua_isnil(L, -1)) {
        if (lua_isstring(L, -1)) {
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
    if (!lua_isnil(L, -1)) {
        if (lua_isnumber(L, -1)) {
            dotSize = to!double(lua_tonumber(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:dotlabel{}. Expected a number for dot_size.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of dot_size item
    //
    int fontSize = 10;
    lua_getfield(L, 1, "font_size");
    if (!lua_isnil(L, -1)) {
        if (lua_isnumber(L, -1)) {
            fontSize = to!int(lua_tointeger(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:text{}. Expected a number for font_size.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of font_size item
    //
    string colour = "black";
    lua_getfield(L, 1, "colour");
    if (!lua_isnil(L, -1)) {
        if (lua_isstring(L, -1)) {
            colour = to!string(lua_tostring(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:text{}. Expected a string for colour.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of colour item
    //
    string fontFamily = "sanserif";
    lua_getfield(L, 1, "font_family");
    if (!lua_isnil(L, -1)) {
        if (lua_isstring(L, -1)) {
            fontFamily = to!string(lua_tostring(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:text{}. Expected a string for font_family.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of font_family item
    //
    my_sketch.dotlabel(point, label, anchor, dotSize, fontSize, colour, fontFamily);
    return 0;
} // end dotlabelSketch()


extern(C) int ruleSketch(lua_State* L)
{
    auto my_sketch = checkObj!(Sketch, SketchMT)(L, 1); // first argument "this"
    lua_remove(L, 1); // remove first argument "this"
    // Now, I'm expecting a table of arguments as the only item on the stack.
    int narg = lua_gettop(L);
    if (narg == 0 || !lua_istable(L, 1)) {
        string errMsg = "Error in call to Sketch:rule{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    //
    string direction = "x";
    lua_getfield(L, 1, "direction");
    if (lua_isstring(L, -1)) {
        direction = to!string(lua_tostring(L, -1));
    } else {
        string errMsg = "Error in call to Sketch:rule{}. Expected a string for direction.";
        luaL_error(L, errMsg.toStringz());
    }
    lua_pop(L, 1); // dispose of direction item
    //
    double vmin = 0.0;
    lua_getfield(L, 1, "vmin");
    if (!lua_isnil(L, -1)) {
        if (lua_isnumber(L, -1)) {
            vmin = to!double(lua_tonumber(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:rule{}. Expected a number for vmin.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of vmin item
    //
    double vmax = 1.0;
    lua_getfield(L, 1, "vmax");
    if (!lua_isnil(L, -1)) {
        if (lua_isnumber(L, -1)) {
            vmax = to!double(lua_tonumber(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:rule{}. Expected a number for vmax.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of vmax item
    //
    double vtic = 0.2;
    lua_getfield(L, 1, "vtic");
    if (!lua_isnil(L, -1)) {
        if (lua_isnumber(L, -1)) {
            vtic = to!double(lua_tonumber(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:rule{}. Expected a number for vtic.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of vtic item
    //
    Vector3 anchorPoint = Vector3(0.0,0.0,0.0);
    lua_getfield(L, 1, "anchor_point");
    if (!lua_isnil(L, -1)) {
        anchorPoint = toVector3(L, -1);
    }
    lua_pop(L, 1); // dispose of anchorPoint item
    //
    double ticMarkSize = 0.02;
    lua_getfield(L, 1, "tic_mark_size");
    if (!lua_isnil(L, -1)) {
        if (lua_isnumber(L, -1)) {
            ticMarkSize = to!double(lua_tonumber(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:rule{}. Expected a number for tic_mark_size.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of tic_mark_size item
    //
    string numberFormat = "%.1f";
    lua_getfield(L, 1, "number_format");
    if (!lua_isnil(L, -1)) {
        if (lua_isstring(L, -1)) {
            numberFormat = to!string(lua_tostring(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:rule{}. Expected a string for number_format.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of number_format item
    //
    double textOffset = 0.06;
    lua_getfield(L, 1, "text_offset");
    if (!lua_isnil(L, -1)) {
        if (lua_isnumber(L, -1)) {
            textOffset = to!double(lua_tonumber(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:rule{}. Expected a number for text_offset.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of text_offset item
    //
    double textAngle = 0.0;
    lua_getfield(L, 1, "text_angle");
    if (!lua_isnil(L, -1)) {
        if (lua_isnumber(L, -1)) {
            textAngle = to!double(lua_tonumber(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:rule{}. Expected a number for text_angle.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of text_angle item
    //
    int fontSize = 10;
    lua_getfield(L, 1, "font_size");
    if (!lua_isnil(L, -1)) {
        if (lua_isnumber(L, -1)) {
            fontSize = to!int(lua_tointeger(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:rule{}. Expected a number for font_size.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of font_size item
    //
    my_sketch.rule(direction, vmin, vmax, vtic, anchorPoint,
                   ticMarkSize, numberFormat, textOffset, textAngle, fontSize);
    return 0;
} // end ruleSketch()

extern(C) int renderSketch(lua_State* L)
{
    auto my_sketch = checkObj!(Sketch, SketchMT)(L, 1); // first argument "this"
    lua_remove(L, 1); // remove first argument "this"
    // Now, I'm expecting a table of arguments as the only item on the stack.
    int narg = lua_gettop(L);
    if (narg == 0 || !lua_istable(L, 1)) {
        string errMsg = "Error in call to Sketch:render{}.; " ~
            "A table containing arguments is expected, but no table was found.";
        luaL_error(L, errMsg.toStringz);
    }
    //
    Path my_path;
    lua_getfield(L, 1, "path");
    if (!lua_isnil(L, -1)) { my_path = checkPath(L, -1); }
    lua_pop(L, 1); // dispose of path item
    ParametricSurface my_surf;
    lua_getfield(L, 1, "surf");
    if (!lua_isnil(L, -1)) { my_surf = checkSurface(L, -1); }
    lua_pop(L, 1); // dispose of surf item
    ParametricVolume my_vol;
    lua_getfield(L, 1, "volume");
    if (!lua_isnil(L, -1)) { my_vol = checkVolume(L, -1); }
    lua_pop(L, 1); // dispose of volume item
    if ((my_path is null) && (my_surf is null) && (my_vol is null)) {
        string errMsg = "Error in call to Sketch:render{}. " ~
            "Could not find a path nor surf nor volume element.";
        luaL_error(L, errMsg.toStringz());
    }
    //
    size_t n = 30;
    lua_getfield(L, 1, "n");
    if (!lua_isnil(L, -1)) {
        if (lua_isnumber(L, -1)) {
            n = to!size_t(lua_tointeger(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:render{}. Expected a number for n.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of n item
    //
    bool fill = true;
    lua_getfield(L, 1, "fill");
    if (!lua_isnil(L, -1)) {
        if (lua_isboolean(L, -1)) {
            fill = to!bool(lua_toboolean(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:render{}. Expected a bool for fill.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of fill item
    //
    bool facets = false;
    lua_getfield(L, 1, "facets");
    if (!lua_isnil(L, -1)) {
        if (lua_isboolean(L, -1)) {
            facets = to!bool(lua_toboolean(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:render{}. Expected a bool for facets.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of facets item
    //
    bool stroke = true;
    lua_getfield(L, 1, "stroke");
    if (!lua_isnil(L, -1)) {
        if (lua_isboolean(L, -1)) {
            stroke = to!bool(lua_toboolean(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:render{}. Expected a bool for stroke.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of stroke item
    //
    bool dashed = false;
    lua_getfield(L, 1, "dashed");
    if (!lua_isnil(L, -1)) {
        if (lua_isboolean(L, -1)) {
            dashed = to!bool(lua_toboolean(L, -1));
        } else {
            string errMsg = "Error in call to Sketch:render{}. Expected a bool for dashed.";
            luaL_error(L, errMsg.toStringz());
        }
    }
    lua_pop(L, 1); // dispose of dashed item
    //
    if (my_path) { my_sketch.render(my_path, dashed, n); }
    if (my_surf) { my_sketch.render(my_surf, fill, stroke, dashed, n, facets); }
    if (my_vol) { my_sketch.render(my_vol, fill, stroke, dashed, n, facets); }
    return 0;
} // end renderSketch()


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
    lua_pushcfunction(L, &arcSketch);
    lua_setfield(L, -2, "arc");
    lua_pushcfunction(L, &textSketch);
    lua_setfield(L, -2, "text");
    lua_pushcfunction(L, &dotlabelSketch);
    lua_setfield(L, -2, "dotlabel");
    lua_pushcfunction(L, &polylineSketch);
    lua_setfield(L, -2, "polyline");
    lua_pushcfunction(L, &polygonSketch);
    lua_setfield(L, -2, "polygon");
    lua_pushcfunction(L, &ruleSketch);
    lua_setfield(L, -2, "rule");
    lua_pushcfunction(L, &renderSketch);
    lua_setfield(L, -2, "render");

    lua_setglobal(L, SketchMT.toStringz);
} // end registerSketch()
