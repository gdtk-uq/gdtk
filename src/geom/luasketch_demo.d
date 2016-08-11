/**
 * luasketch_demo.d  Demonstrate the use of the Sketch class from Lua.
 *
 * Authors: Peter J. and Rowan G.
 * Version: 2016-08-11
 */

import std.stdio;
import std.conv;
import std.string;
import util.lua;
import geom;
import gpath;
import surface;
import sketch;
import luageom;
import luagpath;
import luasurface;
import luasketch;

void main()
{
    writeln("Begin demonstration of LuaD connection to Sketch class.");
    auto L = luaL_newstate();
    luaL_openlibs(L);
    registerVector3(L);
    registerPaths(L);
    registerSurfaces(L);
    registerSketch(L);
    string test_code = `
-- Add a couple of points and alter their data.
s = Sketch:new{renderer="svg", projection="xyortho"}
print("s=", s)
s:set{canvas={0.0,0.0,120.0,120.0}, viewport={-2.0,-2.0,2.0,2.0}}
a = Vector3:new{x=2.0, y=0.0}
b = Vector3:new{x=1.0, y=1.0}
c = Vector3:new{x=1.0, y=0.0}
abc = Arc:new{p0=a, p1=b, centre=c}
s:start{file_name="test.svg"}
s:set{line_width=0.1}
s:line{p0=Vector3:new{x=-1.0,y=1.0,z=0.0}, p1=Vector3:new{x=1.0,y=-1.0,z=0.0}}
s:render{path=abc}
s:set{fill_colour="green"}
s:text{point=Vector3:new{x=-1.0,y=1.5},text="A sample, just to see",font_size=20}
p00 = Vector3:new{x=0.0, y=0.1}
p10 = Vector3:new{x=1.0, y=0.1}
p11 = Vector3:new{x=1.0, y=1.1}
p01 = Vector3:new{x=0.0, y=1.1}
s:set{line_width=0.3}
-- s:polyline{points={p00,p10,p11,p01}}
-- s:polygon{points={p00,p10,p11,p01}, dashed=true}
my_patch = CoonsPatch:new{p00=p00, p10=p10, p11=p11, p01=p01}
s:dotlabel{point=p00, label="p00"}
s:dotlabel{point=p10, label="p10"}
s:dotlabel{point=p11, label="p11"}
s:dotlabel{point=p01, label="p01", anchor="left", dot_size=1.0}
s:render{surf=my_patch};
s:finish{}
    `;
    if ( luaL_dostring(L, toStringz(test_code)) != 0 ) {
	writeln("There was a problem interpreting the test code.");
	writeln(to!string(lua_tostring(L, -1)));
    }
    writeln("Done with luasketch_demo.");
}
