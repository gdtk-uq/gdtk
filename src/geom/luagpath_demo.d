/**
 * luagpath_demo.d Demonstrate some of the behaviour of the Path primitives.
 *
 * Author: Rowan G. and Peter J.
 * Version: 2015-02-22
 */

import std.stdio;
import std.conv;
import std.string;
import util.lua;
import geom;
import gpath;
import luageom;
import luagpath;

void main()
{
    writeln("Begin demonstration of LuaD connection to Paths.");
    auto L = luaL_newstate();
    luaL_openlibs(L);
    registerVector3(L);
    registerPaths(L);
    string test_code = `
-- Add a couple of points and alter their data.
a = Vector3:new{x=1.0, y=2.0}
b = Vector3:new(0.0, 5.0, 4.0)
ab = Line:new{a, b}
print("ab= ", ab)
print("Try evaluating a point midway on the line.")
pt = ab(0.5)
print("pt= ", pt)
print("Or with an eval.")
pt2 = ab:eval(0.5)
print("pt2= ", pt2)
--
print("Arc")
a = Vector3:new(2.0, 2.0, 0.0)
b = Vector3:new(1.0, 2.0, 1.0)
c = Vector3:new(1.0, 2.0, 0.0)
abc = Arc:new{a, b, c}
d = abc(0.5)
print("d=", d, "expected approximately Vector3(1.7071068, 2.0, 0.7071068)")
--
print("Arc3")
a = Vector3:new(2.0, 2.0, 0.0)
b = Vector3:new(1.0, 2.0, 1.0)
m = Vector3:new(1.7071068, 2.0, 0.7071068)
amb = Arc3:new{a, m, b}
dd = amb(0.5)
print("dd=", dd, "expected approximately Vector3(1.7071068, 2.0, 0.7071068)")
--
print("Bezier")
adb = Bezier:new{points={a, d, b}}
e = adb(0.5)
print("e=", e, "expected approximately Vector3(1.60355, 2, 0.603553)")
--
print("Polyline")
polyline = Polyline:new{abc, Line:new{b,c}}
print("polyline= ", polyline)
f = polyline(0.5)
print("polyline(0.25)= ", polyline(0.25))
print("polyline(0.5)= ", f, "expected approximately Vector3(1.28154, 2, 0.95955)")
print("polyline(0.75)= ", polyline(0.75))
--
print("LuaFnPath")
function myLuaFunction(t)
   -- Straight line from 0,0,0 to 1.0,2.0,3.0
   return {t, 2*t, 3*t}
end
myPath = LuaFnPath:new{"myLuaFunction"}
print("myLine= ", myPath)
g = myPath(0.5)
print("myPath(0.5)= ", g)
--
print("ArcLengthParameterizedPath")
p0 = Vector3:new(0.0, 0.0, 0.0)
p1 = Vector3:new(1.0, 1.0, 1.0)
p2 = Vector3:new(4.0, 4.0, 4.0)
alpp = ArcLengthParameterizedPath:new{Bezier:new{points={p0, p1, p2}}}
print("alpp=", alpp)
print("alpp(0.5)=", alpp(0.5))
--
print("SubrangedPath")
srp = SubRangedPath:new{polyline, t0=1.0, t1=0.0} -- effectively reversed
print("srp=", srp)
print("srp(0.25)=", srp(0.25))
print("srp(0.50)=", srp(0.50))
print("srp(0.75)=", srp(0.75))
--
print("Fiddle with parameter limits.")
srp:t1(0.8)
srp2 = srp:copy()
srp2:t0(0.2)
print("srp:t0()= ", srp:t0(), "srp:t1()= ", srp:t1())
print("srp2:t0()= ", srp2:t0(), "srp2:t1()= ", srp2:t1())
--
print("ReversedPath")
rp = ReversedPath:new{polyline}
print("rp=", rp)
print("rp(0.25)=", rp(0.25))
print("rp(0.50)=", rp(0.50))
print("rp(0.75)=", rp(0.75))
--
print("Spline (Polyline)")
pnts = {Vector3:new{0.0, -1.0, 0.0},
	Vector3:new{-1.0, 0.0, 0.0},
	Vector3:new{0.0, 1.0, 0.0},
	Vector3:new{1.0, 0.0, 0.0},
	Vector3:new{0.0, -1.0, 0.0}}
circle = Spline:new{points=pnts}
print("circle=", circle)
print("circle(5.0/8)=", circle(5.0/8))
    `;
    if ( luaL_dostring(L, toStringz(test_code)) != 0 ) {
	writeln("There was a problem interpreting the test code.");
	writeln(to!string(lua_tostring(L, -1)));
    }
    writeln("Done with luagpath_demo.");
}
