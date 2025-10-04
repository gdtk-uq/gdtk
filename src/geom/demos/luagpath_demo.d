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
import geom.luawrap;

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
b = Vector3:new{x=0.0, y=5.0, z=4.0}
ab = Line:new{p0=a, p1=b}
print("ab= ", ab)
print("Try evaluating a point midway on the line.")
pt = ab(0.5)
print("pt= ", pt)
print("Or with an eval.")
pt2 = ab:eval(0.5)
print("pt2= ", pt2)
--
print("Arc")
a = Vector3:new{x=2.0, y=2.0, z=0.0}
b = Vector3:new{x=1.0, y=2.0, z=1.0}
c = {x=1.0, y=2.0, z=0.0} -- should also accept a table with names coordinates
abc = Arc:new{p0=a, p1=b, centre=c}
d = abc(0.5)
print("d=", d, "expected approximately Vector3([1.7071068, 2.0, 0.7071068])")
--
print("Arc3")
a = Vector3:new{x=2.0, y=2.0, z=0.0}
b = Vector3:new{x=1.0, y=2.0, z=1.0}
m = Vector3:new{x=1.7071068, y=2.0, z=0.7071068}
amb = Arc3:new{p0=a, pmid=m, p1=b}
dd = amb(0.5)
print("dd=", dd, "expected approximately Vector3([1.7071068, 2.0, 0.7071068])")
--
print("Helix")
axis0 = Vector3:new{x=0}
axis1 = Vector3:new{x=1}
pstart = Vector3:new{y=1}
pend = Vector3:new{x=1, z=1}
h1 = Helix:new{point_start=pstart, point_end=pend, axis0=axis0, axis1=axis1}
ddd = h1(0.5)
print("ddd=", ddd, "expected approximately Vector3([0.5, 0.7071068, 0.7071068])")
-- Define equivalent Helix using fundamental parameters.
h2 = Helix:new{a0=Vector3:new{x=0.0}, a1=Vector3:new{x=1.0},
               xlocal=Vector3:new{y=1.0},
               r0=1.0, r1=1.0, dtheta=math.pi/2};
ddd2 = h2(0.5);
print("ddd2=", ddd2, "expected approximately Vector3([0.5, 0.7071068, 0.7071068])")
--
print("Bezier")
adb = Bezier:new{points={a, d, b}}
e = adb(0.5)
print("e=", e, "expected approximately Vector3([1.60355, 2, 0.603553])")
--
mypth = Arc3:new{p0=Vector3:new{x=0.0,y=1.0}, pmid=Vector3:new{x=0.5,y=1.2}, p1=Vector3:new{x=1.0,y=1.0}}
myps = Vector3:new{x=0.5,y=0.5}
mydir = Vector3:new{x=0.0,y=1.0}
found, t = mypth:intersect2D{ps=myps, d=mydir, nseg=10}
print("intersection on Arc3: found=", found, "t=", t)
--
print("NURBS")
pts = {Vector3:new{x=-4, y=-4},
       Vector3:new{x=-2, y=4},
       Vector3:new{x=2, y=-4},
       Vector3:new{x=4, y=4},
       Vector3:new{x=3.778, y=1.836, z=2.933},
       Vector3:new{x=2.772, y=-3.875, z=1.736}}
w = {1.0, 1.0, 5.0, 1.0, 1.0, 4.0}
U = {0.0, 0.0, 0.0, 0.375, 0.5, 0.625, 1.0, 1.0, 1.0}
p = 2
nrb = NURBS:new{points=pts, weights=w, knots=U, degree=p}
np = nrb(0.6)
print("np= ", np, "expected appoximately Vector3([3.782, 2.939, 0.435])")
--
print("Polyline")
polyline = Polyline:new{segments={abc, Line:new{p0=b,p1=c}}}
print("polyline= ", polyline)
f = polyline(0.5)
print("polyline(0.25)= ", polyline(0.25))
print("polyline(0.5)= ", f, "expected approximately Vector3([1.28154, 2, 0.95955])")
print("polyline(0.75)= ", polyline(0.75))
--
print("SVGPath")
svgpth = SVGPath:new{path="M3.0,3.0;L4.0,3.0;v1.0;h-1.0;Z"}
f1 = svgpth(0.5)
print("svgpth(0.25)= ", svgpth(0.25), "expected Vector3([4.0, 3.0, 0.0])")
print("svgpth(0.5)= ", f1, "expected Vector3([4.0, 4.0, 0.0])")
print("svgpth(0.75)= ", svgpth(0.75), "expected Vector3([3.0, 4.0, 0.0])")
--
print("LuaFnPath")
function myLuaFunction(t)
   -- Straight line from 0,0,0 to 1.0,2.0,3.0
   return {x=t, y=2*t, z=3*t}
end
myPath = LuaFnPath:new{luaFnName="myLuaFunction"}
print("myLine= ", myPath)
g = myPath(0.5)
print("myPath(0.5)= ", g)
--
print("ArcLengthParameterizedPath")
p0 = Vector3:new{x=0.0, y=0.0, z=0.0}
p1 = Vector3:new{x=1.0, y=1.0, z=1.0}
p2 = Vector3:new{x=4.0, y=4.0, z=4.0}
alpp = ArcLengthParameterizedPath:new{underlying_path=Bezier:new{points={p0, p1, p2}}}
print("alpp=", alpp)
print("alpp(0.5)=", alpp(0.5))
--
print("SubrangedPath")
srp = SubRangedPath:new{underlying_path=polyline, t0=1.0, t1=0.0} -- effectively reversed
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
rp = ReversedPath:new{underlying_path=polyline}
print("rp=", rp)
print("rp(0.25)=", rp(0.25))
print("rp(0.50)=", rp(0.50))
print("rp(0.75)=", rp(0.75))
--
print("Spline (Polyline)")
pnts = {Vector3:new{ x=0.0, y=-1.0},
        Vector3:new{x=-1.0,  y=0.0},
        Vector3:new{ x=0.0,  y=1.0},
        Vector3:new{ x=1.0,  y=0.0},
        Vector3:new{ x=0.0, y=-1.0}}
circle = Spline:new{points=pnts}
print("circle=", circle)
print("circle(5.0/8)=", circle(5.0/8))
--
print("Spline2 (Polyline)")
circle2 = Spline2:new{filename="sample-data/test-spline2.dat"}
print("circle2=", circle2)
print("circle2(5.0/8)=", circle2(5.0/8))
--
print("XSpline")
cspline = XSpline:new{xs={1,2,3,4}, ys={1.5,2.5,1.5,2.5}}
print("cspline=", cspline)
print("cspline(0.5)=", cspline(0.5))
--
print("XSpline2")
cspline2 = XSpline2:new{filename="sample-data/test-xspline2.dat"}
print("cspline2=", cspline2)
print("cspline2(0.5)=", cspline2(0.5))
--
print("XSplineLsq")
xd = {}; yd = {}
for i=0,100 do
   x = 2.0*math.pi*i/100.0
   y = math.sin(x)
   xd[#xd+1] = x
   yd[#yd+1] = y
end
cspline3 = XSplineLsq:new{xd=xd, yd=yd, nseg=10}
print("cspline3=", cspline3)
print("cspline3(0.5)=", cspline3(0.5))
--
print("XSplineLsq2")
cspline4 = XSplineLsq2:new{filename="sample-data/test-xsplinelsq.dat", nseg=10}
print("cspline4=", cspline4)
print("cspline4(0.5)=", cspline4(0.5))
--
print("TranslatedPath")
a = Vector3:new{x=2.0, y=0.0}
b = Vector3:new{x=0.0, y=2.0}
c = Vector3:new{x=0.0, y=0.0}
abc = Arc:new{p0=a, p1=b, centre=c}
abc_translated = TranslatedPath:new{original_path=abc, shift=Vector3:new{z=0.33333}}
print("abc_translated=", abc_translated)
print("abc_translated(0.5)=", abc_translated(0.5))
--
print("MirrorImagePath")
p0 = Vector3:new{x=1, y=0}
p1 = Vector3:new{x=0.7071, y=0.7071}
p2 = Vector3:new{x=0, y=1}
original_path = Bezier:new{points={p0, p1, p2}}
mi_path = MirrorImagePath:new{original_path=original_path,
                              point=Vector3:new{x=1.0, y=0.0},
                              normal=Vector3:new{x=1.0, y=0.0}}
print("mi_path=", mi_path)
print("original_path(0.5)=", original_path(0.5))
print("mi_path(0.5)=", mi_path(0.5))
--
print("RotatedAboutZAxisPath")
a = Vector3:new{x=2.0, y=0.0}
b = Vector3:new{x=0.0, y=2.0}
c = Vector3:new{x=0.0, y=0.0}
abc = Arc:new{p0=a, p1=b, centre=c}
abc_rotated = RotatedAboutZAxisPath:new{original_path=abc, angle=math.pi/4}
print("abc_rotated=", abc_rotated)
print("abc_rotated(1.0)=", abc_rotated(1.0))
    `;
    if ( luaL_dostring(L, toStringz(test_code)) != 0 ) {
        writeln("There was a problem interpreting the test code.");
        writeln(to!string(lua_tostring(L, -1)));
    }
    writeln("Done with luagpath_demo.");
}
