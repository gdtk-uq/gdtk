# path-test.py
#
# PJ, 2020-02-05

import math
from gdtk.geom.vector3 import Vector3
from gdtk.geom.path import *
import numpy

def approxEqual(a, b):
    result = math.isclose(a, b, rel_tol=1.0e-2, abs_tol=1.0e-5)
    # print("a=",a, "b=",b, "rel=",(a-b)/b, "abs=",a-b, "result=",result)
    return result

print("Begin test of Path classes...")
a = Vector3(0.0, 2.0)
print("a=", a)
assert approxEqual(a.x, 0.0), "Vector3 list constructor x-component"
assert approxEqual(a.y, 2.0), "Vector3 list constructor y-component"
assert approxEqual(a.z, 0.0), "Vector3 list constructor z-component"
b = Vector3(2.0, 0.0)
print("b=", b)
assert approxEqual(b.x, 2.0), "Vector3 dict constructor x-component"
assert approxEqual(b.y, 0.0), "Vector3 dict constructor y-component"
assert approxEqual(b.z, 0.0), "Vector3 dict constructor z-component"
line_ab = Line(p0=a, p1=b)
print("line_ab=", line_ab)
c = line_ab(0.5)
print("c=line_ab(0.5)=", c)
assert approxEqual(c.x, 1.0), "Line evaluation x-component"
assert approxEqual(c.y, 1.0), "Line evaluation y-component"
assert approxEqual(c.z, 0.0), "Line evaluation z-component"
print("line_ab.length()=", line_ab.length())
assert approxEqual(line_ab.length(), math.sqrt(8.0)), "Line length()"

c = Vector3(0.0, 0.0, 0.0)
arc_abc = Arc(a, b, c)
print("arc_abc=", arc_abc)
d = arc_abc(0.5)
print("d=line_abc(0.5)=", d)
assert approxEqual(d.x, math.sqrt(2)), "Arc evaluation x-component"
assert approxEqual(d.y, math.sqrt(2)), "Arc evaluation y-component"
assert approxEqual(d.z, 0.0), "Arc evaluation z-component"
print("arc_abc.length()=", arc_abc.length())
assert approxEqual(arc_abc.length(), math.pi), "Arc length()"

adb = Bezier([a, d, b])
adb5 = adb(0.5)
assert approxEqual(adb5.x, 1.2071), "Bezier evaluation x-component"
assert approxEqual(adb5.y, 1.2071), "Bezier evaluation y-component"
assert approxEqual(adb5.z, 0.0), "Bezier evaluation z-component"
print("Bezier adb(0.5)=", adb5)
print("Length of adb=", adb.length())

e = Polyline([Line(p0=[-math.pi,2.0],p1=[0.0,2.0]), arc_abc])
print("e=", e)
f = e(0.0)
assert approxEqual(f.x, -math.pi), "Polyline evaluation x-component"
assert approxEqual(f.y, 2.0), "Polyline evaluation y-component"
assert approxEqual(f.z, 0.0), "Polyline evaluation z-component"
print("f=e(0.0)=", f)
f = e(0.5)
assert approxEqual(f.x, 0.0), "Polyline evaluation x-component"
assert approxEqual(f.y, 2.0), "Polyline evaluation y-component"
assert approxEqual(f.z, 0.0), "Polyline evaluation z-component"
print("f=e(0.5)=", f)
f = e(1.0)
assert approxEqual(f.x, 2.0), "Polyline evaluation x-component"
assert approxEqual(f.y, 0.0), "Polyline evaluation y-component"
assert approxEqual(f.z, 0.0), "Polyline evaluation z-component"
print("f=e(1.0)=", f)

g = ArcLengthParameterizedPath(e)
print("g=", g)
h = g(0.5)
print("h=g(0.5)=", h)
assert approxEqual(h.x, 0.0), "Polyline evaluation x-component"
assert approxEqual(h.y, 2.0), "Polyline evaluation y-component"
assert approxEqual(h.z, 0.0), "Polyline evaluation z-component"

pnts = [Vector3(x, math.sin(x))
        for x in numpy.linspace(0.0, 2*math.pi, 10)]
s = Spline(points=pnts)
print("s=", s)
z = s(0.5)
print("s(0.5)=", z)
assert approxEqual(z.x, math.pi), "Spline zero-crossing x-component"
assert approxEqual(z.y, 0.0), "Spline zero-crossing y-component"
assert approxEqual(z.z, 0.0), "Spline zero-crossing z-component"
p = s(0.25)
print("s(0.25)=", p)
assert approxEqual(p.x, math.pi/2), "Spline peak x-component"
assert approxEqual(p.y, 1.0), "Spline peak y-component"
assert approxEqual(p.z, 0.0), "Spline peak z-component"

print("Done.")
