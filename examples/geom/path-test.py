# path-test.py
#
# PJ, 2020-02-05

import math
from eilmer.geom.vector3 import Vector3
from eilmer.geom.path import Line, Arc, Polyline

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

e = Polyline([Line(p0=[-math.pi,2.0],p1=[0.0,2.0]), arc_abc])
print("e=", e)
f = e(0.5)
print("e(0.5)=", e(0.5))

print("Done.")
