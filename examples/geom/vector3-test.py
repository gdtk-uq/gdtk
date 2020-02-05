# vector3-test.py
#
# PJ, 2020-02-05

import math
from eilmer.geom.vector3 import Vector3, cross

def approxEqual(a, b):
    result = math.isclose(a, b, rel_tol=1.0e-2, abs_tol=1.0e-5)
    # print("a=",a, "b=",b, "rel=",(a-b)/b, "abs=",a-b, "result=",result) 
    return result

print("Begin test of Vector3 class...")
a = Vector3([1.0, 2.0])
print("a=", a)
assert approxEqual(a.x, 1.0), "Vector3 list constructor x-component" 
assert approxEqual(a.y, 2.0), "Vector3 list constructor y-component" 
assert approxEqual(a.z, 0.0), "Vector3 list constructor z-component"
print("abs(a)=", abs(a))
assert approxEqual(abs(a), math.sqrt(5.0)), "Vector3 magnitude"
b = Vector3({'x':1.0, 'z':2.0})
print("b=", b)
assert approxEqual(b.x, 1.0), "Vector3 dict constructor x-component" 
assert approxEqual(b.y, 0.0), "Vector3 dict constructor y-component" 
assert approxEqual(b.z, 2.0), "Vector3 dict constructor z-component" 
c = Vector3(z=3.0, x=1, y=2.0)
print("c=", c)
assert approxEqual(c.x, 1.0), "Vector3 numbers constructor x-component" 
assert approxEqual(c.y, 2.0), "Vector3 numbers constructor y-component" 
assert approxEqual(c.z, 3.0), "Vector3 numbers constructor z-component" 
d = a+b
print("d=a+b=", d)
assert approxEqual(d.x, 2.0), "Vector3 __add__ x-component" 
assert approxEqual(d.y, 2.0), "Vector3 __add__ y-component" 
assert approxEqual(d.z, 2.0), "Vector3 __add__ z-component" 
e = 2*a
print("e=2*a=", e)
assert approxEqual(e.x, 2.0), "Vector3 __rmul__ x-component" 
assert approxEqual(e.y, 4.0), "Vector3 __rmul__ y-component" 
assert approxEqual(e.z, 0.0), "Vector3 __rmul__ z-component" 
f = a*3
print("f=a*3=", f)
assert approxEqual(f.x, 3.0), "Vector3 __mul__ x-component" 
assert approxEqual(f.y, 6.0), "Vector3 __mul__ y-component" 
assert approxEqual(f.z, 0.0), "Vector3 __mul__ z-component" 
d += a
print("d+=a; d=", d)
assert approxEqual(d.x, 3.0), "Vector3 __iadd__ x-component" 
assert approxEqual(d.y, 4.0), "Vector3 __iadd__ y-component" 
assert approxEqual(d.z, 2.0), "Vector3 __iadd__ z-component" 
d *= 2
print("d*=2; d=", d)
assert approxEqual(d.x, 6.0), "Vector3 __imul__ x-component" 
assert approxEqual(d.y, 8.0), "Vector3 __imul__ y-component" 
assert approxEqual(d.z, 4.0), "Vector3 __imul__ z-component" 
d /= 2
print("d/=2; d=", d)
assert approxEqual(d.x, 3.0), "Vector3 __itruediv__ x-component" 
assert approxEqual(d.y, 4.0), "Vector3 __itruediv__ y-component" 
assert approxEqual(d.z, 2.0), "Vector3 __itruediv__ z-component"
a.normalize()
print("after normalize() call, a=", a)
assert approxEqual(a.x, 1.0/math.sqrt(5)), "Vector3 normalize() x-component" 
assert approxEqual(a.y, 2.0/math.sqrt(5)), "Vector3 normalize() y-component" 
assert approxEqual(a.z, 0.0), "Vector3 normalize() z-component"
g = cross(Vector3(1.0, 0.0, 0.0), Vector3(0.0, 1.0, 0.0))
print("g=", g)
assert approxEqual(g.x, 0.0), "Vector3 cross product x-component" 
assert approxEqual(g.y, 0.0), "Vector3 cross product y-component" 
assert approxEqual(g.z, 1.0), "Vector3 cross product z-component"
# Make local right-handed system
x = Vector3(1.0, 0.0, 0.0)
y = Vector3(0.0, 1.0, 0.0)
z = Vector3(0.0, 0.0, 1.0)
n = -y
t1 = -x
t2 = -z
print("local RH system")
print("n=", n, "t1=", t1, "t2=", t2)
h = Vector3(3.0, 4.0, 1.0)
print("in global frame, h=", h)
h.transform_to_local_frame(n, t1, t2)
print("in local frame, h=", h)
assert approxEqual(h.x, -4.0), "Vector3 transform to local x-component" 
assert approxEqual(h.y, -3.0), "Vector3 transform to local y-component" 
assert approxEqual(h.z, -1.0), "Vector3 transform to local z-component"
h.transform_to_global_frame(n, t1, t2)
print("in global frame, h=", h)
assert approxEqual(h.x, 3.0), "Vector3 transform to global x-component" 
assert approxEqual(h.y, 4.0), "Vector3 transform to global y-component" 
assert approxEqual(h.z, 1.0), "Vector3 transform to global z-component"
print("Done.")
