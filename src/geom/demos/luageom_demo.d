/**
 * luageom_demo.d
 * Shows the wrapped Vector3 in action.
 *
 * Author: Rowan G. and Peter J.
 * Date: 2015-02-21
 */

import std.stdio;
import std.conv;
import std.string;
import util.lua;
import geom.luawrap;

void main()
{
    writeln("Begin demo of wrapped D Vector3 struct for use in Lua.");
    auto L = luaL_newstate();
    luaL_openlibs(L);
    registerVector3(L);
    string test_code = `
-- Add some points and manipulate them.
print("Test some constructors.")
a = Vector3:new{}
print("a= ", a)
b = Vector3:new{x=7.0, y=3.0, z=-2.5}
print("b= ", b)
print("Get the 'z' value of b")
print("b.z= ", b.z)
bb = Vector3:new(b)
print("bb=", bb)
bbb = Vector3:new{x=7.0}
print("bbb=", bbb)
c = Vector3:new{x=4.0, y=1.2, z=17}
print("c= ", c)
cc = Vector3:new{c}
print("cc=", cc)
d = Vector3:new{x=5, z="9", y=78.6, label="some crap"}
print("d= ", d)
e = Vector3:new{x=1.0, y=0.0}
print("e= ", e)
function bad_fn()
   f = Vector3:new{true}
end
if pcall(bad_fn) then
   print("Sorry, you gave me a value I couldn't use.")
end
f = Vector3:new{}
print("f= ", f)
print("Change the x value of f.")
f.x = 5.4
print("f= ", f)
g = -f
print("g= ", g)
assert(g.x == -f.x); assert(g.y == -f.y); assert(g.z == -f.z)
h = a + b
assert(h.x == a.x + b.x); assert(h.y == a.y + b.y); assert(h.z == a.z + b.z)
i = unit(h)
print("i=", i)
h:normalize()
print("After normalizing, h=", h)
print("Check components")
assert(h.x == i.x); assert(h.y == i.y); assert(h.z == i.z)
print("Check magnitude")
print("vabs(i)=", vabs(i))
assert(math.abs(vabs(i)-1.0) < 1.0e-9)
print("Exercise dot product")
j = dot(g, f)
print("j= ", j)
k = cross(g, f)
print("k= ", k)
m = Vector3:new{x=1.0}
m:rotateAboutZAxis(math.pi/2)
print("rotated m=", m, "expected Vector3([0.0, 1.0, 0.0])")
n = Vector3:new{x=1.0,y=0.0,z=0.0}
n:mirrorImage(Vector3:new{x=2.0}, Vector3:new{x=1.0})
print("mirror-image n=", n, "expected Vector3([3.0, 0.0, 0.0])")
--
print("Transform from one frame to another.")
n = Vector3:new{x=1.0,y=1.0,z=0.0}; n = unit(n)
t1 = Vector3:new{x=-1.0,y=1.0,z=0.0}; t1 = unit(t1)
t2 = cross(n, t1)
h = Vector3:new{x=1.0,y=0.0,z=1.0}
h_ref = Vector3:new{x=h.x,y=h.y,z=h.z}
print("original h = ", h)
h:transformToLocalFrame(n, t1, t2)
print("in local frame h = ", h)
h:transformToGlobalFrame(n, t1, t2)
print("back to global frame h = ", h)
--
print("Try calling cell geometry calculation functions.")
-- 2019-07-07 Allow tables as well as Vector3 objects for the following tests.
t = quadProperties{p0={x=0.0, y=0.0}, p1={x=1.0, y=0.0},
                   p2=Vector3:new{x=1.0, y=1.0}, p3=Vector3:new{x=0.0, y=1.0}}
print("area=", t.area, "centroid=", t.centroid, "n=", t.n, "t1=", t.t1, "t2=", t.t2)
t = hexCellProperties{p0={x=0.0, y=0.0, z=0.0},
                      p1=Vector3:new{x=1.0, y=0.0, z=0.0},
                      p2=Vector3:new{x=1.0, y=1.0, z=0.0},
                      p3=Vector3:new{x=0.0, y=1.0, z=0.0},
                      p4=Vector3:new{x=0.0, y=0.0, z=1.0},
                      p5=Vector3:new{x=1.0, y=0.0, z=1.0},
                      p6=Vector3:new{x=1.0, y=1.0, z=1.0},
                      p7={x=0.0, y=1.0, z=1.0}}
print("volume=", t.volume, "centroid=", t.centroid, "iLen=", t.iLen, "jLen=", t.jLen, "kLen=", t.kLen)
    `;
    if ( luaL_dostring(L, toStringz(test_code)) != 0 ) {
        writeln("There was a problem interpreting the test code.");
        writeln(to!string(lua_tostring(L, -1)));
    }
    writeln("End luageom_demo.");
}

