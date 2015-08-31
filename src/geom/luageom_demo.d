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
import luageom;

void main()
{
    writeln("Begin demo of wrapped D Vector3 struct for use in Lua.");
    auto L = luaL_newstate();
    luaL_openlibs(L);
    registerVector3(L);
    string test_code = `
-- Add some points and manipulate them.
print("Test some constructors.")
a = Vector3:new()
print("a= ", a)
b = Vector3:new(7.0, 3.0, -2.5)
print("b= ", b)
c = Vector3:new{4.0, 1.2, 17}
print("c= ", c)
d = Vector3:new{x=5, z="9", y=78.6, label="some crap"}
print("d= ", d)
e = Vector3:new(1.0)
print("e= ", e)
f = Vector3:new(true, nil, 6.0)
print("f= ", f)
print("Sorry, you gave me values I couldn't use in slots 0 and 1.")
print("Change the x value of f.")
f:x(5.4)
print("f= ", f)
g = -f
print("g= ", g)
assert(g:x() == -f:x()); assert(g:y() == -f:y()); assert(g:z() == -f:z())
h = a + b
assert(h:x() == a:x() + b:x()); assert(h:y() == a:y() + b:y()); assert(h:z() == a:z() + b:z())
i = unit(h)
h:normalize()
print("After normalizing, h=", h)
assert(h:x() == i:x()); assert(h:y() == i:y()); assert(h:z() == i:z())
assert(vabs(i) == 1.0)
j = dot(g, f)
print("j= ", j)
k = cross(g, f)
print("k= ", k)
    `;
    if ( luaL_dostring(L, toStringz(test_code)) != 0 ) {
	writeln("There was a problem interpreting the test code.");
	writeln(to!string(lua_tostring(L, -1)));
    }
    writeln("End luageom_demo.");
}

