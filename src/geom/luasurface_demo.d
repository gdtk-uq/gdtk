/**
 * luasurface_demo.d
 * Demonstrate the wrapped Surface objects.
 *
 * Authors: Rowan G. and Peter J.
 * Version: 2015-02-24
 */

import std.stdio;
import std.conv;
import std.string;
import util.lua;
import luageom;
import luagpath;
import luasurface;

void main()
{
    writeln("Begin demonstration of Lua connection to Surface objects.");
    auto L = luaL_newstate();
    luaL_openlibs(L);
    registerVector3(L);
    registerPaths(L);
    registerSurfaces(L);
    string test_code = `
print("Construct from edges")
a = Vector3:new{0.0, 0.0}
b = Vector3:new{0.0, 1.0}
c = Vector3:new{1.0, 0.0}
d = Vector3:new{1.0, 1.0}
surf = CoonsPatch:new{north=Line:new{b, d}, east=Line:new{c, d},
                      south=Line:new{a, c}, west=Line:new{a, b}}
print("CoonsPatch representation: ", surf)
print("surf(0.5,0.5)= ", surf(0.5, 0.5))
--
print("Try construction using corners")
surf2 = CoonsPatch:new{p00=a, p01=b, p11=c, p10=d}
p = surf2:eval(0.5, 0.5)
print("same point p= ", p)
--
print("AO patch")
p00 = Vector3:new{0.0, 0.1, 3.0}
p10 = Vector3:new{1.0, 0.4, 3.0}
p11 = Vector3:new{1.0, 1.1, 3.0}
p01 = Vector3:new{0.0, 1.1, 3.0}
my_aopatch = AOPatch:new{p00=p00, p10=p10, p11=p11, p01=p01}
p = my_aopatch(0.1, 0.1);
print("my_aopatch(0.1, 0.1)= ", p)
--
print("SubRangedSurface")
srs = SubRangedSurface:new{my_aopatch, r0=0.0, r1=0.5, s0=0.0, s1=0.5}
print("srs(0.2,0.2)=", srs(0.2,0.2))
--
print("ChannelPatch")
cA = Line:new{Vector3:new{0.0,0.0}, Vector3:new{1.0,0.0}}
cB = Line:new{Vector3:new{0.0,0.25}, Vector3:new{1.0,1.0}}
chanp = ChannelPatch:new{south=cA, north=cB}
print("chanp= ", chanp)
print("chanp(0.5,0.5)= ", chanp(0.5, 0.5))
--
print("Utility functions")
print("isSurface(my_aopatch)= ", isSurface(my_aopatch))
print("isSurface(surf2)= ", isSurface(surf2));
print("isSurface(a)= ", isSurface(a));
north = Line:new{b, d}
east = Line:new{c, d}
south = Line:new{a, c}
west = Line:new{a, b}
surf3 = makePatch{north, east, south, west, gridType="ao"}
print("surf3= ", surf3)
print("Done luasurface_demo.")
    `;
    if ( luaL_dostring(L, toStringz(test_code)) != 0 ) {
	writeln("There was a problem interpreting the test code.");
	writeln(to!string(lua_tostring(L, -1)));
    }
    writeln("Done with luageom_demo.");
}

    
