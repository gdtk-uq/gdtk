/**
 * luavolume_demo.d
 * Demonstrate the wrapped ParametricVolume objects.
 *
 * Authors: Peter J. and Rowan G.
 * Version: 2015-04-07
 */

import std.stdio;
import std.conv;
import std.string;
import util.lua;
import geom.luawrap;

void main()
{
    writeln("Begin demonstration of Lua connection to Volume objects.");
    auto L = luaL_newstate();
    luaL_openlibs(L);
    registerVector3(L);
    registerPaths(L);
    registerSurfaces(L);
    registerVolumes(L);
    string test_code = `
p000 = Vector3:new{x=0.0, y=0.1, z=0.0}
p100 = Vector3:new{x=1.0, y=0.1, z=0.0}
p110 = Vector3:new{x=1.0, y=1.1, z=0.0}
p010 = Vector3:new{x=0.0, y=1.1, z=0.0}
p001 = Vector3:new{x=0.0, y=0.1, z=3.0}
p101 = Vector3:new{x=1.0, y=0.1, z=3.0}
p111 = Vector3:new{x=1.0, y=1.1, z=3.0}
p011 = Vector3:new{x=0.0, y=1.1, z=3.0}
-- Beware of Lua indexing for tables starting at 1.
-- By providing a list literal, we side-step the issue.
my_volume = TFIVolume:new{vertices={p000,p100,p110,p010,p001,p101,p111,p011}}
print("my_volume=", my_volume)
p = my_volume(0.1, 0.1, 0.5)
print("my_volume(0.1, 0.1, 0.5)= ", p)
--
print("SweptSurfaceVolume demo")
myface0123 = CoonsPatch:new{p00=p000, p10=p100, p11=p110, p01=p010}
myedge04 = Line:new{p0=p000, p1=p001}
ssv = SweptSurfaceVolume:new{face0123=myface0123, edge04=myedge04}
print("ssv=", ssv)
print("ssv(0.1, 0.1, 0.5)= ", ssv(0.1, 0.1, 0.5))
--
print("SlabVolume demo")
mydz = p001 - p000
slabv = SlabVolume:new{face0123=myface0123, dz=mydz}
print("slabv=", slabv)
print("slabv(0.1, 0.1, 0.5)= ", slabv(0.1, 0.1, 0.5))
--
print("WedgeVolume demo")
wedgev = WedgeVolume:new{face0123=myface0123, dtheta=0.1}
print("wedgev=", wedgev)
print("wedgev(0.1, 0.1, 0.5)= ", wedgev(0.1, 0.1, 0.5))
--
print("TwoSurfaceVolume demo")
myface4567 = CoonsPatch:new{p00=p001, p10=p101, p11=p111, p01=p011}
tsv = TwoSurfaceVolume:new{face0=myface0123, face1=myface4567, ruled_direction="k"}
print("tsv=", tsv)
print("tsv(0.1, 0.1, 0.5)= ", tsv(0.1, 0.1, 0.5))
--
print("LuaFnVolume")
function myLuaFunction(r, s, t)
   -- Simple cube
   return {x=r, y=s, z=t}
end
myFnVol = LuaFnVolume:new{luaFnName="myLuaFunction"}
print("myFnVol= ", myFnVol)
print("myFnVol(0.3, 0.4, 0.5)= ", myFnVol(0.3, 0.4, 0.5))
--
print("SubRangedVolume demo")
srv = SubRangedVolume:new{underlying_pvolume=my_volume,
                          r0=0.0,r1=0.5,s0=0.0,s1=0.5,t0=0.0,t1=0.5}
print("srv(0.2,0.2,1.0)=", srv(0.2,0.2,1.0))
--
print("isVolume(my_volume)= ", isVolume(my_volume))
print("Done luavolume_demo")
    `;
    if ( luaL_dostring(L, toStringz(test_code)) != 0 ) {
        writeln("There was a problem interpreting the test code.");
        writeln(to!string(lua_tostring(L, -1)));
    }
    writeln("Done with luavolume_demo.");
}
