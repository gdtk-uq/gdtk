/**
 * luasgrid_demo.d
 * Demonstrate the wrapped StructuredGrid objects.
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
import luavolume;
import luaunifunction;
import luasgrid;

void main()
{
    writeln("Begin demonstration of Lua connection to StructuredGrid objects.");
    auto L = luaL_newstate();
    luaL_openlibs(L);
    registerVector3(L);
    registerPaths(L);
    registerSurfaces(L);
    registerVolumes(L);
    registerUnivariateFunctions(L);
    registerStructuredGrid(L);
    string test_code = `
print("2D grid")
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=0.0, y=1.0}
c = Vector3:new{x=1.0, y=0.0}
d = Vector3:new{x=1.0, y=1.0}
surf = CoonsPatch:new{north=Line:new{p0=b, p1=d}, east=Line:new{p0=c, p1=d},
                      south=Line:new{p0=a, p1=c}, west=Line:new{p0=a, p1=b}}
print("CoonsPatch representation: ", surf)
print("surf(0.5,0.5)= ", surf(0.5, 0.5))
--
myrf = RobertsFunction:new{end0=true, end1=false, beta=1.01}
grid = StructuredGrid:new{psurface=surf, niv=10, njv=20} 
ni = grid:get_niv()
nj = grid:get_njv()
print("grid size=", ni, nj)
print("upper-right corner=", grid:get_vtx(ni-1,nj-1))
--
print("MeshPatch surface")
surf2 = MeshPatch:new{sgrid=grid}
print("mid-point on MeshPatch surface p= ", surf2(0.5, 0.5))
--
print("3D grid")
pArray = {Vector3:new{x=0.0, y=0.1, z=0.0}, Vector3:new{x=1.0, y=0.1, z=0.0},
          Vector3:new{x=1.0, y=1.1, z=0.0}, Vector3:new{x=0.0, y=1.1, z=0.0},
          Vector3:new{x=0.0, y=0.1, z=3.0}, Vector3:new{x=1.0, y=0.1, z=3.0},
          Vector3:new{x=1.0, y=1.1, z=3.0}, Vector3:new{x=0.0, y=1.1, z=3.0}}
volume = TFIVolume:new{vertices=pArray}
grid3D = StructuredGrid:new{pvolume=volume, niv=11, njv=21, nkv=11}
print("somewhere in the middle=", grid3D:get_vtx(5,10,5))
subgrid3D = grid3D:subgrid(5,3,10,3,5,3)
print("same point in subgrid=", subgrid3D:get_vtx(0,0,0))
--
print("Try 3D grid from user-defined Lua function")
function myLuaFunction(r, s, t)
   -- Simple cube
   return {x=r, y=s, z=t}
end
myFnVol = LuaFnVolume:new{luaFnName="myLuaFunction"}
grid3DLF = StructuredGrid:new{pvolume=myFnVol, niv=11, njv=21, nkv=11}
print("somewhere in the middle=", grid3DLF:get_vtx(5,10,5))
--
print("Try Gridpro import")
grids = importGridproGrid("../../examples/eilmer/3D/gridpro-import/blk.tmp", 0.001)
print("no. of grids read= ", #grids)
print("size of grid 1= ", grids[1]:get_niv(), grids[1]:get_njv())
    `;
    if ( luaL_dostring(L, toStringz(test_code)) != 0 ) {
	writeln("There was a problem interpreting the test code.");
	writeln(to!string(lua_tostring(L, -1)));
    }
    writeln("Done with luasgrid_demo.");
}

    
