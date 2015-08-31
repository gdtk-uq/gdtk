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
a = Vector3:new{0.0, 0.0}
b = Vector3:new{0.0, 1.0}
c = Vector3:new{1.0, 0.0}
d = Vector3:new{1.0, 1.0}
surf = CoonsPatch:new{north=Line:new{b, d}, east=Line:new{c, d},
                      south=Line:new{a, c}, west=Line:new{a, b}}
print("CoonsPatch representation: ", surf)
myrf = RobertsFunction:new{end0=true, end1=false, beta=1.01}
grid = StructuredGrid:new{psurface=surf, niv=10, njv=20} 
ni = grid:get_niv()
nj = grid:get_njv()
print("grid size=", ni, nj)
print("upper-right corner=", grid:get_vtx(ni-1,nj-1))
print("3D grid")
pArray = {Vector3:new{0.0, 0.1, 0.0}, Vector3:new{1.0, 0.1, 0.0},
          Vector3:new{1.0, 1.1, 0.0}, Vector3:new{0.0, 1.1, 0.0},
          Vector3:new{0.0, 0.1, 3.0}, Vector3:new{1.0, 0.1, 3.0},
          Vector3:new{1.0, 1.1, 3.0}, Vector3:new{0.0, 1.1, 3.0}}
volume = TFIVolume:new{vertices=pArray}
grid3D = StructuredGrid:new{pvolume=volume, niv=11, njv=21, nkv=11}
print("somewhere in the middle=", grid3D:get_vtx(5,10,5))
subgrid3D = grid3D:subgrid(5,3,10,3,5,3)
print("same point in subgrid=", subgrid3D:get_vtx(0,0,0))
--
print("Try Gridpro import")
grids = importGridproGrid("../../examples/eilmer3/3D/gridpro-import/blk.tmp", 0.001)
print("no. of grids read= ", #grids)
print("size of grid 1= ", grids[1]:get_niv(), grids[1]:get_njv())
    `;
    if ( luaL_dostring(L, toStringz(test_code)) != 0 ) {
	writeln("There was a problem interpreting the test code.");
	writeln(to!string(lua_tostring(L, -1)));
    }
    writeln("Done with luasgrid_demo.");
}

    
