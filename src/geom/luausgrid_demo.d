/**
 * luausgrid_demo.d
 * Demonstrate the wrapped UnstructuredGrid objects.
 *
 * Authors: Peter J. and Rowan G.
 * Version: 2015-11-06
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
import luausgrid;

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
    registerUnstructuredGrid(L);
    string test_code = `
print("2D grid")
a = Vector3:new{0.0, 0.0}
b = Vector3:new{0.0, 1.0}
c = Vector3:new{1.0, 0.2}
d = Vector3:new{1.0, 1.0}
surf = CoonsPatch:new{north=Line:new{b, d}, east=Line:new{c, d},
                      south=Line:new{a, c}, west=Line:new{a, b}}
print("CoonsPatch representation: ", surf)
my_grid = StructuredGrid:new{psurface=surf, niv=10, njv=20}
my_usg = UnstructuredGrid:new{sgrid=my_grid} 
my_usg:write_to_vtk_file("test_grid.vtk")
my_usg:write_to_gzip_file("test_grid.gz")
my_usg_2 = UnstructuredGrid:new{filename="test_grid.gz"}
my_usg_2:write_to_vtk_file("test_grid_2.vtk")

print("3D grid")
pArray = {Vector3:new{0.0, 0.1, 0.0}, Vector3:new{1.0, 0.1, 0.0},
          Vector3:new{1.0, 1.1, 0.0}, Vector3:new{0.0, 1.1, 0.0},
          Vector3:new{0.0, 0.1, 3.0}, Vector3:new{1.0, 0.1, 3.0},
          Vector3:new{1.0, 1.1, 3.0}, Vector3:new{0.0, 1.1, 3.0}}
volume = TFIVolume:new{vertices=pArray}
my_grid3D = StructuredGrid:new{pvolume=volume, niv=11, njv=21, nkv=11}
my_usg3D = UnstructuredGrid:new{sgrid=my_grid3D} 
my_usg3D:write_to_vtk_file("test_grid3D.vtk")
    `;
    if ( luaL_dostring(L, toStringz(test_code)) != 0 ) {
	writeln("There was a problem interpreting the test code.");
	writeln(to!string(lua_tostring(L, -1)));
    }
    writeln("Done with luausgrid_demo.");
}

    
