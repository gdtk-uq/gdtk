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
import geom.luawrap;

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
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=0.0, y=1.0}
c = Vector3:new{x=1.0, y=0.2}
d = Vector3:new{x=1.0, y=1.0}
surf = CoonsPatch:new{north=Line:new{p0=b, p1=d}, east=Line:new{p0=c, p1=d},
                      south=Line:new{p0=a, p1=c}, west=Line:new{p0=a, p1=b}}
print("CoonsPatch representation: ", surf)
my_grid = StructuredGrid:new{psurface=surf, niv=10, njv=20}
my_grid:set_tags{north="my-special-tag"}
my_usg = UnstructuredGrid:new{sgrid=my_grid}
print("type=", my_usg:get_type(), "dimensions=", my_usg:get_dimensions())
my_usg:write_to_vtk_file("test_grid.vtk")
my_usg:write_to_stl_file("test_usgrid.stl", 25.0)
my_usg:write_to_gzip_file("test_grid.gz")
my_usg:write_to_su2_file("test_grid.su2")
my_usg_2 = UnstructuredGrid:new{filename="test_grid.gz"}
my_usg_2:write_to_vtk_file("test_grid_2.vtk")

print("3D grid")
pArray = {Vector3:new{x=0.0, y=0.1, z=0.0}, Vector3:new{x=1.0, y=0.1, z=0.0},
          Vector3:new{x=1.0, y=1.1, z=0.0}, Vector3:new{x=0.0, y=1.1, z=0.0},
          Vector3:new{x=0.0, y=0.1, z=3.0}, Vector3:new{x=1.0, y=0.1, z=3.0},
          Vector3:new{x=1.0, y=1.1, z=3.0}, Vector3:new{x=0.0, y=1.1, z=3.0}}
volume = TFIVolume:new{vertices=pArray}
my_grid3D = StructuredGrid:new{pvolume=volume, niv=11, njv=21, nkv=11}
my_usg3D = UnstructuredGrid:new{sgrid=my_grid3D}
my_usg3D:write_to_vtk_file("test_grid3D.vtk")
for i=0,my_usg3D:get_nboundaries()-1 do
   print(string.format("boundaryset_tag[%d]=%s", i, my_usg3D:get_boundaryset_tag(i)))
end

print("Try joining grids")
grid_a = StructuredGrid:new{pvolume=volume, niv=3, njv=3, nkv=4}
usg3D_a = UnstructuredGrid:new{sgrid=grid_a, label="usga"}
pArrayb = {Vector3:new{x=1.0, y=0.1, z=0.0}, Vector3:new{x=2.0, y=0.1, z=0.0},
           Vector3:new{x=2.0, y=1.1, z=0.0}, Vector3:new{x=1.0, y=1.1, z=0.0},
           Vector3:new{x=1.0, y=0.1, z=3.0}, Vector3:new{x=2.0, y=0.1, z=3.0},
           Vector3:new{x=2.0, y=1.1, z=3.0}, Vector3:new{x=1.0, y=1.1, z=3.0}}

simple_box_b = TFIVolume:new{vertices=pArrayb}
grid_b = StructuredGrid:new{pvolume=simple_box_b, niv=3, njv=3, nkv=4}
usg3D_b = UnstructuredGrid:new{sgrid=grid_b, label="usgb"}
usg3D_a:writeStats()
usg3D_b:writeStats()
usg3D_a:joinGrid(usg3D_b)
usg3D_a:writeStats()
usg3D_a:write_to_vtk_file("test-join-grid_lua.vtk")
usg3D_a:writeOpenFoamPolyMesh("test_openFoam")
    `;
    if ( luaL_dostring(L, toStringz(test_code)) != 0 ) {
        writeln("There was a problem interpreting the test code.");
        writeln(to!string(lua_tostring(L, -1)));
    }
    writeln("Done with luausgrid_demo.");
}
