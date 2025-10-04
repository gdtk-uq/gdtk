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
    string test_code = `
print("1D grid")
a = Vector3:new{x=0.0, y=0.0}
b = Vector3:new{x=0.0, y=1.0}
grid1D = StructuredGrid:new{path=Line:new{p0=a, p1=b}, niv=11}
print("somewhere in the middle=", grid1D:get_vtx(5))
--
print("2D grid")
c = Vector3:new{x=1.0, y=0.0}
d = Vector3:new{x=1.0, y=1.0}
surf = CoonsPatch:new{north=Line:new{p0=b, p1=d}, east=Line:new{p0=c, p1=d},
                      south=Line:new{p0=a, p1=c}, west=Line:new{p0=a, p1=b}}
print("CoonsPatch representation: ", surf)
print("surf(0.5,0.5)= ", surf(0.5, 0.5))
--
myrf = RobertsFunction:new{end0=true, end1=false, beta=1.01}
grid = StructuredGrid:new{psurface=surf, niv=10, njv=20, interpolation="linear"}
grid:set_tags{north="my-special-tag"}
print("type=", grid:get_type(), "dimensions=", grid:get_dimensions())
ni = grid:get_niv()
nj = grid:get_njv()
print("grid size=", ni, nj)
print("upper-right corner=", grid:get_vtx(ni-1,nj-1))
grid:write_to_stl_file("test_sgrid.stl", 25.0)
--
print("SlabGrid")
grid3 = grid:makeSlabGrid{dz=0.2, symmetric=true, label="slab"}
grid3:write_to_vtk_file("test_grid3-slab.vtk");
print("WedgeGrid")
grid4 = grid:makeWedgeGrid{dtheta=0.2, symmetric=true, label="wedge"}
grid4:write_to_vtk_file("test_grid4-wedge.vtk");
--
print("MeshPatch surface")
surf2 = MeshPatch:new{sgrid=grid}
print("mid-point on MeshPatch surface p= ", surf2(0.5, 0.5))
--

print("StructuredGrid with internal control.")
p00 = Vector3:new{x=0.0, y=0.1}
p10 = Vector3:new{x=1.0, y=0.4}
p11 = Vector3:new{x=1.0, y=1.1}
p01 = Vector3:new{x=0.0, y=1.1}
myPatch = AOPatch:new{p00=p00, p10=p10, p11=p11, p01=p01}
r_grid = {{0.0, 1/3, 2/3, 1.0},
          {0.0, 1/3-0.3, 2/3+0.3, 1.0},
          {0.0, 1/3-0.3, 2/3+0.3, 1.0},
          {0.0, 1/3, 2/3, 1.0}}
s_grid = {{0.0, 0.0, 0.0, 0.0},
          {1.0/3, 1.0/3-0.3, 1.0/3-0.3, 1.0/3},
          {2.0/3, 2.0/3+0.3, 2.0/3+0.3, 2.0/3},
          {1.0, 1.0, 1.0, 1.0}}
gridB = StructuredGrid:new{psurface=myPatch, niv=11, njv=21, r_grid=r_grid, s_grid=s_grid}
gridB:write_to_vtk_file("test_grid_b-2D-lua.vtk")
print("initial badness=", gridB:measure_of_badness())
r_grid, s_grid = gridB:determine_rs_grids(myPatch)
print("final badness=", gridB:measure_of_badness())
print("r_grid=[[", r_grid[1][1], r_grid[1][2], r_grid[1][3], r_grid[1][4], "]")
print("        [", r_grid[2][1], r_grid[2][2], r_grid[2][3], r_grid[2][4], "]")
print("        [", r_grid[3][1], r_grid[3][2], r_grid[3][3], r_grid[3][4], "]")
print("        [", r_grid[4][1], r_grid[4][2], r_grid[4][3], r_grid[4][4], "]]")
print("s_grid=[[", s_grid[1][1], s_grid[1][2], s_grid[1][3], s_grid[1][4], "]")
print("        [", s_grid[2][1], s_grid[2][2], s_grid[2][3], s_grid[2][4], "]")
print("        [", s_grid[3][1], s_grid[3][2], s_grid[3][3], s_grid[3][4], "]")
print("        [", s_grid[4][1], s_grid[4][2], s_grid[4][3], s_grid[4][4], "]]")
gridB:write_to_vtk_file("test_grid_b_better-2D.vtk");
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
