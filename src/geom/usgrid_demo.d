// usgrid_demo.d

import std.stdio;
import std.conv;
import geom;
import gpath;
import surface;
import volume;
import univariatefunctions;
import sgrid;
import usgrid;

void main()
{
    writeln("Begin unstructured-grid demo...");
    auto p00 = Vector3(0.0, 0.1);
    auto p10 = Vector3(1.0, 0.4);
    auto p11 = Vector3(1.0, 1.1);
    auto p01 = Vector3(0.0, 1.1);
    auto my_patch = new AOPatch(p00, p10, p11, p01);
    auto cf = [new LinearFunction(), new LinearFunction(), 
	       new LinearFunction(), new LinearFunction()];
    auto my_grid = new StructuredGrid(my_patch, 11, 21, cf);
    writeln("grid point 5 5 at x=", my_grid[5,5].x, " y=", my_grid[5,5].y);
    auto usg = new UnstructuredGrid(my_grid, 2);
    usg.write_to_vtk_file("test_grid.vtk");
    usg.write_to_gzip_file("test_grid.gz");

    writeln("3D grid");
    Vector3[8] p;
    p[0] = Vector3(0.0, 0.1, 0.0);
    p[1] = Vector3(1.0, 0.1, 0.0);
    p[2] = Vector3(1.0, 1.1, 0.0);
    p[3] = Vector3(0.0, 1.1, 0.0);
    //
    p[4] = Vector3(0.0, 0.1, 3.0);
    p[5] = Vector3(1.0, 0.1, 3.0);
    p[6] = Vector3(1.0, 1.1, 3.0);
    p[7] = Vector3(0.0, 1.1, 3.0);
    //
    auto simple_box = new TFIVolume(p);
    auto my_3Dgrid = new StructuredGrid(simple_box, 11, 21, 11, cf);
    writeln("grid point 5 5 5 at p=", my_3Dgrid[5,5,5]);
    auto usg3D = new UnstructuredGrid(my_3Dgrid, 3);
    usg3D.write_to_vtk_file("test_3Dgrid.vtk");

    writeln("Done.");
}
