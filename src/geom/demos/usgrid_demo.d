// usgrid_demo.d

import std.stdio;
import std.conv;
import geom;

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
    auto usg = new UnstructuredGrid(my_grid);
    usg.write_to_vtk_file("test_grid.vtk");
    usg.write_to_stl_file("test_grid.stl", 25.4); // inches to mm scale

    usg.write_to_gzip_file("test_grid.gz");
    auto usg2 = new UnstructuredGrid("test_grid.gz", "gziptext", true);
    usg2.write_to_vtk_file("test_grid_2.vtk");

    usg.write_to_raw_binary_file("test_grid.bin");
    auto usg2b = new UnstructuredGrid("test_grid.bin", "rawbinary", true);
    usg2b.write_to_vtk_file("test_grid_2b.vtk");

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
    writeln("grid point 5 5 5 at p=", *my_3Dgrid[5,5,5]);
    auto usg3D = new UnstructuredGrid(my_3Dgrid);
    usg3D.sort_cells_into_bins(10, 10, 10);
    Vector3 my_point = 0.5*(p[0] + p[7]);
    size_t cell_indx = 0; bool found = false;
    usg3D.find_enclosing_cell(my_point, cell_indx, found);
    writeln("Search for cell enclosing my_point= ", my_point);
    if (found) {
        writeln("    cell found, index= ", cell_indx);
        writeln("    cell barycentre= ", usg3D.cell_barycentre(cell_indx));
    } else {
        writeln("    cell not found");
    }
    usg3D.write_to_vtk_file("test_3Dgrid.vtk");
    usg3D.write_to_gzip_file("test_3Dgrid.gz");
    auto usg3 = new UnstructuredGrid("test_3Dgrid.gz", "gziptext", true);
    usg3.write_to_vtk_file("test_3Dgrid_2.vtk");
    usg3.write_to_vtk_file("test_3Dgrid_2.su2");

    writeln("su2 2D grid -- triangles");
    auto su2grid = new UnstructuredGrid("sample-data/square-mesh.su2", "su2text", true);
    su2grid.write_to_vtk_file("test_su2-square-mesh.vtk");
    su2grid.write_to_su2_file("test_su2-square-mesh.su2");
    writeln("su2 2D grid -- quadrangles");
    auto su2grid2 = new UnstructuredGrid("sample-data/square-mesh-quads.su2", "su2text", true);
    su2grid2.write_to_vtk_file("test_su2-square-mesh-quads.vtk");
    su2grid2.write_to_su2_file("test_su2-square-mesh-quads.su2");
    writeln("su2 3D grid -- hexagons");
    auto su2grid3 = new UnstructuredGrid("sample-data/cube-mesh-hex.su2", "su2text", true);
    su2grid3.write_to_vtk_file("test_su2-cube-mesh-hex.vtk");
    su2grid3.write_to_su2_file("test_su2-cube-mesh-hex.su2");
    su2grid3.write_openFoam_polyMesh("test_openFoam");

    writeln("Try joining grids");
    auto grid_a = new StructuredGrid(simple_box, 3, 3, 4, cf);
    auto usg3D_a = new UnstructuredGrid(grid_a, "usga");
    Vector3[8] pb;
    pb[0] = Vector3(1.0, 0.1, 0.0);
    pb[1] = Vector3(2.0, 0.1, 0.0);
    pb[2] = Vector3(2.0, 1.1, 0.0);
    pb[3] = Vector3(1.0, 1.1, 0.0);
    //
    pb[4] = Vector3(1.0, 0.1, 3.0);
    pb[5] = Vector3(2.0, 0.1, 3.0);
    pb[6] = Vector3(2.0, 1.1, 3.0);
    pb[7] = Vector3(1.0, 1.1, 3.0);
    //
    auto simple_box_b = new TFIVolume(pb);
    auto grid_b = new StructuredGrid(simple_box_b, 3, 3, 4, cf);
    writeln("grid point 1 1 2 at p=", *my_3Dgrid[1,1,2]);
    auto usg3D_b = new UnstructuredGrid(grid_b, "usgb");
    usg3D_a.writeStats();
    usg3D_b.writeStats();
    usg3D_a.joinGrid(usg3D_b);
    usg3D_a.writeStats();
    usg3D_a.write_to_vtk_file("test-join-grid.vtk");

    writeln("Done.");
}
