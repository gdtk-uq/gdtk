// paver2_demo.d
// 
// 15/4/2016 Heather Muir - Undergraduate Thesis
// 2018-03-05, PJ integrated into geometry package.

import std.stdio;
import std.conv;
import geom;


void main()
{
    /*notes:
        - define boundary by points connected by paths from gpath
        - boundary must be anti-clockwise
        - boundary must contain an even number of segments [TODO] check assertion
        - try to avoid large jumps in cell size
      P3--------<--------P2
       |                 |
       |                 |
       v   (clockwise)   ^
       |                 |
       |                 |
      PO-------->--------P1 
    */

    Vector3 P0 = Vector3(-0.5,0);
    Vector3 P1 = Vector3(0,0);
    Vector3 P2 = Vector3(0.5, 1.0);
    Vector3 P3 = Vector3(3, 3);
    Vector3 P4 = Vector3(3, 5.5);
    Vector3 P5 = Vector3(0, 1.5);
    Vector3 M0 = Vector3(0.108, 0.54);
    Vector3 M1 = Vector3(-0.3, 1.04);

    auto pths = [new Line(P0, P1), new Arc3(P1, M0, P2), new Line(P2, P3),
                 new Line(P3, P4), new Line(P4, P5), new Arc3(P5, M1, P0)];
    size_t[] nn = [11, 31, 51, 31, 71, 31];

    // Discretize the bounding paths.
    Vector3[][] bndryArray;
    auto cf = new LinearFunction();
    foreach (i, pth; pths) { bndryArray ~= discretize_path(pth, nn[i], cf); }
    // use the new paver constructor to make an unstructured grid.
    auto grid2 = new UnstructuredGrid(bndryArray, "PavedGrid2");
    grid2.write_to_vtk_file("paved_grid2.vtk");
} // end main()
    



