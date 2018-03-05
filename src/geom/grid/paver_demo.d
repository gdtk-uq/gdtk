// Paver_demo.d
/* 
15/4/2016
Heather Muir - Undergraduate Thesis
*/

import std.stdio;
import std.conv;
import geom;


Vector3[] construct_boundary(Path[] edges, size_t[] n)
{
    assert(edges.length == n.length);
    Vector3[] boundary;
    size_t i = 0;
    double[] t_list;
    double t;
    double t_incr;
    foreach(edge;edges){
        t = 0;
        t_incr = 1.0/(n[i]-1);
        t_list = [];
        for(int j; j<(n[i]-1); ++j){
            t_list ~= t;
            t += t_incr;
        }
        foreach(val; t_list){
            boundary ~= edge.opCall(val);
        }
        ++i;
    }
    return boundary;
} 

BoundaryFaceSet[] construct_boundary_faces(string[] BC_list, size_t[] n)
{
    assert(BC_list.length == n.length);
    BoundaryFaceSet[] boundaries = [];
    size_t index = 0;
    size_t list_index = 0;
    foreach(BC; BC_list){
        boundaries ~= new BoundaryFaceSet(BC);
        size_t[] face_ID_list = [];
        int[] outsign_list = [];
        for(int i; i<n[list_index]; ++i){
            face_ID_list ~= index;
            outsign_list ~= -1;
            ++index;
        }
        boundaries[$-1].face_id_list = face_ID_list;
        boundaries[$-1].outsign_list = outsign_list;
        ++list_index;
    }
    return boundaries;
}


void main()
{
    /*notes:
        - define boundary by points connected by paths from gpath
        - boundary must be anti-clockwise
        - boundary must contain an even number of nodes
        - try to avoid large jumps in cell size
    /*

      P3--------<--------P2
       |                 |
       |                 |
       v   (clockwise)   ^
       |                 |
       |                 |
      PO-------->--------P1 

    */

    //points:
    Vector3 P0 = Vector3(-0.5,0);
    Vector3 P1 = Vector3(0,0);
    Vector3 P2 = Vector3(0.5, 1.0);
    Vector3 P3 = Vector3(3, 3);
    Vector3 P4 = Vector3(3, 5.5);
    Vector3 P5 = Vector3(0, 1.5);

    Vector3 M0 = Vector3(0.108, 0.54);
    Vector3 M1 = Vector3(-0.3, 1.04);


    //paths:                            node counts:            boundary conditions:
    auto P0P1 = new Line(P0, P1);       size_t n0 = 10;         string BC0 = "slip-wall";
    auto P1P2 = new Arc3(P1, M0, P2);   size_t n1 = 30;         string BC1 = "slip-wall";
    auto P2P3 = new Line(P2, P3);       size_t n2 = 50;         string BC2 = "slip-wall";
    auto P3P4 = new Line(P3, P4);       size_t n3 = 30;         string BC3 = "outflow-boundary";
    auto P4P5 = new Line(P4, P5);       size_t n4 = 70;         string BC4 = "inflow-boundary";
    auto P5P0 = new Arc3(P5, M1, P0);   size_t n5 = 30;         string BC5 = "inflow-boundary";


    // construct the list of boundary points:
    Vector3[] boundary_points = construct_boundary([P0P1, P1P2, P2P3, P3P4, P4P5, P5P0],
                                            [n0, n1, n2, n3, n4, n5]);

    // construct the boundary condition sets:
    BoundaryFaceSet[] boundaries = construct_boundary_faces([BC0, BC1, BC2, BC3, BC4, BC5], 
                                                        [n0, n1, n2, n3, n4, n5]);


    // use the old paver constructor to make an unstructured grid:
    auto grid = new UnstructuredGrid(boundary_points, boundaries, "PavedGrid1");
    grid.write_to_vtk_file("paved_grid.vtk");
} // end main()
    



