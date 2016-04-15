// Paver_demo.d
/* 
15/4/2016
Heather Muir - Undergraduate Thesis
*/

import std.stdio;
import std.conv;
import geom;
import gpath;
import surface;
import volume;
import univariatefunctions;
import sgrid;
import usgrid;

import paver;

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


void main()
{
    /*notes:
	- define boundary by points connected by paths from gpath
	- boundary must be anti-clockwise
	- boundary contain an even number of points
	- try to avoid large differences in adjacent edge cell counts
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
    Vector3 P0 = Vector3(0,0);
    Vector3 P1 = Vector3(5,0);
    Vector3 P2 = Vector3(5,5);
    Vector3 P3 = Vector3(0,5);

    //paths:
    auto P0P1 = new Line(P0, P1); size_t n0 = 30;
    auto P1P2 = new Line(P1, P2); size_t n1 = 30;
    auto P2P3 = new Line(P2, P3); size_t n2 = 30;
    auto P3P0 = new Line(P3, P0); size_t n3 = 30;

    Vector3[] boundary_points = construct_boundary([P0P1, P1P2, P2P3, P3P0],
					    [n0, n1, n2, n3]);

    BoundaryFaceSet[] boundaries = [];

    auto grid = new UnstructuredGrid(boundary_points, boundaries, "PavedGrid1");
    grid.write_to_vtk_file("paved_grid.vtk");

}
    



