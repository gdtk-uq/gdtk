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
    Vector3 P0 = Vector3(-2,-1.5);
    Vector3 P1 = Vector3(3,-1.5);
    Vector3 P2 = Vector3(3,0);
    Vector3 P3 = Vector3(3,1.5);
    Vector3 P4 = Vector3(-2,1.5);
    Vector3 P5 = Vector3(1,0);
    Vector3 P6 = Vector3(0,0);

    Vector3 B1 = Vector3(0, 0.0361); 		Vector3 B1_ = Vector3(0, -0.0361);
    Vector3 B2 = Vector3(0.0994, 0.062);	Vector3 B2_ = Vector3(0.0994, -0.062);
    Vector3 B3 = Vector3(0.2494, 0.0771);	Vector3 B3_ = Vector3(0.2494, -0.0771);
    Vector3 B4 = Vector3(0.3976, 0.0620);	Vector3 B4_ = Vector3(0.3976, -0.0620);
    Vector3 B5 = Vector3(0.5518, 0.0241);	Vector3 B5_ = Vector3(0.5518, -0.0241);
    Vector3 B6 = Vector3(0.7018, 0.04096);	Vector3 B6_ = Vector3(0.7018, -0.04096);
    Vector3 B7 = Vector3(0.8542, 0.0241);	Vector3 B7_ = Vector3(0.8542, -0.0241);

    //paths:
    auto P0P1 = new Line(P0, P1); size_t n0 = 60;
    auto P1P2 = new Line(P1, P2); size_t n1 = 20;
    auto P2P5 = new Line(P2, P5); size_t n2 = 25;
    auto P5P6 = new Bezier([P5, B7_, B6_, B5_, B4_, B3_, B2_, B1_, P6]); size_t n3 = 30;
    auto P6P5 = new Bezier([P6, B1, B2, B3, B4, B5, B6, B7, P5]); size_t n4 = 30;
    auto P5P2 = new Line(P5, P2); size_t n5 = 25;
    auto P2P3 = new Line(P2, P3); size_t n6 = 20;
    auto P3P4 = new Line(P3, P4); size_t n7 = 60;
    auto P4P0 = new Line(P4, P0); size_t n8 = 41;

    Vector3[] boundary_points = construct_boundary([P0P1, P1P2, P2P5, P5P6, P6P5, P5P2, P2P3, P3P4, P4P0],
					    [n0, n1, n2, n3, n4, n5, n6, n7, n8]);

    BoundaryFaceSet[] boundaries = [];

    auto grid = new UnstructuredGrid(boundary_points, boundaries, "PavedGrid1");
    grid.write_to_vtk_file("paved_grid.vtk");

}
    



