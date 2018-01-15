/**
 * Authors: Rowan G. and Peter J.
 * Date: 2018-01-15
 *
 */

import std.stdio;
import std.math; 

import geom;

void main()
{
    //writeln("Demo of gpath_utils functions.");

    //writeln("Using optimiseBezierPoints.");
    Bezier bezier = optimiseBezierPoints("sample-data/cowl-spline-pts.dat", 7);
    
    //writeln("Optimised control points are:");
    foreach (p; bezier.B) {
        writeln(p.x, " ", p.y);
    }
}
