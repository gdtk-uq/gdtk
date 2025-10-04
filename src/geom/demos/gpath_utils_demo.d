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
    writeln("Demo of gpath_utils functions.");
    writeln("Using optimiseBezierPoints.");
    Bezier b_init_guess;
    double[] ts;
    bool success;
    Bezier bezier = optimiseBezierPoints("sample-data/cowl-spline-pts.dat", 7, b_init_guess, ts, success);

    writeln("Optimised control points are:");
    foreach (p; bezier.B) {
        writeln(p.x, " ", p.y);
    }
    foreach (t; ts) {
        writeln(t);
    }
}
