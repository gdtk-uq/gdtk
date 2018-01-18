/**
 * Authors: Rowan G. and Peter J.
 * Date: 2018-01-14
 *
 */

module geom.gpath.gpath_utils;

import std.conv;
import std.stdio;
import std.string;
import std.range : iota;
import geom;
import nm.nelmin;

void readPointsFromFile(string fileName, ref Vector3[] points)
{
    auto f = File(fileName, "r");
    foreach (line; f.byLine) {
        auto tokens = line.strip().split();
        if (tokens.length == 0) continue; // ignore blank lines
        if (tokens[0] == "#") continue; // ignore comment lines
        double x = to!double(tokens[0]);
        double y = (tokens.length > 1) ? to!double(tokens[1]) : 0.0;
        double z = (tokens.length > 2) ? to!double(tokens[2]) : 0.0;
        points ~= Vector3(x, y, z);
    }
    f.close();
}

Bezier optimiseBezierPoints(string fileName, int nCtrlPts, ref double[] ts, int dim=2)
{
    // Read in points.
    Vector3[] points;
    readPointsFromFile(fileName, points);
    // call bezier optimiser
    Bezier myBez = optimiseBezierPoints(points, nCtrlPts, ts);
    return myBez;
}

Bezier optimiseBezierPoints(Vector3[] points, int nCtrlPts, ref double[] ts, int dim=2)
{
    double t;
    // Check that there are more data points than
    // desired control points
    if (points.length < nCtrlPts) {
        string errMsg = "Error in bestFitBezier: there are fewer data points than the desired number of control points.\n";
        errMsg ~= format("No. of data points: %d\n", points.length);
        errMsg ~= format("No. of desired control points: %d", nCtrlPts);
        throw new Error(errMsg);
    }
    ts.length = points.length;
    // ------------------------------------------------------------------
    // Establish an initial guess for the location of the control points.
    // ------------------------------------------------------------------
    // Distribute Bezier control points at roughly equal
    // spacing using the given data points.
    Vector3[] bezPts;
    // Add starting point as first control point.
    bezPts ~= points[0];
    // Distribute points in line joining between start and end points
    auto line = new Line(points[0], points[$-1]);
    double dt = 1.0/(nCtrlPts-1);
    foreach (i; 1 .. nCtrlPts-1) {
        bezPts ~= line(i*dt);
    }
    // Add end point as final contrl point.
    bezPts ~= points[$-1];
    Bezier myBez = new Bezier(bezPts);
    double[] d;
    if (dim == 2)
        d.length = 2*(nCtrlPts-2);
    else
        d.length = 3*(nCtrlPts-2);
    foreach (i; 1 .. nCtrlPts-1) {
        if (dim == 2) {
            d[2*i - 2] = myBez.B[i].x;
            d[2*i - 1] = myBez.B[i].y;
        }
        else {
            d[3*i - 3] = myBez.B[i].x;
            d[3*i - 2] = myBez.B[i].y;
            d[3*i - 1] = myBez.B[i].z;
        }
    }
    // --------------- Done establishing guess ---------------------------
    
    // --------------------------------------------------------------------
    // Build cost function to be minimised.
    // --------------------------------------------------------------------
    double fMin(double[] d)
    {
        // Adjust Bezier points based on supplied design values 'd'
        foreach (i; 1 .. nCtrlPts-1) {
            if (dim == 2) {
                myBez.B[i].refx = d[2*i - 2];
                myBez.B[i].refy = d[2*i - 1];
            }
            else {
                myBez.B[i].refx = d[3*i - 3];
                myBez.B[i].refy = d[3*i - 2];
                myBez.B[i].refz = d[3*i - 1];
            }
        }
        // Evaluate error between data points and Bezier curve
        double err = 0.0;
        foreach (p; points) {
            err += myBez.closestDistance(p, t);
        }
        return err;
    }
    // ------------- Done building cost function -----------------------

    // -----------------------------------------------------------------
    // Optimise placement of control points using Nelder-Mead minimiser
    // -----------------------------------------------------------------
    double f_min;
    int n_fe, n_restart;
    double[] dx;
    dx.length = d.length;
    // Make the initial perturbations 1/100th of the arc length.
    double dp = myBez.length()/100.0;
    dx[] = dp;
    auto success = minimize!fMin(d, f_min, n_fe, n_restart, dx, 1.0e-6, 10000);
    if (success) {
        foreach (i; 1 .. nCtrlPts-1) {
            if (dim == 2) {
                myBez.B[i].refx = d[2*i - 2];
                myBez.B[i].refy = d[2*i - 1];
            }
            else {
                myBez.B[i].refx = d[3*i - 3];
                myBez.B[i].refy = d[3*i - 2];
                myBez.B[i].refz = d[3*i - 1];
            }
        }
	// Populate ts vector with all the t-values associated with the data points.
	ts[0] = 0.0;
	foreach (i; 1 .. points.length-1) {
            myBez.closestDistance(points[i], t);
	    ts[i] = t;
        }
	ts[$-1] = 1.0;
        return myBez;
    }
    // Otherwise, we failed in the optimiser.
    string errMsg = "Error in bestFitBezier while using optimiser.\n";
    errMsg ~= "Optimiser stats:\n";
    errMsg ~= format("  number of function evaluations: %d\n", n_fe);
    errMsg ~= format("  number of restarts: %d\n", n_restart);
    errMsg ~= format("  minimum of function found: %12.6e\n", f_min);
    throw new Exception(errMsg);
    // ------------- Done calling optimiser -----------------------
}
