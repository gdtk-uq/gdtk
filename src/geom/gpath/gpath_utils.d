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
import ntypes.complex;
import nm.number;


Vector3[] discretize_path(const Path pth, size_t niv,
                          const(UnivariateFunction) clusterf)
{
    // First, set up clustered parameter values.
    double[] r = clusterf.distribute_parameter_values(niv);
    // Now, accumulate points, one at a time.
    Vector3[] points;
    foreach (i; 0 .. niv) { points ~= pth(r[i]); }
    return points;
} // end discretize_path()

void readPointsFromFile(string fileName, ref Vector3[] points)
{
    auto f = File(fileName, "r");
    foreach (line; f.byLine) {
        auto tokens = line.strip().split();
        if (tokens.length == 0) continue; // ignore blank lines
        if (tokens[0] == "#") continue; // ignore comment lines
        number x = to!double(tokens[0]);
        number y = (tokens.length > 1) ? to!double(tokens[1]) : 0.0;
        number z = (tokens.length > 2) ? to!double(tokens[2]) : 0.0;
        points ~= Vector3(x, y, z);
    }
    f.close();
}

Bezier optimiseBezierPoints(string fileName, int nCtrlPts, Bezier initGuess, ref double[] ts, ref bool success, int dim=2)
{
    // Read in points.
    Vector3[] points;
    readPointsFromFile(fileName, points);
    // call bezier optimiser
    Bezier myBez = optimiseBezierPoints(points, nCtrlPts, initGuess, ts, success);
    return myBez;
}

Bezier optimiseBezierPoints(Vector3[] points, int nCtrlPts, Bezier initGuess, ref double[] ts, ref bool success,
                            double tol=1.0e-06, int max_steps=10000, int dim=2)
{
    double t;
    // Check that there are more data points than
    // desired control points
    if (points.length < nCtrlPts) {
        string errMsg = "Error in optimiseBezierPoints: there are fewer data points than the desired number of control points.\n";
        errMsg ~= format("No. of data points: %d\n", points.length);
        errMsg ~= format("No. of desired control points: %d", nCtrlPts);
        throw new Exception(errMsg);
    }
    ts.length = points.length;
    // ------------------------------------------------------------------
    // Establish an initial guess for the location of the control points.
    // ------------------------------------------------------------------

    Bezier myBez;
    if (initGuess) {
        if (initGuess.B.length != nCtrlPts) {
            string errMsg = "Error in optimiseBezierPoints: the supplied initial guess Bezier does not have the same number of controls points\n";
            errMsg ~= "as the desired number of control points for the optimised Bezier.\n";
            errMsg ~= format("No. of control points desired: %d\n", nCtrlPts);
            errMsg ~= format("No. of control points in supplied initial guess: %d", initGuess.B.length);
            throw new Exception(errMsg);
        }
        myBez = new Bezier(initGuess);
    }
    else {
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
        myBez = new Bezier(bezPts);
    }
    number[] d;
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
    number funMin(number[] d)
    {
        // Adjust Bezier points based on supplied design values 'd'
        foreach (i; 1 .. nCtrlPts-1) {
            if (dim == 2) {
                myBez.B[i].x = d[2*i - 2];
                myBez.B[i].y = d[2*i - 1];
            }
            else {
                myBez.B[i].x = d[3*i - 3];
                myBez.B[i].y = d[3*i - 2];
                myBez.B[i].z = d[3*i - 1];
            }
        }
        // Evaluate error between data points and Bezier curve
        number err = 0.0;
        foreach (p; points) {
            err += myBez.closestDistance(p, t);
        }
        return err;
    }
    // ------------- Done building cost function -----------------------

    // -----------------------------------------------------------------
    // Optimise placement of control points using Nelder-Mead minimiser
    // -----------------------------------------------------------------
    number f_min;
    int n_fe, n_restart;
    number[] dx;
    dx.length = d.length;
    // Make the initial perturbations 1/100th of the arc length.
    number dp = myBez.length()/100.0;
    dx[] = dp;
    success = nm.nelmin.minimize!(funMin,number)(d, f_min, n_fe, n_restart, dx, tol, 1, max_steps);
    foreach (i; 1 .. nCtrlPts-1) {
        if (dim == 2) {
            myBez.B[i].x = d[2*i - 2];
            myBez.B[i].y = d[2*i - 1];
        }
        else {
            myBez.B[i].x = d[3*i - 3];
            myBez.B[i].y = d[3*i - 2];
            myBez.B[i].z = d[3*i - 1];
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
    // ------------- Done calling optimiser -----------------------
}

Bezier optimiseBezierPoints2(Vector3[] points, int nCtrlPts, Bezier initGuess, Vector3 startSlope, Vector3 endSlope,
                             ref double[] ts, ref bool success,
                             double tol=1.0e-06, int max_steps=10000, int dim=2)
{
    double t;
    // Check that there are more data points than
    // desired control points
    if (points.length < nCtrlPts) {
        string errMsg = "Error in optimiseBezierPoints: there are fewer data points than the desired number of control points.\n";
        errMsg ~= format("No. of data points: %d\n", points.length);
        errMsg ~= format("No. of desired control points: %d", nCtrlPts);
        throw new Exception(errMsg);
    }
    ts.length = points.length;
    // ------------------------------------------------------------------
    // Establish an initial guess for the location of the control points.
    // ------------------------------------------------------------------

    Bezier myBez;
    if (initGuess) {
        if (initGuess.B.length != nCtrlPts) {
            string errMsg = "Error in optimiseBezierPoints: the supplied initial guess Bezier does not have the same number of controls points\n";
            errMsg ~= "as the desired number of control points for the optimised Bezier.\n";
            errMsg ~= format("No. of control points desired: %d\n", nCtrlPts);
            errMsg ~= format("No. of control points in supplied initial guess: %d", initGuess.B.length);
            throw new Exception(errMsg);
        }
        myBez = new Bezier(initGuess);
    }
    else {
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
        // Fix first interior points
        auto scale0 = (bezPts[1].x - bezPts[0].x)/startSlope.x;
        bezPts[1] = bezPts[0] + scale0*startSlope;
        auto scale1 = (bezPts[nCtrlPts-1].x - bezPts[nCtrlPts-2].x)/endSlope.x;
        bezPts[nCtrlPts-2] = bezPts[nCtrlPts-1] - scale1*endSlope;

        myBez = new Bezier(bezPts);
    }
    number[] d;
    if (dim == 2)
        d.length = 2*(nCtrlPts-2) - 2;
    else
        d.length = 3*(nCtrlPts-2) - 4;
    foreach (i; 2 .. nCtrlPts-2) {
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
    d[$-2] = myBez.B[1].x;
    d[$-1] = myBez.B[nCtrlPts-2].x;
    // --------------- Done establishing guess ---------------------------

    // --------------------------------------------------------------------
    // Build cost function to be minimised.
    // --------------------------------------------------------------------
    number funMin(number[] d)
    {
        // Adjust Bezier points based on supplied design values 'd'
        foreach (i; 2 .. nCtrlPts-2) {
            if (dim == 2) {
                myBez.B[i].x = d[2*i - 2];
                myBez.B[i].y = d[2*i - 1];
            }
            else {
                myBez.B[i].x = d[3*i - 3];
                myBez.B[i].y = d[3*i - 2];
                myBez.B[i].z = d[3*i - 1];
            }
        }
        // Then set y and z coordinates of first points in from ends based on slope.
        auto scale0 = (d[$-2] - myBez.B[0].x)/startSlope.x;
        myBez.B[1] = myBez.B[0] + scale0*startSlope;
        auto scale1 = (myBez.B[nCtrlPts-1].x - d[$-1])/endSlope.x;
        myBez.B[nCtrlPts-2] = myBez.B[nCtrlPts-1] - scale1*endSlope;
        // Evaluate error between data points and Bezier curve
        number err = 0.0;
        foreach (p; points) {
            err += myBez.closestDistance(p, t);
        }
        return err;
    }
    // ------------- Done building cost function -----------------------

    // -----------------------------------------------------------------
    // Optimise placement of control points using Nelder-Mead minimiser
    // -----------------------------------------------------------------
    number f_min;
    int n_fe, n_restart;
    number[] dx;
    dx.length = d.length;
    // Make the initial perturbations 1/100th of the arc length.
    number dp = myBez.length()/100.0;
    dx[] = dp;
    success = nm.nelmin.minimize!(funMin,number)(d, f_min, n_fe, n_restart, dx, tol, 1, max_steps);
    foreach (i; 2 .. nCtrlPts-2) {
        if (dim == 2) {
            myBez.B[i].x = d[2*i - 2];
            myBez.B[i].y = d[2*i - 1];
        }
        else {
            myBez.B[i].x = d[3*i - 3];
            myBez.B[i].y = d[3*i - 2];
            myBez.B[i].z = d[3*i - 1];
        }
    }
    // Then set y and z coordinates of first points in from ends based on slope.
    auto scale0 = (d[$-2] - myBez.B[0].x)/startSlope.x;
    myBez.B[1] = myBez.B[0] + scale0*startSlope;
    auto scale1 = (myBez.B[nCtrlPts-1].x - d[$-1])/endSlope.x;
    myBez.B[nCtrlPts-2] = myBez.B[nCtrlPts-1] - scale1*endSlope;
    // Evaluate error between data points and Bezier curve
    // Populate ts vector with all the t-values associated with the data points.
    ts[0] = 0.0;
    foreach (i; 1 .. points.length-1) {
        myBez.closestDistance(points[i], t);
        ts[i] = t;
    }
    ts[$-1] = 1.0;
    return myBez;
    // ------------- Done calling optimiser -----------------------
}
