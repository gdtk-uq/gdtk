/*
 * Authors: Rowan G. and Peter J.
 * Date: 2019-04-02
 *
 * The Bezier triangle patch is not part of the ParametricSurface
 * family. Although the Bezier triangle is bivariate (in u and v),
 * the u,v ranges do not *independently* extend from 0 <= u,v <= 1.
 * Rather the u,v parameters are part of a triple (u, v, w) that
 * form barycentric coordinates. These have the constraint:
 *    u + v + w = 1.
 * We can treat the surface as bivariate because given u and v,
 * then w is implied:
 *    w = 1 - u - v
 *
 * We provide this Bezier triangle patch because it is useful as
 * a construction patch and for parameterisation of surfaces for
 * design optimization.
 *
 * The implementation here follows closesly the description in
 * Farin (2001) and Hansford (2002).
 *
 * References
 * ----------
 * Farin (2001)
 * Curves and Surfaces for CAGD: A Practical Guide.
 * Elsevier Science & Technology
 *
 * Hansford (2002)
 * Bezier Techniques, Ch.4 in Handbook of Computer Aided Geometric Design.
 * Editors: Farin et al.
 * Elsevier Science & Technology
 *
 */

module geom.surface.beziertrianglepatch;

import std.stdio;
import std.stdio : File;
import std.conv : text, to;
import std.range : iota;
import std.string : format;
import std.algorithm.comparison : max;
import ntypes.complex;
import nm.number;

import nm.nelmin;

import geom;

@nogc
size_t toSingleIndex()(size_t i, size_t j, size_t k)
{
    auto n = i + j + k;
    auto m = n - i;
    return (m*(m+5)/2) - 2*j - k;
}

class BezierTrianglePatch {
public:
    Vector3[] B; // control points
    int n; // degree of patch

    this(const Vector3[] _B, int _n)
    {
        n = _n;
        auto nCtrlPts = (n+1)*(n+2)/2;
        if (nCtrlPts != _B.length) {
            throw new Error(text("BezierTrianglePatch: incorrect number of control points supplied.",
                                 nCtrlPts, " expected, ", _B.length, " supplied."));
        }
        B.length = nCtrlPts;
        B[] = _B[];
        initialiseWorkingSpace();
    }

    void initialiseWorkingSpace()
    {
        _b.length = n + 1;
        foreach (i; 0 .. n + 1) {
            auto n_local = n - i;
            _b[i].length = (n_local + 1)*(n_local + 2) / 2;
        }
        // Copy B into _b[0].
        // Yes, we are duplicating the storage of those control points,
        // but it makes the implementation of de Casteljau's algorithm
        // much nicer. Typically, the number of control points will
        // be relatively small, so this shouldn't be a big storage penalty.
        _b[0][] = B[];
    }

    void updateWorkingSpace()
    {
        // After changing control points, we need to reflect this
        // in the working space.
        _b[0][] = B[];
    }

    Vector3 opCall(double u, double v)
    {
        auto w = 1.0 - u - v;
        // We may allow u or v to get slightly larger than 1.0
        // if we are perturbing to get sensitivity.
        // In which case, we won't allow w to go negative.
        if (w < 0.0) w = 0.0;

        // de Casteljau algorithm as described in:
        // Farin (2001), pp. 310--311
        foreach (r; 1 .. n + 1) {
            // Loop over points via triangle indices.
            auto n_local = n - r;
            foreach (i; iota(n_local, -1, -1)) {
                foreach (j; iota(n_local-i, -1, -1)) {
                    auto k = n_local - (i+j);
                    auto idx = toSingleIndex!()(i, j, k);
                    auto uIdx = toSingleIndex!()(i+1, j, k);
                    auto vIdx = toSingleIndex!()(i, j+1, k);
                    auto wIdx = toSingleIndex!()(i, j, k+1);
                    // Do Vector3 operation component-wise to
                    // avoid allocation of temporaries.
                    _b[r][idx].x = u*_b[r-1][uIdx].x + v*_b[r-1][vIdx].x + w*_b[r-1][wIdx].x;
                    _b[r][idx].y = u*_b[r-1][uIdx].y + v*_b[r-1][vIdx].y + w*_b[r-1][wIdx].y;
                    _b[r][idx].z = u*_b[r-1][uIdx].z + v*_b[r-1][vIdx].z + w*_b[r-1][wIdx].z;
                }
            }
        }
        return _b[n][0];
    }

    number projectPoint(Vector3 p, ref Vector3 q, ref double uFound, ref double vFound, bool returnOnFailure=false)
    {
        // Establish an initial guess for u and v
        double[] samples = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
                            0.6, 0.7, 0.8, 0.9, 1.0];
        q = opCall(0.0, 0.0);
        number minDist = distance_between(p, q);
        number maxDist = minDist;
        double min_u = 0.0;
        double min_v = 0.0;
        foreach (u; samples) {
            foreach (v; samples) {
                if (u + v > 1.0) break;
                q = opCall(u, v);
                number dist = distance_between(p, q);
                if (dist < minDist) {
                    minDist = dist;
                    min_u = u;
                    min_v = v;
                }
                maxDist = max(maxDist, dist);
            }
        }

        // Define function to minimize
        number penalty = 1000 * maxDist;
        number fMin(number[] x)
        {
            double u = x[0].re;
            double v = x[1].re;
            // Check that u and v are feasible.
            // Apply penalty if not.
            // We let the values of u and v wander off the
            // edge of the patch a little. This produces
            // a better behaved approach to the answer
            // than applying a harsh penalty right at the edge.
            if (u < -0.1) return penalty;
            if (v < -0.1) return penalty;
            if (u + v > 1.1) return penalty;
            // Now proceed to evaluate distance.
            q = opCall(u, v);
            return to!number(distance_between(p, q));
        }

        // Set up initial guess
        number[] x = [to!number(min_u), to!number(min_v)];

        // Call optimizer
        number f_min;
        int n_fe, n_restart;
        number[] dx = [to!number(0.01), to!number(0.01)];
        double tol = 1.0e-8;
        int max_steps = 10000;

        bool success = nm.nelmin.minimize!(fMin, number)(x, f_min, n_fe, n_restart, dx, tol, 1, max_steps);

        if (success || returnOnFailure) {
            uFound = x[0].re;
            vFound = x[1].re;
            q = opCall(uFound, vFound);
            return to!number(distance_between(p, q));
        }

        // Otherwise, we failed in the minimize step.
        string errMsg = "Error in BezierTrianglePatch.projectPoint()\n";
        errMsg ~= "Optimizer stats:\n";
        errMsg ~= format("  number of function evaluations: %d\n", n_fe);
        errMsg ~= format("  number of restarts: %d\n", n_restart);
        errMsg ~= format("  minimum of function found: %12.6e\n", f_min);
        throw new Exception(errMsg);
    }

private:
    Vector3[][] _b; // working space for de Casteljau algorithm

} // end class BezierTrianglePatch

void writeBezierTriangleCtrlPtsAsText(BezierTrianglePatch btp, string fileName)
{
    auto f = File(fileName, "w");
    foreach (ref p; btp.B) {
        f.writefln("%20.16e %20.16e %20.16e", p.x.re, p.y.re, p.z.re);
    }
    f.close();
}

void writeBezierTriangleCtrlPtsAsVtkXml(BezierTrianglePatch btp, string fileName)
{
    auto nPtsTotal = btp.B.length;
    auto n = btp.n;
    auto nCells = n^^2;

    auto f = File(fileName, "w");
    f.writeln("<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" header_type=\"UInt64\">");
    f.writeln("  <UnstructuredGrid>");
    f.writefln("    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">",  nPtsTotal, nCells);
    f.writeln("       <Points>");
    f.writeln("         <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">");
    foreach (ref p; btp.B) {
        f.writefln("       %20.16e %20.16e %20.16e", p.x.re, p.y.re, p.z.re);
    }
    f.writeln("        </DataArray>");
    f.writeln("      </Points>");
    f.writeln("      <Cells>");
    f.writeln("        <DataArray type=\"Int32\" Name=\"connectivity\">");
    foreach (i; iota(n, 0, -1)) {
            foreach (j; iota(n-i, -1, -1)) {
                auto k = n - (i+j);
                // Triangle pattern A
                auto idx0 = toSingleIndex!()(i,j,k);
                auto idx1 = toSingleIndex!()(i-1,j+1,k);
                auto idx2 = toSingleIndex!()(i-1,j,k+1);
                f.writefln("       %d %d %d", idx0, idx1, idx2);
                // Triangle pattern B
                if (j == 0) continue;
                idx1 = toSingleIndex!()(i-1,j,k+1);
                idx2 = toSingleIndex!()(i,j-1,k+1);
                f.writefln("       %d %d %d", idx0, idx1, idx2);
            }
    }
    f.writeln("         </DataArray>");
    f.writeln("         <DataArray type=\"Int32\" Name=\"offsets\">");
    int offset = 0;
    foreach (i; 0 .. nCells) {
        offset += 3;
        f.writefln("        %d", offset);
    }
    f.writeln("         </DataArray>");
    f.writeln("         <DataArray type=\"UInt8\" Name=\"types\">");
    foreach (i; 0 .. nCells) {
        f.writeln("          5");
    }
    f.writeln("         </DataArray>");
    f.writeln("       </Cells>");
    f.writeln("    </Piece>");
    f.writeln("  </UnstructuredGrid>");
    f.writeln("</VTKFile>");
    f.close();
}

void writeBezierTriangleAsVtkXml(BezierTrianglePatch btp, string fileName, int nEdgePts)
{
    int nPtsTotal = nEdgePts*(nEdgePts+1)/2;
    int n = nEdgePts - 1;
    int nCells = (nEdgePts-1)^^2;

    auto f = File(fileName, "w");
    f.writeln("<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" header_type=\"UInt64\">");
    f.writeln("  <UnstructuredGrid>");
    f.writefln("    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">",  nPtsTotal, nCells);
    f.writeln("       <Points>");
    f.writeln("         <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">");
    foreach (i; iota(n, -1, -1)) {
        foreach (j; iota(n-i, -1, -1)) {
            double u = double(i)/n;
            double v = double(j)/n;
            auto p = btp(u, v);
            f.writefln("       %20.16e %20.16e %20.16e", p.x, p.y, p.z);
        }
    }
    f.writeln("        </DataArray>");
    f.writeln("      </Points>");
    f.writeln("      <Cells>");
    f.writeln("        <DataArray type=\"Int32\" Name=\"connectivity\">");
    foreach (i; iota(n, 0, -1)) {
            foreach (j; iota(n-i, -1, -1)) {
                auto k = n - (i+j);
                // Triangle pattern A
                auto idx0 = toSingleIndex!()(i,j,k);
                auto idx1 = toSingleIndex!()(i-1,j+1,k);
                auto idx2 = toSingleIndex!()(i-1,j,k+1);
                f.writefln("       %d %d %d", idx0, idx1, idx2);
                // Triangle pattern B
                if (j == 0) continue;
                idx1 = toSingleIndex!()(i-1,j,k+1);
                idx2 = toSingleIndex!()(i,j-1,k+1);
                f.writefln("       %d %d %d", idx0, idx1, idx2);
            }
    }
    f.writeln("         </DataArray>");
    f.writeln("         <DataArray type=\"Int32\" Name=\"offsets\">");
    int offset = 0;
    foreach (i; 0 .. nCells) {
        offset += 3;
        f.writefln("        %d", offset);
    }
    f.writeln("         </DataArray>");
    f.writeln("         <DataArray type=\"UInt8\" Name=\"types\">");
    foreach (i; 0 .. nCells) {
        f.writeln("          5");
    }
    f.writeln("         </DataArray>");
    f.writeln("       </Cells>");
    f.writeln("    </Piece>");
    f.writeln("  </UnstructuredGrid>");
    f.writeln("</VTKFile>");
    f.close();
}

void writeBezierTriangleAsDat(BezierTrianglePatch btp, string fileName, int nEdgePts)
{
    int nPtsTotal = nEdgePts*(nEdgePts+1)/2;
    int n = nEdgePts - 1;
    int nCells = (nEdgePts-1)^^2;

    auto f = File(fileName, "w");
    f.writefln("%d",  nPtsTotal);
    foreach (i; iota(n, -1, -1)) {
        foreach (j; iota(n-i, -1, -1)) {
            double u = double(i)/n;
            double v = double(j)/n;
            auto p = btp(u, v);
            f.writefln("       %20.16e %20.16e %20.16e", p.x, p.y, p.z);
        }
    }
    f.close();
}


BezierTrianglePatch bezierTriangleFromPointCloud(Vector3[] points, Vector3 p0, Vector3 p1, Vector3 p2, int n, BezierTrianglePatch initGuess)
{
    // NOTE: This function has no constraint at the edges. It's quite possible to find a patch that fits the point
    //       cloud well, but the patch overlaps at the edges. It's not usually a big deal, just something to be
    //       aware of.
    double tol = 1.0e-3;
    int maxSteps = 10000;
    // Corners are located at
    //
    //            p1
    //             + b(0,n,0)
    //            / \
    //           /   \
    //       p2 +_____+ p0
    //   b(0,0,n)        b(n,0,0)
    //

    auto nCtrlPts = (n+1)*(n+2)/2;
    Vector3[] ctrlPts;
    ctrlPts.length = nCtrlPts;
    number[] d;

    // -----------------------------------------------------------------
    // Establish an initial guess for the location of control points
    // -----------------------------------------------------------------
    // We know the corner points. Put them in the list.
    // These points are not to be touched by the optimiser.
    ctrlPts[toSingleIndex!()(n,0,0)] = p0;
    ctrlPts[toSingleIndex!()(0,n,0)] = p1;
    ctrlPts[toSingleIndex!()(0,0,n)] = p2;
    // The user might have supplied an initial guess.
    // We can use this if the degree is correct.
    if (initGuess !is null && initGuess.n == n) {
        foreach (i; iota(n, -1, -1)) {
            foreach (j; iota(n-i, -1, -1)) {
                auto k = n - (i+j);
                // Skip corner points
                if (i == n || j == n || k == n) continue;
                auto idx = toSingleIndex!()(i,j,k);
                ctrlPts[idx] = initGuess.B[idx];
                d ~= ctrlPts[idx].x;
                d ~= ctrlPts[idx].y;
                d ~= ctrlPts[idx].z;
            }
        }
    }
    else {
        // Distribute the remaining control points roughly equally using barycentric coordinates
        double du = 1.0/n;
        double dv = du;
        foreach (i; iota(n, -1, -1)) {
            foreach (j; iota(n-i, -1, -1)) {
                auto k = n - (i+j);
                // Skip corner points
                if (i == n || j == n || k == n) continue;
                double u = i*du;
                double v = j*dv;
                double w = 1.0 - u - v;
                auto idx = toSingleIndex!()(i,j,k);
                ctrlPts[idx] = u*p0 + v*p1 + w*p2;
                d ~= ctrlPts[idx].x;
                d ~= ctrlPts[idx].y;
                d ~= ctrlPts[idx].z;
            }
        }
    }
    auto myBezTriPatch = new BezierTrianglePatch(ctrlPts, n);
    // -------------- Done: establishing starting guess ------------------------

    // --------------------------------------------------------------
    // Build cost function to be minimized.
    // --------------------------------------------------------------
    number fMin(number[] d)
    {
        // Adjust control points based on supplied design values 'd'.
        size_t pos = 0;
        foreach (i; iota(n, -1, -1)) {
            foreach (j; iota(n-i, -1, -1)) {
                auto k = n - (i+j);
                // Skip corner points
                if (i == n || j == n || k == n) continue;
                auto idx = toSingleIndex!()(i,j,k);
                myBezTriPatch.B[idx].x = d[pos]; pos++;
                myBezTriPatch.B[idx].y = d[pos]; pos++;
                myBezTriPatch.B[idx].z = d[pos]; pos++;
            }
        }
        // Remember to update working space after changing control
        // point locations.
        myBezTriPatch.updateWorkingSpace();
        // Evaluate error between point cloud and triangle batch
        number err = 0.0;
        Vector3 q;
        double uFound, vFound;
        foreach (p; points) {
            err += myBezTriPatch.projectPoint(p, q, uFound, vFound, true);
        }
        return err;
    }
    // ----------- Done: building cost function ----------------------------

    // ---------------------------------------------------------------------
    // Optimize the placement of control points using Nelder-Mead minimiser
    // ---------------------------------------------------------------------
    number f_min;
    int n_fe, n_restart;
    number[] dx;
    dx.length = d.length;
    // Make the initial perturbations of control points 1/100th of the longest side length
    number dp = 0.01*max(distance_between(p0, p1), distance_between(p0, p2), distance_between(p1, p2));
    dx[] = dp;
    bool success = nm.nelmin.minimize!(fMin, number)(d, f_min, n_fe, n_restart, dx, tol, 1, maxSteps);
    if (success) {
        size_t pos = 0;
        foreach (i; iota(n, -1, -1)) {
            foreach (j; iota(n-i, -1, -1)) {
                auto k = n - (i+j);
                // Skip corner points
                if (i == n || j == n || k == n) continue;
                auto idx = toSingleIndex!()(i,j,k);
                myBezTriPatch.B[idx].x = d[pos]; pos++;
                myBezTriPatch.B[idx].y = d[pos]; pos++;
                myBezTriPatch.B[idx].z = d[pos]; pos++;
            }
        }
        // Remember to update working space after changing control
        // point locations.
        myBezTriPatch.updateWorkingSpace();
        return myBezTriPatch;
    }
    // Otherwise, we failed in the optimizer.
    string errMsg = "Error in bezierTriangleFromPointCloud().\n";
    errMsg ~= "Optimizer stats:\n";
    errMsg ~= format("  number of function evaluations: %d\n", n_fe);
    errMsg ~= format("  number of restarts: %d\n", n_restart);
    errMsg ~= format("  minimum of function found: %12.6e\n", f_min);
    throw new Exception(errMsg);
}

BezierTrianglePatch bezierTriangleFromPointCloud(Vector3[] points, Bezier b0, Bezier b1, Bezier b2, int n, BezierTrianglePatch initGuess, ref bool success)
{
    double tol = 0.01;
    int maxSteps = 500;
    // Edges and their directions are defined as:
    // (_,_,_) = (u,v,w)
    // b_ = Bezier curve
    //
    //                     + (0,1,0)
    //              _0_   / \   *1
    //               /   /   \   \
    //          b0--/   /     \   \--b2
    //             /   /       \   \
    //            /   /         \   \
    //          1*   /           \  -0-
    //              +_____________+
    //           (0,0,1)         (1,0,0)
    //
    //              1*-----------|0
    //                     |
    //                     b1
    //
    // b0 : defines Bezier triangle curve u=0
    //    : runs from v=1 --> w=1
    //
    // b1 : defines Bezier triangle curve v=0
    //    : runs from u=1 --> w=1
    //
    // b2 : defines Bezier triangle curve w=0
    //    : runs from u=1 --> v=1
    //

    // Perform some checks that Bezier end points are (approximately) coincident
    if (!approxEqualVectors(b0.B[0], b2.B[$-1])) {
        throw new Error("Error: b0 and b2 end points are not coincident.\n");
    }
    if (!approxEqualVectors(b0.B[$-1], b1.B[$-1])) {
        throw new Error("Error: b0 and b1 end points are not coincident.\n");
    }
    if (!approxEqualVectors(b1.B[0], b2.B[0])) {
        throw new Error("Error: b1 and b2 end points are not coincident.\n");
    }

    auto nCtrlPts = (n+1)*(n+2)/2;
    Vector3[] ctrlPts;
    ctrlPts.length = nCtrlPts;
    number[] d;

    // -----------------------------------------------------------------
    // Establish an initial guess for the location of control points
    // -----------------------------------------------------------------
    // We know the edge points. Put them in the list.
    // These points are not to be touched by the optimiser.
    // Work down b0 curve
    foreach (i, b; b0.B) {
        ctrlPts[$-n-1+i] = b;
    }
    // Work down b1 curve
    foreach (i, b; b1.B) {
        auto idx = (i+1)*(i+2)/2 - 1;
        ctrlPts[idx] = b;
    }
    // Work down b2 curve
    foreach (i, b; b2.B) {
        auto idx = i*(i+1)/2;
        ctrlPts[idx] = b;
    }
    // The user might have supplied an initial guess.
    // We can use this if the degree is correct.
    if (initGuess !is null && initGuess.n == n) {
        foreach (i; iota(n, -1, -1)) {
            foreach (j; iota(n-i, -1, -1)) {
                auto k = n - (i+j);
                // Skip edge points
                if (i == 0 || j == 0 || k == 0) continue;
                auto idx = toSingleIndex!()(i,j,k);
                ctrlPts[idx] = initGuess.B[idx];
                d ~= ctrlPts[idx].x;
                d ~= ctrlPts[idx].y;
                d ~= ctrlPts[idx].z;
            }
        }
    }
    else {
        // Distribute the remaining control points roughly equally using barycentric coordinates
        double du = 1.0/n;
        double dv = du;
        foreach (i; iota(n, -1, -1)) {
            foreach (j; iota(n-i, -1, -1)) {
                auto k = n - (i+j);
                // Skip edge points
                if (i == 0 || j == 0 || k == 0) continue;
                double u = i*du;
                double v = j*dv;
                double w = 1.0 - u - v;
                auto idx = toSingleIndex!()(i,j,k);
                ctrlPts[idx] = u*ctrlPts[0] + v*ctrlPts[$-n-1] + w*ctrlPts[$-1];
                d ~= ctrlPts[idx].x;
                d ~= ctrlPts[idx].y;
                d ~= ctrlPts[idx].z;
            }
        }
    }
    auto myBezTriPatch = new BezierTrianglePatch(ctrlPts, n);
    // -------------- Done: establishing starting guess ------------------------

    // --------------------------------------------------------------
    // Build cost function to be minimized.
    // --------------------------------------------------------------
    int count = 0;
    number fMin(number[] d)
    {
        // Adjust control points based on supplied design values 'd'.
        size_t pos = 0;
        foreach (i; iota(n, -1, -1)) {
            foreach (j; iota(n-i, -1, -1)) {
                auto k = n - (i+j);
                // Skip edge points
                if (i == 0 || j == 0 || k == 0) continue;
                auto idx = toSingleIndex!()(i,j,k);
                //writeln("i= ", i, " j= ", j, " k= ", k);
                //writeln("idx= ", idx);
                myBezTriPatch.B[idx].x = d[pos]; pos++;
                myBezTriPatch.B[idx].y = d[pos]; pos++;
                myBezTriPatch.B[idx].z = d[pos]; pos++;
            }
        }
        // Remember to update working space after changing control
        // point locations.
        myBezTriPatch.updateWorkingSpace();
        // Evaluate error between point cloud and triangle batch
        number err = 0.0;
        Vector3 q;
        double uFound, vFound;
        foreach (p; points) {
            err += myBezTriPatch.projectPoint(p, q, uFound, vFound, true);
        }
        count++;
        //import std.stdio;
        //writefln("count= %d, fMin= %e", count, err);
        //import core.stdc.stdlib : exit;
        //exit(1);
        return err;
    }
    // ----------- Done: building cost function ----------------------------

    // ---------------------------------------------------------------------
    // Optimize the placement of control points using Nelder-Mead minimiser
    // ---------------------------------------------------------------------
    number f_min;
    int n_fe, n_restart;
    number[] dx;
    dx.length = d.length;
    // Make the initial perturbations of control points 1/100th of the longest side length
    number dp = 0.01*max(distance_between(b0.B[0], b0.B[$-1]),
                         distance_between(b1.B[0], b1.B[$-1]),
                         distance_between(b2.B[0], b2.B[$-1]));
    dx[] = dp;
    success = nm.nelmin.minimize!(fMin, number)(d, f_min, n_fe, n_restart, dx, tol, 1, maxSteps);

    size_t pos = 0;
    foreach (i; iota(n, -1, -1)) {
        foreach (j; iota(n-i, -1, -1)) {
            auto k = n - (i+j);
            // Skip edge points
            if (i == 0 || j == 0 || k == 0) continue;
            auto idx = toSingleIndex!()(i,j,k);
            myBezTriPatch.B[idx].x = d[pos]; pos++;
            myBezTriPatch.B[idx].y = d[pos]; pos++;
            myBezTriPatch.B[idx].z = d[pos]; pos++;
        }
    }
    // Remember to update working space after changing control
    // point locations.
    myBezTriPatch.updateWorkingSpace();
    //writeln("Final ctrl pts:");
    //writeln(myBezTriPatch.B);
    return myBezTriPatch;
}

version(beziertrianglepatch_test) {
    import util.msg_service;
    import std.stdio;
    import std.math;
    int main()
    {
        // Example on 17.1 on p. 312 in Farin (2001)
        Vector3[] B = [Vector3(6, 0, 9),
                       Vector3(3, 3, 6), Vector3(3, 0, 0),
                       Vector3(0, 6, 0), Vector3(0, 3, 0), Vector3(0, 0, 0)];
        auto btp = new BezierTrianglePatch(B, 2);
        double u = 1./3;
        double v = 1./3;
        auto p = btp(u, v);
        // Farin gives results as (2, 2, 7/3);
        assert(approxEqualVectors(p, Vector3(2, 2, 7./3)), failedUnitTest());

        // Test projection of point that is ON surface
        Vector3 q;
        auto dist = btp.projectPoint(p, q, u, v);
        assert(dist < 1.0e-6, failedUnitTest());
        assert(approxEqualVectors(q, Vector3(2, 2, 7./3)), failedUnitTest());
        assert(isClose(u, 1./3, 1.0e-6), failedUnitTest());
        assert(isClose(v, 1./3, 1.0e-6), failedUnitTest());

        // Test optimizer with this surface.
        // 1. Generate a bunch of points
        Vector3[] points;
        double du = 1.0/4;
        double dv = du;
        foreach (uu; iota(0.0, 1.0+0.5*du, du)) {
            foreach (vv; iota(0.0, (1.0-uu)+0.5*dv, dv)) {
                if ((uu + vv) > 1.0) continue;
                points ~= btp(uu, vv);
            }
        }
        // 2. Find Bezier triangle to fit this point cloud.
        auto p0 = Vector3(6, 0, 9);
        auto p1 = Vector3(0, 6, 0);
        auto p2 = Vector3(0, 0, 0);
        // Let's perturb the points a little as an initial guess.
        B[1].x = 3.1; B[1].y = 2.9; B[1].z = 6.2;
        B[2].x = 2.8; B[2].y = 0.15; B[2].z = -0.05;
        B[4].x = -0.3; B[4].y = 2.85; B[4].z = -0.1;
        auto guess = new BezierTrianglePatch(B, 2);
        auto testPatch = bezierTriangleFromPointCloud(points, p0, p1, p2, 2, guess);

        p = testPatch(0.2, 0.3);
        q = btp(0.2, 0.3);
        assert(distance_between(p, q) < 0.1);

        return 0;
    }

}
