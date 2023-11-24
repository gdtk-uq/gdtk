/** path.d
 * Geometry-building elements for our 3D world -- one-parameter elements.
 * Note that these are geometric paths, as distinct from the file-system paths.
 *
 * Author: Peter J and Rowan G.
 * Version: 2015-02-19 first code
 *          2015-04-21 added Arc, Bezier, Polyline
 *          2015-07-02 simplify classes
 *          2017-11-26 refactor package
 */

module geom.gpath.path;

import std.conv;
import std.math;
import std.stdio;
import std.string;
import geom.elements;
import geom.gpath.line;
import geom.gpath.arc;
import geom.gpath.bezier;
import geom.gpath.nurbs;
import geom.gpath.polyline;
import geom.surface;
import ntypes.complex;
import nm.number;
import nm.bbla;
import nm.linesearch;

class Path {
public:
    abstract Path dup() const;
    abstract Vector3 opCall(double t) const;
    Vector3 dpdt(double t) const
    {
        // Obtain the derivative approximately, via a finite-difference.
        double dt = 0.001;
        Vector3 p0 = this.opCall(t);
        Vector3 derivative;
        if ( t+dt > 1.0 ) {
            // t is close to the t=1.0 boundary, use a one-sided difference.
            Vector3 pminus1 = this.opCall(t-dt);
            derivative = (p0 - pminus1) / dt;
        } else if ( t-dt < 0.0 ) {
            // s is close to the s=0 boundary, use a one-sided difference.
            Vector3 pplus1 = this.opCall(t+dt);
            derivative = (pplus1 - p0) / dt;
        } else {
            // Not near a boundary, use central-difference.
            Vector3 pminus1 = this.opCall(t-dt);
            Vector3 pplus1 = this.opCall(t+dt);
            derivative = (pplus1 - pminus1) / (2.0 * dt);
        }
        return derivative;
    }
    Vector3 d2pdt2(double t) const
    {
        // Obtain the derivative approximately, via a finite-difference.
        double dt = 0.001;
        Vector3 p0 = this.opCall(t);
        Vector3 derivative;
        if ( t+dt > 1.0 ) {
            // t is close to the t=1.0 boundary, use a one-sided difference.
            Vector3 pminus1 = this.opCall(t-dt);
            Vector3 pminus2 = this.opCall(t-2*dt);
            derivative = (p0 - 2*pminus1 + pminus2) / (dt*dt);
        } else if ( t-dt < 0.0 ) {
            // s is close to the s=0 boundary, use a one-sided difference.
            Vector3 pplus1 = this.opCall(t+dt);
            Vector3 pplus2 = this.opCall(t+2*dt);
            derivative = (pplus2 - 2*pplus1 + p0) / (dt*dt);
        } else {
            // Not near a boundary, use central-difference.
            Vector3 pminus1 = this.opCall(t-dt);
            Vector3 pplus1 = this.opCall(t+dt);
            derivative = (pplus1 - 2*p0 + pminus1) / (dt*dt);
        }
        return derivative;
    }
    number partial_length(double ta, double tb) const
    {
        if( tb < ta ) {
            double tmp = ta; ta = tb; tb = tmp;
        }
        number L = 0.0;
        int n = 100;
        double dt = (tb - ta) / n;
        Vector3 p0 = this.opCall(ta);
        Vector3 p1, dp;
        foreach (i; 1 .. n+1) {
            p1 = this.opCall(ta + dt * i);
            dp = p1 - p0;
            L += geom.abs(dp);
            p0 = p1;
        }
        return L;
    }
    number length() const
    {
        return partial_length(0.0, 1.0);
    }
    Vector3 point_from_length(number length, out double t) const
    {
        number L = 0.0;
        int n = 1000;
        double dt = 1.0 / n;
        Vector3 p0 = this.opCall(0.0);
        Vector3 p1, dp;
        foreach (i; 1 .. n+1) {
            p1 = this.opCall(dt * i);
            dp = p1 - p0;
            L += geom.abs(dp);
            p0 = p1;
            if(L > length) {
                t = dt * i;
                return p1;
            }
        }
        t = dt * n;
        return p1;
    }
    abstract override string toString() const;
    abstract string classString() const;
    bool intersect2D(const Vector3 ps, const Vector3 d, out double t, int nseg=20) const
    // Determine the intersection of a projected line on the Path.
    // Input:
    //     ps starting point for projected line
    //     d direction of projected line
    // Output:
    //     t parametric position of intersection along the Path
    // Returns:
    //     true, if the intersection point was located;
    //     false, if the intersection was not found
    // See PJ's workbook page 34, 2017-06-24 for notation and derivation.
    {
        if (cast(Line)this !is null) { nseg = 1; } // straight Line
        double delt = 1.0/nseg;
        double t0 = 0.0; Vector3 p0 = this.opCall(0.0);
        foreach (i; 0 .. nseg) {
            double t1 = delt*(i+1); Vector3 p1 = this.opCall(t1);
            double tOnSegment;
            bool intersectionOK = geom.intersect2D(p0, p1, ps, d, tOnSegment);
            if (intersectionOK && tOnSegment >= 0.0 && tOnSegment <= 1.0) {
                t = t0 + tOnSegment*delt;
                return true;
            }
            t0 = t1; p0 = p1; // for next segment
        }
        return false;
    } // end intersect2D()
    double closestDistance(const Vector3 p, ref double t) const
    {
        // This is the function we'll try to minimise with a line-search method.
        double distanceToPoint(double t)
        {
            auto ptOnPath = opCall(t);
            return distance_between(p, ptOnPath);
        }
        double tL = 0.0;
        double tR = 1.0;
        double tol = 1.0e-6;
        minimize!(distanceToPoint,double)(tL, tR, tol);

        Vector3 ptOnPath;
        // Handle the special cases of being right at the end points.
        if (tL == 0.0) {
            t = 0.0;
        } else if (tR == 1.0) {
            t = 1.0;
        } else {
            t = 0.5*(tL+tR);
        }
        ptOnPath = opCall(t);
        return distance_between(p, ptOnPath);
    } // end closestDistance

} // end class Path
