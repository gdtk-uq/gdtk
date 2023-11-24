// projection.d
// Various operations and tests that are some form of projection or mapping
// in our Vector3 world.

module geom.elements.projection;

import std.math;
import std.conv;
import ntypes.complex;
import nm.number;
import geom.geometry_exception;
import geom.elements.vector3;
import geom.elements.properties;

@nogc
bool intersect2D(const Vector3 p0, const Vector3 p1, const Vector3 ps, const Vector3 d,
                 out double t)
// Determine the intersection of a projected line with a line segment.
// Input:
//     p0->p1 line segment
//     ps starting point for projected line
//     d direction of projected line
// Output:
//     t parametric position of intersection along the line segment
// Returns:
//     true, if the lines cross;
//     false, if the lines are parallel
// See PJ's workbook page 34, 2017-06-24 for notation and derivation.
{
    number delx = p1.x-p0.x;
    number dely = p1.y-p0.y;
    number denom = delx*d.y - dely*d.x;
    number segmentSize = fabs(delx)+fabs(dely);
    if (fabs(denom) < 1.0e-16*segmentSize) { return false; } // d is parallel to line segment
    t = ((ps.x-p0.x)*d.y - (ps.y-p0.y)*d.x).re/denom.re;
    return true;
} // end intersect2D()

//------------------------------------------------------------------------
// Geometry functions projection and mapping.

/**
 * Project the point q onto the plane, along the vector qr.
 */
int project_onto_plane(ref Vector3 q, ref const(Vector3) qr,
                       ref const(Vector3) a, ref const(Vector3) b, ref const(Vector3) c)
{
    // See Section 7.2 in
    // J. O'Rourke (1998)
    // Computational Geometry in C (2nd Ed.)
    // Cambridge Uni Press
    // See 3D CFD workbook p17, 25-Jan-2006

    // Define a plane Ax + By + Cz = D using the corners of the triangle abc.
    Vector3 N = cross(a-c, b-c); // N = Vector3(A, B, C)
    number D = dot(a, N);

    number numer = D - dot(q, N);
    number denom = dot(qr, N);

    double tol = 1.0e-12;  // floating point tolerance
    if ( fabs(denom) < tol ) {
        if ( fabs(numer) < tol ) {
            return 1;  // qr is parallel to the plane and q is on the plane
        } else {
            return 2;  // qr is parallel to the plane and q is off the plane
        }
    } else {
        q = q + (numer/denom) * qr;
        return 0;  // point q has been projected onto the plane.
    }
} // end project_onto_plane()

/// Map space so that a neutral plane wraps onto a cylinder of radius H.
void map_neutral_plane_to_cylinder(ref Vector3 p, number H)
{
    // The axis of the hypothetical cylinder coincides with the x-axis thus
    // H is also the distance of the neutral plane above the x-axis.
    // For Hannes Wojciak and Paul Petrie-Repar's turbomachinery grids.
    if ( H > 0.0 ) {
        number theta = p.y / H;
        p.set(p.x, p.z*sin(theta), p.z*cos(theta));
    }
    return;
}

/**
 * Find Barycentric Coordinates of point p in triangle,
 * considering (x,y)-plane components, only.
 *  p2
 *   |\
 *   | \
 *   |  \
 *   |   \
 *  p0----p1
*/
@nogc
number[3] barycentricCoords(ref const(Vector3) p, ref const(Vector3) p0,
                            ref const(Vector3) p1, ref const(Vector3) p2)
{
    // Transcribe equations from Wikipedia article "Barycentric coordinate system",
    // subsection on "Conversion between barycentric and Cartesian coordinates".
    number[3] lmbda;
    number numer0 = (p1.y-p2.y)*(p.x-p2.x) + (p2.x-p1.x)*(p.y-p2.y);
    number denom = (p1.y-p2.y)*(p0.x-p2.x) + (p2.x-p1.x)*(p0.y-p2.y);
    lmbda[0] = numer0/denom;
    number numer1 = (p2.y-p0.y)*(p.x-p2.x) + (p0.x-p2.x)*(p.y-p2.y);
    lmbda[1] = numer1/denom;
    lmbda[2] = 1.0 - lmbda[0] - lmbda[1];
    return lmbda;
} // end barycentricCoords()

@nogc
bool is_outside_triangle(number[3] lmbda, double tol=1.0e-10)
{
    bool is_outside = false;
    if (lmbda[0] < -tol) { is_outside = true; }
    if (lmbda[1] < -tol) { is_outside = true; }
    if (lmbda[2] < -tol) { is_outside = true; }
    // Now, we have determined that the coordinates are for a point
    // inside or on an edge or at a vertex, to within the tolerance.
    return is_outside;
} // end is_outside_triangle()

/**
 * Find Barycentric Coordinates of point p in a quadrilateral polygon,
 * considering (x,y)-plane components, only.
 *  p3-----p2
 *   |\   /|
 *   | \ / |
 *   |  p  |
 *   | / \ |
 *   |/   \|
 *  p0-----p1
*/
@nogc
number[4] barycentricCoords(ref const(Vector3) p,
                            ref const(Vector3) p0, ref const(Vector3) p1,
                            ref const(Vector3) p2, ref const(Vector3) p3)
{
    // Compute the normalized barycentric coordinates for a point within
    // a convex quadrilateral polygon, as described in Michael Floater's
    // 2003 paper "Mean Value Coordinates" in the journal
    // Computer Aided Geometric Design Vol 20, issue 1, pages 19-27.
    //
    // Signed areas of triangular corners.
    number[4] Acorner;
    Acorner[0] = xyplane_area(p3, p0, p1);
    Acorner[1] = xyplane_area(p0, p1, p2);
    Acorner[2] = xyplane_area(p1, p2, p3);
    Acorner[3] = xyplane_area(p2, p3, p0);
    if ((Acorner[0] < 0.0) || (Acorner[1] < 0.0) ||
        (Acorner[2] < 0.0) || (Acorner[3] < 0.0)) {
        throw new GeometryException("Quadrilateral is not convex.");
    }
    // Signed areas of triangles associated with edges and point p.
    number[4] Amid;
    Amid[0] = xyplane_area(p0, p1, p);
    Amid[1] = xyplane_area(p1, p2, p);
    Amid[2] = xyplane_area(p2, p3, p);
    Amid[3] = xyplane_area(p3, p0, p);
    // Wachpress coordinates computed as per the final remarks in Floater's paper.
    number[4] lmbda;
    lmbda[0] = Acorner[0] * Amid[1]*Amid[2];
    lmbda[1] = Acorner[1] * Amid[2]*Amid[3];
    lmbda[2] = Acorner[2] * Amid[3]*Amid[0];
    lmbda[3] = Acorner[3] * Amid[0]*Amid[1];
    number scale = 1.0/(lmbda[0]+lmbda[1]+lmbda[2]+lmbda[3]);
    foreach (i; 0 .. 4) { lmbda[i] *= scale; }
    return lmbda;
} // end barycentricCoords()

@nogc
bool is_outside_quad(number[4] lmbda, double tol=1.0e-10)
{
    bool is_outside = false;
    foreach (i; 0 .. 4) {
        if (lmbda[i] < -tol) { is_outside = true; }
    }
    // Now, we have determined that the coordinates are for a point
    // inside or on an edge or at a vertex, to within the tolerance.
    return is_outside;
} // end is_outside_quad()


version(projection_test) {
    import util.msg_service;
    int main() {
        double t_intersection;
        bool foundIntersection = intersect2D(Vector3(0.0,1.0), Vector3(1.0,1.0),
                                             Vector3(0.5,0.5), Vector3(0.0,1.0),
                                             t_intersection);
        assert(foundIntersection, failedUnitTest());
        assert(isClose(t_intersection, 0.5), failedUnitTest());

        // Projection onto a plane.
        Vector3 a = Vector3(1.0, 0.0, 0.0); // plane through a,b,c
        Vector3 b = Vector3(1.0, 1.0, 0.0);
        Vector3 c = Vector3(0.5, 0.0, 0.0);
        Vector3 qr = Vector3(3.0, 3.0, -3.0); // direction
        Vector3 q = Vector3(0.0, 0.0, 1.0); // start point
        int flag =  project_onto_plane(q, qr, a, b, c);
        assert(approxEqualVectors(q, Vector3(1.0,1.0,0.0)), failedUnitTest());

        // Projection onto a plane - complex step derivative test
        version(complex_numbers) {
            // store original projected point q
            Vector3 q0 = Vector3(q);
            // reset q vector
            q = Vector3(0.0, 0.0, 1.0); // start point

            // Complex Step
            number hIm = complex(0.0, 1.0e-20); // complex step-size
            qr.z += hIm; // perturb in complex plane
            flag =  project_onto_plane(q, qr, a, b, c);
            double[3] qDerivCmplx;
            qDerivCmplx[0] = q.x.im/hIm.im;
            qDerivCmplx[1] = q.y.im/hIm.im;
            qDerivCmplx[2] = q.z.im/hIm.im;

            // reset q & qr vectors
            qr = Vector3(3.0, 3.0, -3.0); // direction
            q = Vector3(0.0, 0.0, 1.0); // start point

            // Real Step
            double hRe = 1.0e-06; // real step-size
            qr.z += hRe; // perturb in real plane
            flag =  project_onto_plane(q, qr, a, b, c);
            double[3] qDerivReal;
            qDerivReal[0] = (q.x.re-q0.x.re)/hRe;
            qDerivReal[1] = (q.y.re-q0.y.re)/hRe;
            qDerivReal[2] = (q.z.re-q0.z.re)/hRe;

            foreach( idx; 0..3) {
                // debug { import std.stdio; writeln("cmplx=", qDerivCmplx[idx], " real=", qDerivReal[idx]); }
                assert(std.math.isClose(qDerivCmplx[idx], qDerivReal[idx], 1.0e-5, 1.0e-9), failedUnitTest());
            }
        }

        // projection onto a cylinder.
        Vector3 myp = Vector3(1.0, 1.0, 1.0);
        map_neutral_plane_to_cylinder(myp, to!number(1.0));
        assert(approxEqualVectors(myp, Vector3(1.0, sin(1.0), cos(1.0))), failedUnitTest());

        // projection onto a cylinder - complex step derivative test
        version(complex_numbers) {
            // store original projected point myp
            Vector3 myp0 = Vector3(myp);

            // reset myp vector
            myp = Vector3(1.0, 1.0, 1.0); // start point

            // Complex Step
            myp.z += hIm; // perturb in complex plane, reuse hIm
            map_neutral_plane_to_cylinder(myp, to!number(1.0));
            qDerivCmplx[0] = myp.x.im/hIm.im;
            qDerivCmplx[1] = myp.y.im/hIm.im;
            qDerivCmplx[2] = myp.z.im/hIm.im;

            // reset myp vector
            myp = Vector3(1.0, 1.0, 1.0); // start point

            // Real Step
            myp.z += hRe; // perturb in real plane, reuse hRe
            map_neutral_plane_to_cylinder(myp, to!number(1.0));
            qDerivReal[0] = (myp.x.re-myp0.x.re)/hRe;
            qDerivReal[1] = (myp.y.re-myp0.y.re)/hRe;
            qDerivReal[2] = (myp.z.re-myp0.z.re)/hRe;
            foreach( idx; 0..3) assert(std.math.isClose(qDerivCmplx[idx], qDerivReal[idx]), failedUnitTest());
        }
        //
        // Barycentric coordinates for the midpoint of an equilateral triangle.
        Vector3 p0 = Vector3(0.0, 0.0);
        Vector3 p1 = Vector3(1.0, 0.0);
        Vector3 p2 = Vector3(0.5, std.math.sin(60.0*std.math.PI/180.0));
        Vector3 p = (p0+p1+p2)/3.0;
        number[3] bcc = barycentricCoords(p, p0, p1, p2);
        assert(isClose(bcc[0].re, 1.0/3.0), failedUnitTest());
        assert(isClose(bcc[1].re, 1.0/3.0), failedUnitTest());
        assert(isClose(bcc[2].re, 1.0/3.0), failedUnitTest());
        assert(!is_outside_triangle(bcc), failedUnitTest());
        // Now move the point down to the bottom edge of the triangle.
        p = Vector3(0.5, 0.0);
        bcc = barycentricCoords(p, p0, p1, p2);
        assert(!is_outside_triangle(bcc), failedUnitTest());
        // Finally, move the point just outside the triangle.
        p = Vector3(0.5, -0.01);
        bcc = barycentricCoords(p, p0, p1, p2);
        assert(is_outside_triangle(bcc), failedUnitTest());
        //
        // Barycentric coordinates for the midpoint of a square.
        p0 = Vector3(0.0, 0.0);
        p1 = Vector3(1.0, 0.0);
        p2 = Vector3(1.0, 1.0);
        Vector3 p3 = Vector3(0.0, 1.0);
        p = (p0+p1+p2+p3)/4.0;
        number[4] bcc4 = barycentricCoords(p, p0, p1, p2, p3);
        assert(isClose(bcc4[0].re, 0.25), failedUnitTest());
        assert(isClose(bcc4[1].re, 0.25), failedUnitTest());
        assert(isClose(bcc4[2].re, 0.25), failedUnitTest());
        assert(isClose(bcc4[3].re, 0.25), failedUnitTest());
        assert(!is_outside_quad(bcc4), failedUnitTest());
        // Now move the point down to the bottom edge of the square.
        p = Vector3(0.5, 0.0);
        bcc4 = barycentricCoords(p, p0, p1, p2, p3);
        assert(!is_outside_quad(bcc4), failedUnitTest());
        // Finally, move the point just outside the square.
        p = Vector3(0.5, -0.01);
        bcc4 = barycentricCoords(p, p0, p1, p2, p3);
        assert(is_outside_quad(bcc4), failedUnitTest());
        //
        return 0;
    }
} // end projection_test
