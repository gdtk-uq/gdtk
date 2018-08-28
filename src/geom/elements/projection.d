// projection.d
// Various operations and tests that are some form of projection or mapping
// in our Vector3 world.

module geom.elements.projection;

import std.math;
import std.conv;
import nm.complex;
import nm.number;
import geom.elements.vector3;

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
ref Vector3 map_neutral_plane_to_cylinder(ref Vector3 p, number H)
{
    // The axis of the hypothetical cylinder coincides with the x-axis thus
    // H is also the distance of the neutral plane above the x-axis.
    // For Hannes Wojciak and Paul Petrie-Repar's turbomachinery grids.
    if ( H > 0.0 ) {
        number theta = p.y / H;
        p.set(p.x, p.z*sin(theta), p.z*cos(theta));
    }
    return p;
}

/**
 * Find Barycentric Coordinates of point P in triangle 
 *  p2
 *   |\
 *   | \
 *   |  \
 *   |   \
 *  p0----p1
*/
@nogc
void P_barycentricCoords(ref const(Vector3) p, ref const(Vector3) p0, 
                         ref const(Vector3) p1, ref const(Vector3) p2, 
                         ref Vector3 Coords, 
                         double tol=1.0e-12, double area_tol=1.0e-20)
{
    number numer0 = (p1.y-p2.y)*(p.x -p2.x) + (p2.x-p1.x)*(p.y -p2.y);
    number denom = (p1.y-p2.y)*(p0.x-p2.x) + (p2.x-p1.x)*(p0.y-p2.y);
    number lambda0 = numer0 / denom;
    if (abs(lambda0) < tol) { lambda0 = 0; }
    number numer1 = (p2.y-p0.y)*(p.x -p2.x) + (p0.x-p2.x)*(p.y -p2.y);
    number lambda1 = numer1 / denom;
    if (abs(lambda1) < tol) { lambda1 = 0; }
    number lambda2 = 1 - lambda0 - lambda1;
    if (abs(lambda2) < tol) { lambda2 = 0; }
    // set Barycentric coordinates.
    Coords.set(lambda0, lambda1, lambda2);

} // end P_barycentricCoords()




version(projection_test) {
    import util.msg_service;
    int main() {
        double t_intersection;
        bool foundIntersection = intersect2D(Vector3(0.0,1.0), Vector3(1.0,1.0),
                                             Vector3(0.5,0.5), Vector3(0.0,1.0),
                                             t_intersection);
        assert(foundIntersection, failedUnitTest());
        assert(approxEqual(t_intersection, 0.5), failedUnitTest());

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
            qr.refz += hIm; // perturb in complex plane
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
            qr.refz += hRe; // perturb in real plane
            flag =  project_onto_plane(q, qr, a, b, c);
            double[3] qDerivReal;
            qDerivReal[0] = (q.x.re-q0.x.re)/hRe;
            qDerivReal[1] = (q.y.re-q0.y.re)/hRe;
            qDerivReal[2] = (q.z.re-q0.z.re)/hRe;
            
            foreach( idx; 0..3) assert(std.math.approxEqual(qDerivCmplx[idx], qDerivReal[idx]), failedUnitTest());
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
            myp.refz += hIm; // perturb in complex plane, reuse hIm
            map_neutral_plane_to_cylinder(myp, to!number(1.0));
            qDerivCmplx[0] = myp.x.im/hIm.im;
            qDerivCmplx[1] = myp.y.im/hIm.im;
            qDerivCmplx[2] = myp.z.im/hIm.im;

            // reset myp vector
            myp = Vector3(1.0, 1.0, 1.0); // start point
                        
            // Real Step
            myp.refz += hRe; // perturb in real plane, reuse hRe
            map_neutral_plane_to_cylinder(myp, to!number(1.0));
            qDerivReal[0] = (myp.x.re-myp0.x.re)/hRe;
            qDerivReal[1] = (myp.y.re-myp0.y.re)/hRe;
            qDerivReal[2] = (myp.z.re-myp0.z.re)/hRe;
            foreach( idx; 0..3) assert(std.math.approxEqual(qDerivCmplx[idx], qDerivReal[idx]), failedUnitTest());
        }

        return 0;
    }
} // end projection_test

