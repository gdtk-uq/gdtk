// vector3.cu
// Include file for chicken, adapted from the Dlang modules geom/vector3 and geom/properties.
// PJ 2022-09-11

#ifndef VECTOR3_INCLUDED
#define VECTOR3_INCLUDED

#include <string>
#include "number.cu"

using namespace std;

struct Vector3 {
    number x;
    number y;
    number z;

    __host__ __device__
    void set(number _x, number _y, number _z=0.0)
    {
        x = _x; y = _y; z = _z;
    }

    __host__ __device__
    void set(const Vector3& other)
    {
        x = other.x; y = other.y; z = other.z;
    }

    string toString() {
        return "[x=" + to_string(x) + ", y=" + to_string(y) + ", z=" + to_string(z) + "]";
    }

    __host__ __device__
    void add(const Vector3& other)
    {
        x += other.x; y += other.y; z += other.z;
    }

    __host__ __device__
    void sub(const Vector3& other)
    {
        x -= other.x; y -= other.y; z -= other.z;
    }

    __host__ __device__
    void mul(number scalar)
    {
        x *= scalar; y *= scalar; z *= scalar;
    }

    __host__ __device__
    void div(number scalar)
    {
        x /= scalar; y /= scalar; z /= scalar;
    }

    __host__ __device__
    number dot(const Vector3& other) const
    {
        return x*other.x + y*other.y + z*other.z;
    }

    __host__ __device__
    void normalize()
    // Scales the vector to unit magnitude.
    {
        number magnitude = sqrt(x*x + y*y + z*z);
        if (magnitude > 0.0) {
            x /= magnitude; y /= magnitude; z /= magnitude;
        } else {
            // Clean up, in case dot product underflows.
            x = y = z = 0.0;
        }
        // Flush small components to zero.
        constexpr double small = 1.0e-30;
        if (fabs(x) < small) { x = 0.0; }
        if (fabs(y) < small) { y = 0.0; }
        if (fabs(z) < small) { z = 0.0; }
    }

    // Transform functions used to reorient vector values in the CFD codes.

    __host__ __device__
    void transform_to_local_frame(const Vector3& n, const Vector3& t1, const Vector3& t2)
    // Rotate v from the global xyz coordinate system into the local frame
    // defined by the orthogonal unit vectors n,t1,t2.
    //
    // We assume, without checking, that these vectors do nicely define
    // such a local system.
    {
        number v_x = dot(n); // normal component
        number v_y = dot(t1); // tangential component 1
        number v_z = dot(t2); // tangential component 2
        x = v_x; y = v_y; z = v_z;
    }

    __host__ __device__
    void transform_to_global_frame(const Vector3& n, const Vector3& t1, const Vector3& t2)
    // Rotate v back into the global (xyz) coordinate system.
    {
        number v_x = x*n.x + y*t1.x + z*t2.x; // global-x
        number v_y = x*n.y + y*t1.y + z*t2.y; // global-y
        number v_z = x*n.z + y*t1.z + z*t2.z; // global-z
        x = v_x; y = v_y; z = v_z;
    }

    // 2D flavour for change of coordinate system functions.

    __host__ __device__
    void transform_to_local_frame(const Vector3& n, const Vector3& t1)
    {
        number v_x = x*n.x + y*n.y;   // normal component
        number v_y = x*t1.x + y*t1.y; // tangential component 1
        x = v_x; y = v_y; z = 0.0;
    }

    __host__ __device__
    void transform_to_global_frame(const Vector3& n, const Vector3& t1)
    // Rotate v back into the global (xy) coordinate system.
    {
        number v_x = x*n.x + y*t1.x; // global-x
        number v_y = x*n.y + y*t1.y; // global-y
        x = v_x; y = v_y; z = 0.0;
    }

}; // end Vector3

__host__ __device__
void cross(Vector3& v3, const Vector3& v1, const Vector3& v2)
// Vector cross product for use in a single statement that will not make temporaries.
{
    v3.x = v1.y * v2.z - v2.y * v1.z;
    v3.y = v2.x * v1.z - v1.x * v2.z;
    v3.z = v1.x * v2.y - v2.x * v1.y;
}

__host__ __device__
Vector3 cross(const Vector3& v1, const Vector3& v2)
{
    Vector3 v3;
    cross(v3, v1, v2);
    return v3;
}

/**
 * Component forms
 */

__host__ __device__
number dot_product(number ax, number ay, number az,
                   number bx, number by, number bz)
{
    return ax*bx + ay*by + az*bz;
}

__host__ __device__
void cross_product(number ax, number ay, number az,
                   number bx, number by, number bz,
                   number& cx, number& cy, number& cz)
{
    cx = ay*bz - az*by;
    cy = az*bx - ax*bz;
    cz = ax*by - ay*bx;
    return;
}

__host__ __device__
number dot(number ax, number ay, number az, const Vector3& other)
{
    return ax*other.x + ay*other.y + az*other.z;
}

__host__ __device__
void quad_properties(const Vector3& p0, const Vector3& p1,
                     const Vector3& p2, const Vector3& p3,
                     Vector3& centroid,
                     Vector3& n, Vector3& t1, Vector3& t2,
                     number& area,
                     number tol=1.0e-12, number area_tol=1.0e-20)
// Quadrilateral properties of centroid, associated unit normals and area.
//   p3-----p2
//   |      |
//   |      |
//   p0-----p1
// Resultant normal vector is up, toward you.
// Assume that all points are in the one plane.
{
    centroid.set(0.25*(p0.x + p1.x + p2.x + p3.x),
                 0.25*(p0.y + p1.y + p2.y + p3.y),
                 0.25*(p0.z + p1.z + p2.z + p3.z));
    // Compute areas via the cross products.
    // Vector3 vector_area = 0.25 * cross(p1-p0+p2-p3, p3-p0+p2-p1);
    number p01x = p1.x-p0.x+p2.x-p3.x;
    number p01y = p1.y-p0.y+p2.y-p3.y;
    number p01z = p1.z-p0.z+p2.z-p3.z;
    number p03x = p3.x-p0.x+p2.x-p1.x;
    number p03y = p3.y-p0.y+p2.y-p1.y;
    number p03z = p3.z-p0.z+p2.z-p1.z;
    number vector_area_x, vector_area_y, vector_area_z;
    cross_product(p01x, p01y, p01z, p03x, p03y, p03z,
                  vector_area_x, vector_area_y, vector_area_z);
    vector_area_x *= 0.25; vector_area_y *= 0.25; vector_area_z *= 0.25;
    // unit-normal and area
    // area = abs(vector_area);
    area = sqrt(vector_area_x*vector_area_x + vector_area_y*vector_area_y + vector_area_z*vector_area_z);
    if (area > area_tol) {
        // n = unit(vector_area);
        n.set(vector_area_x/area, vector_area_y/area, vector_area_z/area);
        // Tangent unit-vectors:
        // t1 is parallel to side01 and side32,
        // t2 is normal to n and t1
        // t1 = unit((p1-p0)+(p2-p3)); // Works even if one edge has zero length.
        t1.set(p01x, p01y, p01z);
        t1.normalize();
        // t2 = unit(cross(n, t1)); // Calling unit() to tighten up the magnitude.
        cross(t2, n, t1);
        t2.normalize();
    } else {
        // We have nothing meaningful to put into the unit vectors,
        // so, put in an arbitrary but orthogonal set.
        n.set(1.0,0.0,0.0); t1.set(0.0,1.0,0.0); t2.set(0.0,0.0,1.0);
    }
} // end quad_properties()


// Because of the way we lose precision when reading and writing files,
// it may be that the vertices are not quite in their ideal position.
// We need a couple of finite, but small, tolerances to deal with
// collapsed volumes.
constexpr number smallButSignificantVolume = 1.0e-12;
constexpr number verySmallVolume = 1.0e-20;

// For the tetrahedron geometry, we consider p0,p1,p2 the base.
// Looking from p3 back toward the base, a counter-clockwise cycle
// p0->p1->p2->p0 gives a positive volume.

number tetrahedron_volume(const Vector3& p0, const Vector3& p1,
                          const Vector3& p2, const Vector3& p3)
{
    // Was return dot(p3-p0, cross(p1-p0, p2-p0)) / 6.0;
    // Now, expand up to components to avoid allocation (in Dlang).
    number dx01=p1.x-p0.x; number dy01=p1.y-p0.y; number dz01=p1.z-p0.z;
    number dx02=p2.x-p0.x; number dy02=p2.y-p0.y; number dz02=p2.z-p0.z;
    number cx, cy, cz;
    cross_product(dx01, dy01, dz01, dx02, dy02, dz02, cx, cy, cz);
    number dx03=p3.x-p0.x; number dy03=p3.y-p0.y; number dz03=p3.z-p0.z;
    return dot_product(dx03, dy03, dz03, cx, cy, cz) / 6.0;
} // end tetrahedron_volume()

void tetrahedron_properties(const Vector3& p0, const Vector3& p1,
                            const Vector3& p2, const Vector3& p3,
                            Vector3& centroid, number& volume)
// Variant without L_min calculation.
{
    centroid.set(0.25*(p0.x + p1.x + p2.x + p3.x),
                 0.25*(p0.y + p1.y + p2.y + p3.y),
                 0.25*(p0.z + p1.z + p2.z + p3.z));
    volume = tetrahedron_volume(p0, p1, p2, p3);
} // end tetrahedron_properties()

number tetragonal_dipyramid_volume(const Vector3& p0, const Vector3& p1,
                                   const Vector3& p2, const Vector3& p3,
                                   const Vector3& pb, const Vector3& pc)
// J. Grandy (1997) Efficient Computation of Volume of Hexahedral Cells UCRL-ID-128886.
// Base of each dipyramid is specified clockwise from the outside.
// pc is apex
// pb is barycentre of base quad.
// A base quad cycle p0->p1->p2->p3->p0 that is counterclockwise when looking
// towards it from pc will result in a positive volume.
// A negative volume indicates that the cycle is clockwise when looking from pc.
{
    // number volume = dot(pc-pb, cross(p1-p0+p2-p3, p3-p0+p2-p1)) / 12.0;
    number p01x = 0.5*(p1.x-p0.x+p2.x-p3.x);
    number p01y = 0.5*(p1.y-p0.y+p2.y-p3.y);
    number p01z = 0.5*(p1.z-p0.z+p2.z-p3.z);
    number p03x = 0.5*(p3.x-p0.x+p2.x-p1.x);
    number p03y = 0.5*(p3.y-p0.y+p2.y-p1.y);
    number p03z = 0.5*(p3.z-p0.z+p2.z-p1.z);
    number vector_area_x, vector_area_y, vector_area_z;
    cross_product(p01x, p01y, p01z, p03x, p03y, p03z,
                  vector_area_x, vector_area_y, vector_area_z);
    number bcx = pc.x-pb.x; number bcy = pc.y-pb.y; number bcz = pc.z-pb.z;
    number volume = dot_product(bcx,bcy,bcz,vector_area_x,vector_area_y,vector_area_z)/3.0;
    return volume;
} // end tetragonal_dipyramid_volume()

void pyramid_properties(const Vector3& p0, const Vector3& p1,
                        const Vector3& p2, const Vector3& p3,
                        const Vector3& p4, bool true_centroid,
                        Vector3& centroid, number& volume)
{
    // p0-p1-p2-p3 is the quadrilateral base, p4 is the peak.
    // cycle of base vertices is counterclockwise, viewed from p4.
    //
    // Split into 4 tetrahedra and sum contributions to volume and moment.
    Vector3 pmB; // Mid-point of quadrilateral base.
    pmB.set(0.25*(p0.x+p1.x+p2.x+p3.x),
            0.25*(p0.y+p1.y+p2.y+p3.y),
            0.25*(p0.z+p1.z+p2.z+p3.z));
    //
    volume = 0.0; Vector3 moment = Vector3{0.0, 0.0, 0.0};
    number tet_volume; Vector3 tet_centroid;
    tetrahedron_properties(p0, p1, pmB, p4, tet_centroid, tet_volume);
    volume += tet_volume; tet_centroid.mul(tet_volume); moment.add(tet_centroid);
    tetrahedron_properties(p1, p2, pmB, p4, tet_centroid, tet_volume);
    volume += tet_volume; tet_centroid.mul(tet_volume); moment.add(tet_centroid);
    tetrahedron_properties(p2, p3, pmB, p4, tet_centroid, tet_volume);
    volume += tet_volume; tet_centroid.mul(tet_volume); moment.add(tet_centroid);
    tetrahedron_properties(p3, p0, pmB, p4, tet_centroid, tet_volume);
    volume += tet_volume; tet_centroid.mul(tet_volume); moment.add(tet_centroid);
    //
    if (fabs(volume) > 0.0) { moment.div(volume); } // to get overall centroid
    if (true_centroid) {
        centroid = moment;
    } else {
        // approximating the centroid via a simple averaging of vertex positions
        // has shown to be more robust when importing an unstructurd grid and
        // also appears to provide a better distribution of points for the
        // least-squares gradient estimation [KAD 2022-08-03].
        centroid.set(0.2*(p0.x+p1.x+p2.x+p3.x+p4.x),
                     0.2*(p0.y+p1.y+p2.y+p3.y+p4.y),
                     0.2*(p0.z+p1.z+p2.z+p3.z+p4.z));
    }
    //
    return;
} // end pyramid_properties()

void hex_cell_properties(const Vector3& p0, const Vector3& p1,
                         const Vector3& p2, const Vector3& p3,
                         const Vector3& p4, const Vector3& p5,
                         const Vector3& p6, const Vector3& p7,
                         bool true_centroid,
                         Vector3& centroid, number& volume,
                         number& iLen, number& jLen, number& kLen)
{
    // PJ 10-Sep-2012
    // When computing the volume of Rolf's thin, warped cells, we have to do
    // something better than splitting our cell into six tetrahedra, so we do that
    // by dividing the hex cell into six tetragonal dipyramids with the original
    // faces as the pyramid bases.
    //
    // Estimate the centroid so that we can use it as the peak
    // of each of the pyramid sub-volumes.
    centroid.set(0.125*(p0.x+p1.x+p2.x+p3.x+p4.x+p5.x+p6.x+p7.x),
                 0.125*(p0.y+p1.y+p2.y+p3.y+p4.y+p5.y+p6.y+p7.y),
                 0.125*(p0.z+p1.z+p2.z+p3.z+p4.z+p5.z+p6.z+p7.z));
    // Mid-points of faces.
    Vector3 pmN;
    pmN.set(0.25*(p3.x+p2.x+p6.x+p7.x),
            0.25*(p3.y+p2.y+p6.y+p7.y),
            0.25*(p3.z+p2.z+p6.z+p7.z));
    Vector3 pmE;
    pmE.set(0.25*(p1.x+p2.x+p6.x+p5.x),
            0.25*(p1.y+p2.y+p6.y+p5.y),
            0.25*(p1.z+p2.z+p6.z+p5.z));
    Vector3 pmS;
    pmS.set(0.25*(p0.x+p1.x+p5.x+p4.x),
            0.25*(p0.y+p1.y+p5.y+p4.y),
            0.25*(p0.z+p1.z+p5.z+p4.z));
    Vector3 pmW;
    pmW.set(0.25*(p0.x+p3.x+p7.x+p4.x),
            0.25*(p0.y+p3.y+p7.y+p4.y),
            0.25*(p0.z+p3.z+p7.z+p4.z));
    Vector3 pmT;
    pmT.set(0.25*(p4.x+p5.x+p6.x+p7.x),
            0.25*(p4.y+p5.y+p6.y+p7.y),
            0.25*(p4.z+p5.z+p6.z+p7.z));
    Vector3 pmB;
    pmB.set(0.25*(p0.x+p1.x+p2.x+p3.x),
            0.25*(p0.y+p1.y+p2.y+p3.y),
            0.25*(p0.z+p1.z+p2.z+p3.z));
    // Lengths between mid-points of faces.
    // Note that we are assuming that the hexahedron is not very skewed
    // when we later use these values as the widths of the hex cell.
    number dx, dy, dz;
    dx = pmE.x - pmW.x; dy = pmE.y - pmW.y; dz = pmE.z - pmW.z;
    iLen = sqrt(dx*dx + dy*dy + dz*dz);
    dx = pmN.x - pmS.x; dy = pmN.y - pmS.y; dz = pmN.z - pmS.z;
    jLen = sqrt(dx*dx + dy*dy + dz*dz);
    dx = pmT.x - pmB.x; dy = pmT.y - pmB.y; dz = pmT.z - pmB.z;
    kLen = sqrt(dx*dx + dy*dy + dz*dz);
    // writeln("Single hexahedron divided into six tetragonal dipyramids.");
    // J. Grandy (1997) Efficient Computation of Volume of Hexahedral Cells UCRL-ID-128886.
    // Base of each dipyramid is specified clockwise from the outside.
    number sub_volume; Vector3 sub_centroid;
    volume = 0.0; Vector3 moment = Vector3{0.0, 0.0, 0.0};
    pyramid_properties(p6, p7, p3, p2, centroid, true, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid.mul(sub_volume); moment.add(sub_centroid);
    pyramid_properties(p5, p6, p2, p1, centroid, true, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid.mul(sub_volume); moment.add(sub_centroid);
    pyramid_properties(p4, p5, p1, p0, centroid, true, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid.mul(sub_volume); moment.add(sub_centroid);
    pyramid_properties(p7, p4, p0, p3, centroid, true, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid.mul(sub_volume); moment.add(sub_centroid);
    pyramid_properties(p7, p6, p5, p4, centroid, true, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid.mul(sub_volume); moment.add(sub_centroid);
    pyramid_properties(p0, p1, p2, p3, centroid, true, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid.mul(sub_volume); moment.add(sub_centroid);
    //
    if ( (volume < 0.0 && fabs(volume) < smallButSignificantVolume) ||
         (volume >= 0.0 && volume < verySmallVolume) ) {
        // We assume that we have a collapsed hex cell;
        // no real problem here but it may be a problem for the client code.
        // That code should test the value of volume, on return.
        volume = 0.0;
    }
    //
    if (fabs(volume) > 0.0) { moment.div(volume); } // to get overall centroid
    if (true_centroid) { centroid = moment; }
    return;
} // end hex_cell_properties()

#endif
