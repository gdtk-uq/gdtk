// properties.d
// Geometric properties of various geometric elements in our Vector3 world.

module geom.elements.properties;

import std.math;
import std.conv;
import geom.geometry_exception;
import geom.elements.nomenclature;
import geom.elements.vector3;
import ntypes.complex;
import nm.number;

@nogc
bool inside_bounding_box(ref const(Vector3) p,
                         ref const(Vector3) bb0, ref const(Vector3) bb1,
                         int dimensions)
{
    // Returns true if we do not determine that the point is outside the box.
    // So it may be on the surface.
    bool is_inside = true;
    if (p.x < bb0.x) is_inside = false;
    if (p.y < bb0.y) is_inside = false;
    if (p.x > bb1.x) is_inside = false;
    if (p.y > bb1.y) is_inside = false;
    if (dimensions == 3) {
        if (p.z < bb0.z) is_inside = false;
        if (p.z > bb1.z) is_inside = false;
    }
    return is_inside;
} // end inside_bounding_box()


//------------------------------------------------------------------------
// Utility functions for cell properties in the finite-volume code.

/** Triangle properties of centroid, associated unit normals and area.
 *  p2
 *   |\
 *   | \
 *   |  \
 *   |   \
 *  p0----p1
 * Resultant normal vector is up, toward you.
 * Assume that all points are in the one plane.
 */
@nogc
void triangle_properties(ref const(Vector3) p0, ref const(Vector3) p1,
                         ref const(Vector3) p2,
                         ref Vector3 centroid,
                         ref Vector3 n, ref Vector3 t1, ref Vector3 t2,
                         ref number area,
                         double tol=1.0e-12, double area_tol=1.0e-20)
{
    centroid.set((p0.x + p1.x + p2.x)/3.0,
                 (p0.y + p1.y + p2.y)/3.0,
                 (p0.z + p1.z + p2.z)/3.0);
    // Compute areas via the cross products.
    number p01x=p1.x-p0.x; number p01y=p1.y-p0.y; number p01z=p1.z-p0.z;
    number p02x=p2.x-p0.x; number p02y=p2.y-p0.y; number p02z=p2.z-p0.z;
    number vector_area_x; number vector_area_y; number vector_area_z;
    // Vector3 vector_area = 0.5 * cross(p1-p0, p2-p0);
    cross_product(p01x, p01y, p01z, p02x, p02y, p02z,
                  vector_area_x, vector_area_y, vector_area_z);
    vector_area_x *= 0.5; vector_area_y *= 0.5; vector_area_z *= 0.5;
    // unit-normal and area
    // area = abs(vector_area);
    area = sqrt(vector_area_x^^2 + vector_area_y^^2 + vector_area_z^^2);
    if (area > area_tol) {
        // n = unit(vector_area);
        n.set(vector_area_x/area, vector_area_y/area, vector_area_z/area);
        // Tangent unit-vectors:
        // t1 is parallel to side01
        // t2 is normal to n and t1
        number abs_p01 = sqrt(p01x^^2 + p01y^^2 + p01z^^2);
        // t1 = unit(p1-p0);
        t1.set(p01x/abs_p01, p01y/abs_p01, p01z/abs_p01);
        // t2 = unit(cross(n, t1)); // Calling unit() to tighten up the magnitude.
        cross(t2, n, t1);
        t2.normalize();
    } else {
        // We have nothing meaningful to put into the unit vectors,
        // so, put in an arbitrary but orthogonal set.
        n.set(1.0,0.0,0.0); t1.set(0.0,1.0,0.0); t2.set(0.0,0.0,1.0);
    }
} // end triangle_properties()

/** Quadrilateral properties of centroid, associated unit normals and area.
 *   p3-----p2
 *   |      |
 *   |      |
 *   p0-----p1
 * Resultant normal vector is up, toward you.
 * Assume that all points are in the one plane.
 */
@nogc
void quad_properties(ref const(Vector3) p0, ref const(Vector3) p1,
                     ref const(Vector3) p2, ref const(Vector3) p3,
                     ref Vector3 centroid,
                     ref Vector3 n, ref Vector3 t1, ref Vector3 t2,
                     ref number area,
                     double tol=1.0e-12, double area_tol=1.0e-20)
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
    area = sqrt(vector_area_x^^2 + vector_area_y^^2 + vector_area_z^^2);
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

@nogc
void quad_properties(ref const(Vector3)[] p, ref Vector3 centroid,
                     ref Vector3 n, ref Vector3 t1, ref Vector3 t2,
                     ref number area,
                     double tol=1.0e-12, double area_tol=1.0e-20)
{
    quad_properties(p[0], p[1], p[2], p[3], centroid, n, t1, t2, area, tol, area_tol);
}

@nogc
number xyplane_area(ref const(Vector3) p0, ref const(Vector3) p1, ref const(Vector3) p2)
// 2
// |\
// | \
// |  \
// 0---1
{
    return 0.5 * ((p2.x + p1.x) * (p2.y - p1.y) +
                  (p0.x + p2.x) * (p0.y - p2.y) +
                  (p1.x + p0.x) * (p1.y - p0.y));
}

@nogc
number xyplane_area(ref const(Vector3) p0, ref const(Vector3) p1,
                    ref const(Vector3) p2, ref const(Vector3) p3)
// 3-----2
// |     |    j
// |  c  |    ^
// |     |    |
// 0-----1    O-->i
{
    return 0.5 * ((p2.x + p1.x) * (p2.y - p1.y) + (p3.x + p2.x) * (p3.y - p2.y) +
                  (p0.x + p3.x) * (p0.y - p3.y) + (p1.x + p0.x) * (p1.y - p0.y));
}

@nogc
void xyplane_quad_cell_properties(ref const(Vector3) p0, ref const(Vector3) p1,
                                  ref const(Vector3) p2, ref const(Vector3) p3,
                                  ref Vector3 centroid, ref number xyplane_area,
                                  ref number iLen, ref number jLen, ref number minLen)
// Cell layout goes back to cns4u notation.
// C-----B     3-----2
// |     |     |     |    j
// |  c  |     |  c  |    ^
// |     |     |     |    |
// D-----A     0-----1    O-->i
{
    // These are the corners.
    number xA = p1.x; number yA = p1.y;
    number xB = p2.x; number yB = p2.y;
    number xC = p3.x; number yC = p3.y;
    number xD = p0.x; number yD = p0.y;
    //
    xyplane_area = 0.5 * ((xB + xA) * (yB - yA) + (xC + xB) * (yC - yB) +
                          (xD + xC) * (yD - yC) + (xA + xD) * (yA - yD));
    //
    number x = 1.0 / (xyplane_area * 6.0) *
        ((yB - yA) * (xA * xA + xA * xB + xB * xB) +
         (yC - yB) * (xB * xB + xB * xC + xC * xC) +
         (yD - yC) * (xC * xC + xC * xD + xD * xD) +
         (yA - yD) * (xD * xD + xD * xA + xA * xA));
    number y = -1.0 / (xyplane_area * 6.0) *
        ((xB - xA) * (yA * yA + yA * yB + yB * yB) +
         (xC - xB) * (yB * yB + yB * yC + yC * yC) +
         (xD - xC) * (yC * yC + yC * yD + yD * yD) +
         (xA - xD) * (yD * yD + yD * yA + yA * yA));
    centroid.set(x, y, to!number(0.0));
    //
    // Check cell length scale using North and East boundaries.
    // Also, save the minimum length for later use in the CFL checking routine.
    number dxN = xC - xB; number dyN = yC - yB;
    number dxE = xA - xB; number dyE = yA - yB;
    number lengthN = sqrt(dxN * dxN + dyN * dyN);
    number lengthE = sqrt(dxE * dxE + dyE * dyE);
    number length_cross = xyplane_area / fmax(lengthN, lengthE);
    // estimate of minimum width of cell
    minLen = fmin(lengthN, lengthE);
    if (length_cross < minLen) minLen = length_cross;

    // Record the cell widths in the i- and j-index directions.
    // The widths are measured between corresponding midpoints of the bounding interfaces.
    // This data is used by the high-order reconstruction.
    number xN = 0.5 * (xC + xB);
    number yN = 0.5 * (yC + yB);
    number xS = 0.5 * (xD + xA);
    number yS = 0.5 * (yD + yA);
    number xE = 0.5 * (xA + xB);
    number yE = 0.5 * (yA + yB);
    number xW = 0.5 * (xD + xC);
    number yW = 0.5 * (yD + yC);
    number dx = xN - xS;
    number dy = yN - yS;
    jLen = sqrt(dx * dx + dy * dy);
    dx = xE - xW;
    dy = yE - yW;
    iLen = sqrt(dx * dx + dy * dy);
} // end xyplane_quad_cell_properties()

@nogc
number xyplane_quad_squareness(ref const(Vector3) p0, ref const(Vector3) p1,
                               ref const(Vector3) p2, ref const(Vector3) p3)
{
    Vector3 centroid;
    number xyplane_area, iLen, jLen, minLen;
    xyplane_quad_cell_properties(p0, p1, p2, p3, centroid, xyplane_area, iLen, jLen, minLen);
    number perimeter = distance_between(p0, p1) + distance_between(p1, p2) +
        distance_between(p2, p3) + distance_between(p3, p0);
    return 16*xyplane_area/(perimeter^^2);
}

@nogc
void xyplane_triangle_cell_properties(ref const(Vector3) p0, ref const(Vector3) p1,
                                      ref const(Vector3) p2,
                                      ref Vector3 centroid, ref number xyplane_area,
                                      ref number iLen, ref number jLen, ref number minLen)
// p2
//  |\
//  | \
//  |  \
// p0--p1  counter-clockwise cycle when looking down at the x,y-plane
{
    // These are the corners.
    number x0 = p0.x; number y0 = p0.y;
    number x1 = p1.x; number y1 = p1.y;
    number x2 = p2.x; number y2 = p2.y;
    //
    xyplane_area = 0.5*((x1+x0)*(y1-y0) + (x2+x1)*(y2-y1) + (x0+x2)*(y0-y2));
    //
    number x = 1.0/(xyplane_area*6.0) *
        ((y1-y0)*(x0*x0 + x0*x1 + x1*x1) +
         (y2-y1)*(x1*x1 + x1*x2 + x2*x2) +
         (y0-y2)*(x2*x2 + x2*x0 + x0*x0));
    number y = -1.0/(xyplane_area*6.0) *
        ((x1-x0)*(y0*y0 + y0*y1 + y1*y1) +
         (x2-x1)*(y1*y1 + y1*y2 + y2*y2) +
         (x0-x2)*(y2*y2 + y2*y0 + y0*y0));
    centroid.set(x, y, to!number(0.0));
    //
    // Also, save the minimum length for later use in the CFL checking routine.
    minLen = sqrt(xyplane_area);
    // i- and j-directions don't have much significance in unstructured grids.
    // Just fill in the data with something
    jLen = minLen;
    iLen = minLen;
    return;
}

// For the tetrahedron geometry, we consider p0,p1,p2 the base.
// Looking from p3 back toward the base, a counter-clockwise cycle
// p0->p1->p2->p0 gives a positive volume.

@nogc
number tetrahedron_volume(ref const(Vector3) p0, ref const(Vector3) p1,
                          ref const(Vector3) p2, ref const(Vector3) p3)
{
    // Was return dot(p3-p0, cross(p1-p0, p2-p0)) / 6.0;
    // Now, expand up to components to avoid allocation.
    number dx01=p1.x-p0.x; number dy01=p1.y-p0.y; number dz01=p1.z-p0.z;
    number dx02=p2.x-p0.x; number dy02=p2.y-p0.y; number dz02=p2.z-p0.z;
    number cx, cy, cz;
    cross_product(dx01, dy01, dz01, dx02, dy02, dz02, cx, cy, cz);
    number dx03=p3.x-p0.x; number dy03=p3.y-p0.y; number dz03=p3.z-p0.z;
    return dot_product(dx03, dy03, dz03, cx, cy, cz) / 6.0;
} // end tetrahedron_volume()

@nogc
void tetrahedron_properties(ref const(Vector3) p0, ref const(Vector3) p1,
                            ref const(Vector3) p2, ref const(Vector3) p3,
                            ref Vector3 centroid, ref number volume)
// Variant without L_min calculation.
{
    centroid.set(0.25*(p0.x + p1.x + p2.x + p3.x),
                 0.25*(p0.y + p1.y + p2.y + p3.y),
                 0.25*(p0.z + p1.z + p2.z + p3.z));
    volume = tetrahedron_volume(p0, p1, p2, p3);
} // end tetrahedron_properties()

@nogc
void tetrahedron_properties(ref const(Vector3)[] p,
                            ref Vector3 centroid, ref number volume)
{
    tetrahedron_properties(p[0], p[1], p[2], p[3], centroid, volume);
}

@nogc
void tetrahedron_properties(ref const(Vector3) p0, ref const(Vector3) p1,
                            ref const(Vector3) p2, ref const(Vector3) p3,
                            ref Vector3 centroid, ref number volume, ref number L_min)
// Variant including L_min calculation.
{
    centroid.set(0.25*(p0.x + p1.x + p2.x + p3.x),
                 0.25*(p0.y + p1.y + p2.y + p3.y),
                 0.25*(p0.z + p1.z + p2.z + p3.z));
    volume = tetrahedron_volume(p0, p1, p2, p3);
    number third = 1.0/3.0;
    L_min = pow(volume, third);
    Vector3 pmid = Vector3(third*(p0.x+p1.x+p2.x), third*(p0.y+p1.y+p2.y), third*(p0.z+p1.z+p2.z));
    L_min = fmin(L_min, 2.0*distance_between(pmid, centroid));
    pmid = Vector3(third*(p0.x+p1.x+p3.x), third*(p0.y+p1.y+p3.y), third*(p0.z+p1.z+p3.z));
    L_min = fmin(L_min, 2.0*distance_between(pmid, centroid));
    pmid = Vector3(third*(p1.x+p2.x+p3.x), third*(p1.y+p2.y+p3.y), third*(p1.z+p2.z+p3.z));
    L_min = fmin(L_min, 2.0*distance_between(pmid, centroid));
    pmid = Vector3(third*(p2.x+p0.x+p3.x), third*(p2.y+p0.y+p3.y), third*(p2.z+p0.z+p3.z));
    L_min = fmin(L_min, 2.0*distance_between(pmid, centroid));
} // end tetrahedron_properties()

@nogc
void tetrahedron_properties(ref const(Vector3)[] p,
                            ref Vector3 centroid, ref number volume, ref number L_min)
{
    tetrahedron_properties(p[0], p[1], p[2], p[3], centroid, volume, L_min);
}

// Because of the way we lose precision when reading and writing files,
// it may be that the vertices are not quite in their ideal position.
// We need a couple of finite, but small, tolerances to deal with
// collapsed volumes.
double smallButSignificantVolume = 1.0e-12;
double verySmallVolume = 1.0e-20;


@nogc
number tetragonal_dipyramid_volume(ref const(Vector3) p0, ref const(Vector3) p1,
                                   ref const(Vector3) p2, ref const(Vector3) p3,
                                   ref const(Vector3) pb, ref const(Vector3) pc)
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

@nogc
void pyramid_properties(ref const(Vector3) p0, ref const(Vector3) p1,
                        ref const(Vector3) p2, ref const(Vector3) p3,
                        ref const(Vector3) p4, bool true_centroid,
                        ref Vector3 centroid, ref number volume)
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
    volume = 0.0; Vector3 moment = Vector3(0.0, 0.0, 0.0);
    number tet_volume; Vector3 tet_centroid;
    tetrahedron_properties(p0, p1, pmB, p4, tet_centroid, tet_volume);
    volume += tet_volume; tet_centroid *= tet_volume; moment.add(tet_centroid);
    tetrahedron_properties(p1, p2, pmB, p4, tet_centroid, tet_volume);
    volume += tet_volume; tet_centroid *= tet_volume; moment.add(tet_centroid);
    tetrahedron_properties(p2, p3, pmB, p4, tet_centroid, tet_volume);
    volume += tet_volume; tet_centroid *= tet_volume; moment.add(tet_centroid);
    tetrahedron_properties(p3, p0, pmB, p4, tet_centroid, tet_volume);
    volume += tet_volume; tet_centroid *= tet_volume; moment.add(tet_centroid);
    //
    if (fabs(volume) > 0.0) { moment /= volume; } // to get overall centroid
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

@nogc
void pyramid_properties(ref const(Vector3)[] p, bool true_centroid,
                        ref Vector3 centroid, ref number volume)
{
    pyramid_properties(p[0], p[1], p[2], p[3], p[4], true_centroid, centroid, volume);
}

@nogc
void pyramid_properties(ref const(Vector3) p0, ref const(Vector3) p1,
                        ref const(Vector3) p2, ref const(Vector3) p3,
                        ref const(Vector3) p4, bool true_centroid,
                        ref Vector3 centroid, ref number volume, ref number L_min)
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
    volume = 0.0; Vector3 moment = Vector3(0.0, 0.0, 0.0);
    number tet_volume; Vector3 tet_centroid;
    tetrahedron_properties(p0, p1, pmB, p4, tet_centroid, tet_volume);
    volume += tet_volume; tet_centroid *= tet_volume; moment.add(tet_centroid);
    tetrahedron_properties(p1, p2, pmB, p4, tet_centroid, tet_volume);
    volume += tet_volume; tet_centroid *= tet_volume; moment.add(tet_centroid);
    tetrahedron_properties(p2, p3, pmB, p4, tet_centroid, tet_volume);
    volume += tet_volume; tet_centroid *= tet_volume; moment.add(tet_centroid);
    tetrahedron_properties(p3, p0, pmB, p4, tet_centroid, tet_volume);
    volume += tet_volume; tet_centroid *= tet_volume; moment.add(tet_centroid);
    //
    if (fabs(volume) > 0.0) { moment /= volume; } // to get overall centroid
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
    double third = 1.0/3.0;
    L_min = fmin(pow(volume, third), 2.0*distance_between(pmB, centroid));
    Vector3 pmid = Vector3(third*(p0.x+p1.x+p4.x), third*(p0.y+p1.y+p4.y), third*(p0.z+p1.z+p4.z));
    L_min = fmin(L_min, 2.0*distance_between(pmid, centroid));
    pmid = Vector3(third*(p1.x+p2.x+p4.x), third*(p1.y+p2.y+p4.y), third*(p1.z+p2.z+p4.z));
    L_min = fmin(L_min, 2.0*distance_between(pmid, centroid));
    pmid = Vector3(third*(p2.x+p3.x+p4.x), third*(p2.y+p3.y+p4.y), third*(p2.z+p3.z+p4.z));
    L_min = fmin(L_min, 2.0*distance_between(pmid, centroid));
    pmid = Vector3(third*(p3.x+p0.x+p4.x), third*(p3.y+p0.y+p4.y), third*(p3.z+p0.z+p4.z));
    L_min = fmin(L_min, 2.0*distance_between(pmid, centroid));
    //
    return;
} // end pyramid_properties()

@nogc
void pyramid_properties(ref const(Vector3)[] p, bool true_centroid,
                        ref Vector3 centroid, ref number volume, ref number L_min)
{
    pyramid_properties(p[0], p[1], p[2], p[3], p[4], true_centroid, centroid, volume, L_min);
}

@nogc
void wedge_properties(ref const(Vector3) p0, ref const(Vector3) p1,
                      ref const(Vector3) p2, ref const(Vector3) p3,
                      ref const(Vector3) p4, ref const(Vector3) p5,
                      bool true_centroid,
                      ref Vector3 centroid, ref number volume)
{
    // Use the average of the vertex points to get a rough centroid of the wedge.
    centroid.set(1.0/6.0*(p0.x+p1.x+p2.x+p3.x+p4.x+p5.x),
                 1.0/6.0*(p0.y+p1.y+p2.y+p3.y+p4.y+p5.y),
                 1.0/6.0*(p0.z+p1.z+p2.z+p3.z+p4.z+p5.z));
    // Split the wedge into three pyramids and two tetrahedra
    // using this centroid as the peak of each sub-volume.
    number sub_volume; Vector3 sub_centroid;
    volume = 0.0; Vector3 moment = Vector3(0.0, 0.0, 0.0);
    pyramid_properties(p3, p5, p2, p0, centroid, true, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
    pyramid_properties(p1, p2, p5, p4, centroid, true, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
    pyramid_properties(p0, p1, p4, p3, centroid, true, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
    tetrahedron_properties(p0, p2, p1, centroid, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
    tetrahedron_properties(p3, p4, p5, centroid, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
    //
    if (fabs(volume) > 0.0) { moment /= volume; } // to get overall centroid
    if (true_centroid) { centroid = moment; }
    //
    return;
} // end wedge_properties()

@nogc
void wedge_properties(ref const(Vector3)[] p, bool true_centroid,
                      ref Vector3 centroid, ref number volume)
{
    wedge_properties(p[0], p[1], p[2], p[3], p[4], p[5], true_centroid, centroid, volume);
}

@nogc
void wedge_properties(ref const(Vector3) p0, ref const(Vector3) p1,
                      ref const(Vector3) p2, ref const(Vector3) p3,
                      ref const(Vector3) p4, ref const(Vector3) p5,
                      bool true_centroid,
                      ref Vector3 centroid, ref number volume, ref number L_min)
{
    // Use the average of the vertex points to get a rough centroid of the wedge.
    centroid.set(1.0/6.0*(p0.x+p1.x+p2.x+p3.x+p4.x+p5.x),
                 1.0/6.0*(p0.y+p1.y+p2.y+p3.y+p4.y+p5.y),
                 1.0/6.0*(p0.z+p1.z+p2.z+p3.z+p4.z+p5.z));
    // Split the wedge into three pyramids and two tetrahedra
    // using this centroid as the peak of each sub-volume.
    number sub_volume; Vector3 sub_centroid;
    volume = 0.0; Vector3 moment = Vector3(0.0, 0.0, 0.0);
    pyramid_properties(p3, p5, p2, p0, centroid, true, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
    pyramid_properties(p1, p2, p5, p4, centroid, true, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
    pyramid_properties(p0, p1, p4, p3, centroid, true, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
    tetrahedron_properties(p0, p2, p1, centroid, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
    tetrahedron_properties(p3, p4, p5, centroid, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
    //
    if (fabs(volume) > 0.0) { moment /= volume; } // to get overall centroid
    if (true_centroid) { centroid = moment; }
    //
    double third = 1.0/3.0;
    L_min = pow(volume, third);
    Vector3 pmid = Vector3(third*(p0.x+p2.x+p1.x), third*(p0.y+p2.y+p1.y), third*(p0.z+p2.z+p1.z));
    L_min = fmin(L_min, 2.0*distance_between(pmid, centroid));
    pmid = Vector3(third*(p3.x+p4.x+p5.x), third*(p3.y+p4.y+p5.y), third*(p3.z+p4.z+p5.z));
    L_min = fmin(L_min, 2.0*distance_between(pmid, centroid));
    pmid = Vector3(0.25*(p3.x+p5.x+p2.x+p0.x), 0.25*(p3.y+p5.y+p2.y+p0.y), 0.25*(p3.z+p5.z+p2.z+p0.z));
    L_min = fmin(L_min, 2.0*distance_between(pmid, centroid));
    pmid = Vector3(0.25*(p1.x+p2.x+p5.x+p4.x), 0.25*(p1.y+p2.y+p5.y+p4.y), 0.25*(p1.z+p2.z+p5.z+p4.z));
    L_min = fmin(L_min, 2.0*distance_between(pmid, centroid));
    pmid = Vector3(0.25*(p0.x+p1.x+p4.x+p3.x), 0.25*(p0.y+p1.y+p4.y+p3.y), 0.25*(p0.z+p1.z+p4.z+p3.z));
    L_min = fmin(L_min, 2.0*distance_between(pmid, centroid));
    //
    return;
} // end wedge_properties()

@nogc
void wedge_properties(ref const(Vector3)[] p, bool true_centroid,
                      ref Vector3 centroid, ref number volume, ref number L_min)
{
    wedge_properties(p[0], p[1], p[2], p[3], p[4], p[5], true_centroid, centroid, volume, L_min);
}

@nogc
void hex_cell_properties(ref const(Vector3) p0, ref const(Vector3) p1,
                         ref const(Vector3) p2, ref const(Vector3) p3,
                         ref const(Vector3) p4, ref const(Vector3) p5,
                         ref const(Vector3) p6, ref const(Vector3) p7,
                         bool true_centroid,
                         ref Vector3 centroid, ref number volume,
                         ref number iLen, ref number jLen, ref number kLen)
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
    iLen = sqrt(dx^^2 + dy^^2 + dz^^2);
    dx = pmN.x - pmS.x; dy = pmN.y - pmS.y; dz = pmN.z - pmS.z;
    jLen = sqrt(dx^^2 + dy^^2 + dz^^2);
    dx = pmT.x - pmB.x; dy = pmT.y - pmB.y; dz = pmT.z - pmB.z;
    kLen = sqrt(dx^^2 + dy^^2 + dz^^2);
    // writeln("Single hexahedron divided into six tetragonal dipyramids.");
    // J. Grandy (1997) Efficient Computation of Volume of Hexahedral Cells UCRL-ID-128886.
    // Base of each dipyramid is specified clockwise from the outside.
    number sub_volume; Vector3 sub_centroid;
    volume = 0.0; Vector3 moment = Vector3(0.0, 0.0, 0.0);
    pyramid_properties(p6, p7, p3, p2, centroid, true, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
    pyramid_properties(p5, p6, p2, p1, centroid, true, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
    pyramid_properties(p4, p5, p1, p0, centroid, true, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
    pyramid_properties(p7, p4, p0, p3, centroid, true, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
    pyramid_properties(p7, p6, p5, p4, centroid, true, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
    pyramid_properties(p0, p1, p2, p3, centroid, true, sub_centroid, sub_volume);
    volume += sub_volume; sub_centroid *= sub_volume; moment.add(sub_centroid);
    //
    if ( (volume < 0.0 && fabs(volume) < smallButSignificantVolume) ||
         (volume >= 0.0 && volume < verySmallVolume) ) {
        // We assume that we have a collapsed hex cell;
        // no real problem here but it may be a problem for the client code.
        // That code should test the value of volume, on return.
        volume = 0.0;
    }
    //
    if (fabs(volume) > 0.0) { moment /= volume; } // to get overall centroid
    if (true_centroid) { centroid = moment; }
    return;
} // end hex_cell_properties()

@nogc
void hex_cell_properties(ref const(Vector3)[] p, bool true_centroid,
                         ref Vector3 centroid, ref number volume,
                         ref number iLen, ref number jLen, ref number kLen)
{
    hex_cell_properties(p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7],
                        true_centroid, centroid, volume, iLen, jLen, kLen);
}


//------------------------------------------------------------------------
// Utility functions for searching cells in the finite-volume code.

@nogc
int inside_triangle(ref const(Vector3) p, ref const(Vector3) a,
                    ref const(Vector3) b, ref const(Vector3) c)
{
    // Vector3 n = cross(a-c, b-c); n.normalize(); // normal to plane of triangle
    Vector3 ca; ca.set(a); ca.add(c, -1.0);
    Vector3 cb; cb.set(b); cb.add(c, -1.0);
    Vector3 n; cross(n, ca, cb); n.normalize();
    // Vector3 varea = cross(p-a, p-b);
    Vector3 ap; ap.set(p); ap.add(a, -1.0);
    Vector3 bp; bp.set(p); bp.add(b, -1.0);
    Vector3 varea; cross(varea, ap, bp);
    number area1 = 0.5 * dot(varea, n); // projected area of triangle pab
    Vector3 cp; cp.set(p); cp.add(c, -1.0);
    cross(varea, bp, cp);
    number area2 = 0.5 * dot(varea, n);
    cross(varea, cp, ap);
    number area3 = 0.5 * dot(varea, n);
    // Only a point inside the triangle will have all areas positive.
    if ( area1 > 0.0 && area2 > 0.0 && area3 > 0.0 ) return 1;
    // However, the point may be damned close to an edge but
    // missing because of floating-point round-off.
    double tol = 1.0e-12;
    if ( (fabs(area1) < tol && area2 > -tol && area3 > -tol) ||
         (fabs(area2) < tol && area1 > -tol && area3 > -tol) ||
         (fabs(area3) < tol && area1 > -tol && area2 > -tol) ) return 2;
    // Otherwise, the point is outside the triangle.
    return 0;
}

// Functions for xy-plane cells.
//
// We split the x,y-plane into half-planes and check which side p is on.
// We also assume that the vertices are numbered counter-clockwise
// around the edges of the polygon.  z-components are ignored.

@nogc
bool on_left_of_xy_line(ref const(Vector3) a, ref const(Vector3) b,
                        ref const(Vector3) p)
// Returns true if p is in the left half-plane (or on the line)
// as we look along the line from a to b.
{
    return (p.x - b.x) * (a.y - b.y) >= (p.y - b.y) * (a.x - b.x);
}

@nogc
bool inside_xy_polygon(ref const(Vector3[]) vtx, ref const(Vector3) p)
// Returns true if the point p is inside or on the polygon boundary.
{
    uint count_on_left = 0;
    foreach (i; 0 .. vtx.length) {
        size_t ip1 = i+1;
        if (ip1 == vtx.length) ip1 = 0; // wrap around
        if (on_left_of_xy_line(vtx[i], vtx[ip1], p)) count_on_left += 1;
    }
    return (count_on_left == vtx.length);
} // end inside_xy_polygon()

@nogc
bool inside_xy_triangle(ref const(Vector3) p0, ref const(Vector3) p1,
                        ref const(Vector3) p2, ref const(Vector3) p)
// Returns true if the point p is inside or on the triangle boundary.
//
// This specialized function avoids the array and loop used
// for the more general polygon.
{
    uint count_on_left = 0;
    if (on_left_of_xy_line(p0, p1, p)) count_on_left += 1;
    if (on_left_of_xy_line(p1, p2, p)) count_on_left += 1;
    if (on_left_of_xy_line(p2, p0, p)) count_on_left += 1;
    return (count_on_left == 3);
} // end inside_xy_triangle()

@nogc
bool inside_xy_quad(ref const(Vector3) p0, ref const(Vector3) p1,
                    ref const(Vector3) p2, ref const(Vector3) p3,
                    ref const(Vector3) p)
// Returns true if the point p is inside or on the quadrangle boundary.
{
    uint count_on_left = 0;
    if (on_left_of_xy_line(p0, p1, p)) count_on_left += 1;
    if (on_left_of_xy_line(p1, p2, p)) count_on_left += 1;
    if (on_left_of_xy_line(p2, p3, p)) count_on_left += 1;
    if (on_left_of_xy_line(p3, p0, p)) count_on_left += 1;
    return (count_on_left == 4);
} // end inside_xy_quad()

@nogc
int winding_number(Vector3[] vertices, ref const(Vector3) p)
// Returns a zero winding number if point p is outside the closed polygon
// defined by the vertices.
// The closing segment is from vertices[$-1] to vertices[0].
{
    int wn = 0;
    Vector3 pr = Vector3(p.x+1.0, p.y);
    foreach (i; 0 .. vertices.length) {
        // Look at each segment of the polygon.
        Vector3* a = (i == 0) ? &(vertices[$-1]) : &(vertices[i-1]);
        Vector3* b = &(vertices[i]);
        if (on_left_of_xy_line(p, pr, *b) &&
            !on_left_of_xy_line(p, pr, *a) &&
            on_left_of_xy_line(*a, *b, p) ) { ++wn; }
        if (!on_left_of_xy_line(p, pr, *b) &&
            on_left_of_xy_line(p, pr, *a) &&
            !on_left_of_xy_line(*a, *b, p) ) { --wn; }
    }
    return wn;
} // end winding_number()

// Functions for 3D cells.

@nogc
bool inside_tetrahedron(ref const(Vector3) p0, ref const(Vector3) p1,
                        ref const(Vector3) p2, ref const(Vector3) p3,
                        ref const(Vector3) p)
// Returns true is the point p is inside or on the tetrahedron surface.
{
    if ((tetrahedron_volume(p0, p2, p1, p)) > 0.0) return false;
    if ((tetrahedron_volume(p0, p1, p3, p)) > 0.0) return false;
    if ((tetrahedron_volume(p1, p2, p3, p)) > 0.0) return false;
    if ((tetrahedron_volume(p2, p0, p3, p)) > 0.0) return false;
    // If we arrive here, we haven't determined that the point is outside...
    return true;
} // end inside_tetrahedron()

@nogc
bool inside_pyramid(ref const(Vector3) p0, ref const(Vector3) p1,
                    ref const(Vector3) p2, ref const(Vector3) p3,
                    ref const(Vector3) p4, ref const(Vector3) p)
// Returns true is the point p is inside or on the pyramid surface.
//
// The test consists of using the quadrilateral faces as the bases of the pyramid
// with the sample point p as the apex.
// The cycle of vertices defining each base is such that the base normal
// will be point out of the cell.
// If any of the pyramid volumes are positive (i.e. p is on the positive
// side of a face) and we assume a convex cell, it means that the point
// is outside the cell and we may say so without further testing.
//
// If the logic looks a bit round-about, there is history, since this function
// was adapted from the inside_hexagon function, below.
{
    // Mid-points of quadrilateral base.
    Vector3 pmid; pmid.set(p0); pmid.add(p1); pmid.add(p2); pmid.add(p3); pmid.scale(0.25);
    // Test the volume of the single pyramid.
    if (tetragonal_dipyramid_volume(p3, p2, p1, p0, pmid, p) > 0.0) return false; // Bottom
    // If we arrive here, we haven't determined that the point is outside...
    // And the tetrahedra formed with the triangular faces.
    if ((tetrahedron_volume(p0, p1, p4, p)) > 0.0) return false;
    if ((tetrahedron_volume(p1, p2, p4, p)) > 0.0) return false;
    if ((tetrahedron_volume(p2, p3, p4, p)) > 0.0) return false;
    if ((tetrahedron_volume(p3, p0, p4, p)) > 0.0) return false;
    return true;
} // end inside_pyramid()

@nogc
bool inside_wedge(ref const(Vector3) p0, ref const(Vector3) p1,
                  ref const(Vector3) p2, ref const(Vector3) p3,
                  ref const(Vector3) p4, ref const(Vector3) p5,
                  ref const(Vector3) p)
// Returns true is the point p is inside or on the wedge surface.
//
// The test consists of using the 3 quadrilateral faces as the bases of pyramids
// with the sample point p as the apex of each.
// The cycle of vertices defining each base is such that the base normal
// will be point out of the cell.
// If any of the pyramid volumes are positive (i.e. p is on the positive
// side of a face) and we assume a convex cell, it means that the point
// is outside the cell and we may say so without further testing.
{
    // Mid-points of quadrilateral faces.
    double quarter = 0.25;
    Vector3 pmA; pmA.set(p0); pmA.add(p3); pmA.add(p4); pmA.add(p1); pmA.scale(quarter); // 0.25(p0+p3+p4+p1)
    Vector3 pmB; pmB.set(p2); pmB.add(p1); pmB.add(p4); pmB.add(p5); pmB.scale(quarter); // 0.25(p2+p1+p4+p5);
    Vector3 pmC; pmC.set(p0); pmC.add(p2); pmC.add(p5); pmC.add(p3); pmC.scale(quarter); // 0.25(p0+p2+p5+p3);
    // Test the volume of each pyramid.
    if (tetragonal_dipyramid_volume(p0, p3, p4, p1, pmA, p) > 0.0) return false;
    if (tetragonal_dipyramid_volume(p2, p1, p4, p5, pmB, p) > 0.0) return false;
    if (tetragonal_dipyramid_volume(p0, p2, p5, p3, pmC, p) > 0.0) return false;
    // And the tetrahedra formed with the triangular faces.
    if ((tetrahedron_volume(p0, p1, p2, p)) > 0.0) return false;
    if ((tetrahedron_volume(p3, p5, p4, p)) > 0.0) return false;
    // If we arrive here, we haven't determined that the point is outside...
    return true;
} // end inside_wedge()

@nogc
bool inside_hexahedron(ref const(Vector3) p0, ref const(Vector3) p1,
                       ref const(Vector3) p2, ref const(Vector3) p3,
                       ref const(Vector3) p4, ref const(Vector3) p5,
                       ref const(Vector3) p6, ref const(Vector3) p7,
                       ref const(Vector3) p)
// Returns true is the point p is inside or on the hexahedron surface.
//
// The test consists of using the 6 cell faces as the bases of pyramids
// with the sample point p as the apex of each.
// The cycle of vertices defining each base is such that the base normal
// will be point out of the cell.
// If any of the pyramid volumes are positive (i.e. p is on the positive
// side of a face) and we assume a convex cell, it means that the point
// is outside the cell and we may say so without further testing.
{
    // Mid-points of faces.
    double quarter = 0.25;
    Vector3 pmN; pmN.set(p3); pmN.add(p2); pmN.add(p6); pmN.add(p7); pmN.scale(quarter); // 0.25(p3+p2+p6+p7)
    Vector3 pmE; pmE.set(p1); pmE.add(p2); pmE.add(p6); pmE.add(p5); pmE.scale(quarter); // 0.25(p1+p2+p6+p5)
    Vector3 pmS; pmS.set(p0); pmS.add(p1); pmS.add(p5); pmS.add(p4); pmS.scale(quarter); // 0.25(p0+p1+p5+p4)
    Vector3 pmW; pmW.set(p0); pmW.add(p3); pmW.add(p7); pmW.add(p4); pmW.scale(quarter); // 0.25(p0+p3+p7+p4)
    Vector3 pmT; pmT.set(p4); pmT.add(p5); pmT.add(p6); pmT.add(p7); pmT.scale(quarter); // 0.25(p4+p5+p6+p7)
    Vector3 pmB; pmB.set(p0); pmB.add(p1); pmB.add(p2); pmB.add(p3); pmB.scale(quarter); // 0.25(p0+p1+p2+p3)
    // Test the volume of each pyramid.
    if (tetragonal_dipyramid_volume(p2, p3, p7, p6, pmN, p) > 0.0) return false; // North
    if (tetragonal_dipyramid_volume(p1, p2, p6, p5, pmE, p) > 0.0) return false; // East
    if (tetragonal_dipyramid_volume(p0, p1, p5, p4, pmS, p) > 0.0) return false; // South
    if (tetragonal_dipyramid_volume(p3, p0, p4, p7, pmW, p) > 0.0) return false; // West
    if (tetragonal_dipyramid_volume(p4, p5, p6, p7, pmT, p) > 0.0) return false; // Top
    if (tetragonal_dipyramid_volume(p1, p0, p3, p2, pmB, p) > 0.0) return false; // Bottom
    // If we arrive here, we haven't determined that the point is outside...
    return true;
} // end inside_hexahedron()


version(properties_test) {
    import util.msg_service;
    int main() {
        Vector3 bb0 = Vector3(1.0, 1.0, 1.0);
        Vector3 bb1 = Vector3(2.0, 2.0, 2.0);
        Vector3 p = to!number(0.5)*(bb0 + bb1);
        assert(inside_bounding_box(p, bb0, bb1, 3), failedUnitTest());
        p += Vector3(0.0, 0.0, 1.0);
        assert(!inside_bounding_box(p, bb0, bb1, 3), failedUnitTest());
        assert(inside_bounding_box(p, bb0, bb1, 2), failedUnitTest());

        Vector3 a = Vector3(1.0, 0.0, 0.0); // plane through a,b,c
        Vector3 b = Vector3(1.0, 1.0, 0.0);
        Vector3 c = Vector3(0.5, 0.0, 0.0);
        Vector3 d = Vector3(0.65, 0.0, 0.0);
        assert(inside_triangle(d, a, b, c) > 0, failedUnitTest());
        Vector3 e = Vector3(0.65, -0.1, 0.0);
        assert(!inside_triangle(e, a, b, c), failedUnitTest());

        auto p0 = Vector3(0.0, 0.0);
        auto p1 = Vector3(1.0, -0.1);
        auto p2 = Vector3(1.0, 1.0);
        auto p3 = Vector3(0.0, 1.0);
        Vector3[] poly = [p0, p1, p2, p3];
        assert(inside_xy_polygon(poly, d), failedUnitTest());
        assert(!inside_xy_polygon(poly, e), failedUnitTest());
        assert(inside_xy_triangle(p0, p1, p2, d), failedUnitTest());
        assert(!inside_xy_triangle(p0, p1, p2, e), failedUnitTest());
        assert(inside_xy_quad(p0, p1, p2, p3, d), failedUnitTest());
        assert(winding_number(poly, d) != 0, failedUnitTest());
        assert(!inside_xy_quad(p0, p1, p2, p3, e), failedUnitTest());
        assert(winding_number(poly, e) == 0, failedUnitTest());

        // Set up a simple triangle for which we can easily compute
        // the properties by hand.
        //      p2
        //    / |
        //   /  |
        // p0---p1
        double alpha = 1.0;
        p0 = Vector3(0.0, 0.0, 1.0);
        p1 = Vector3(alpha, 0.0, 1.0);
        p2 = Vector3(alpha, alpha, 1.0);
        Vector3 centroid, n, t1, t2;
        number area;
        triangle_properties(p0, p1, p2, centroid, n, t1, t2, area);
        assert(approxEqualNumbers(area, to!number(0.5)), failedUnitTest());
        assert(approxEqualNumbers(xyplane_area(p0, p1, p2), to!number(0.5)), failedUnitTest());
        assert(approxEqualVectors(centroid, Vector3(2.0/3,1.0/3,1.0)), failedUnitTest());
        assert(approxEqualVectors(n, Vector3(0.0,0.0,1.0)), failedUnitTest());
        assert(approxEqualVectors(t1, Vector3(1.0,0.0,0.0)), failedUnitTest());
        assert(approxEqualVectors(t2, Vector3(0.0,1.0,0.0)), failedUnitTest());

        version(complex_numbers) {
            // triangle properties - complex derivative test
            double h = 1.0e-20;
            number ih = complex(0,h);
            number zero = 0.0;
            number one = 1.0;
            auto p0_dash = Vector3(zero, zero, one);
            auto p1_dash = Vector3(alpha+ih, zero, one);
            auto p2_dash = Vector3(alpha+ih, alpha+ih, one);
            Vector3 centroid_dash, n_dash, t1_dash, t2_dash;
            number area_dash;
            triangle_properties(p0_dash, p1_dash, p2_dash, centroid_dash,
                                n_dash, t1_dash, t2_dash, area_dash);
            double darea_da = area_dash.im / h;
            Vector3 dcentroid_da = Vector3(centroid_dash.x.im/h, centroid_dash.y.im/h, centroid_dash.z.im/h);
            Vector3 dn_da = Vector3(n_dash.x.im/h, n_dash.y.im/h, n_dash.z.im/h);
            Vector3 dt1_da = Vector3(t1_dash.x.im/h, t1_dash.y.im/h, t1_dash.z.im/h);
            Vector3 dt2_da = Vector3(t2_dash.x.im/h, t2_dash.y.im/h, t2_dash.z.im/h);
            assert(approxEqualVectors(dcentroid_da, Vector3(2.0/3.0*alpha, 1.0/3.0*alpha, 0.0)),
                   failedUnitTest());
            // Because a change in alpha does not change the orientation of the triangle,
            // we expect the unit vectors to have zero sensitivity.
            // [TODO] It would probably be good to set up another case in which alpha
            // does change the orientation of the triangle.
            assert(approxEqualVectors(dn_da, Vector3(0.0, 0.0, 0.0)), failedUnitTest());
            assert(approxEqualVectors(dt1_da, Vector3(0.0, 0.0, 0.0)), failedUnitTest());
            assert(approxEqualVectors(dt2_da, Vector3(0.0, 0.0, 0.0)), failedUnitTest());
            // Since, for this triangle, area = 0.5*alpha^^2, we expect darea/dalpha = alpha
            assert(isClose(darea_da, alpha), failedUnitTest());
        }

        // Set up a simple quadrilateral for which we can easily compute
        // the properties by hand.
        // p3---p2
        // |    |
        // |    |
        // p0---p1
        // For the quadrilateral, we reuse the triangle points but need a fourth point.
        p3 = Vector3(0.0, alpha, 1.0);
        quad_properties(p0, p1, p2, p3, centroid, n, t1, t2, area);
        assert(approxEqualNumbers(area, to!number(1.0)), failedUnitTest());
        assert(approxEqualVectors(centroid, Vector3(0.5,0.5,1.0)), failedUnitTest());
        assert(approxEqualVectors(n, Vector3(0.0,0.0,1.0)), failedUnitTest());
        assert(approxEqualVectors(t1, Vector3(1.0,0.0,0.0)), failedUnitTest());
        assert(approxEqualVectors(t2, Vector3(0.0,1.0,0.0)), failedUnitTest());

        version(complex_numbers) {
            // quad properties - complex derivative test
            auto p3_dash = Vector3(zero, alpha+ih, one);
            quad_properties(p0_dash, p1_dash, p2_dash, p3_dash, centroid_dash,
                            n_dash, t1_dash, t2_dash, area_dash);
            darea_da = area_dash.im / h;
            dcentroid_da = Vector3(centroid_dash.x.im/h, centroid_dash.y.im/h, centroid_dash.z.im/h);
            dn_da = Vector3(n_dash.x.im/h, n_dash.y.im/h, n_dash.z.im/h);
            dt1_da = Vector3(t1_dash.x.im/h, t1_dash.y.im/h, t1_dash.z.im/h);
            dt2_da = Vector3(t2_dash.x.im/h, t2_dash.y.im/h, t2_dash.z.im/h);
            assert(approxEqualVectors(dcentroid_da, Vector3(0.5*alpha, 0.5*alpha, 0.0)),
                   failedUnitTest());
            // Because a change in alpha does not change the orientation of the triangle,
            // we expect the unit vectors to have zero sensitivity.
            // [TODO] It would probably be good to set up another case in which alpha
            // does change the orientation of the triangle.
            assert(approxEqualVectors(dn_da, Vector3(0.0, 0.0, 0.0)), failedUnitTest());
            assert(approxEqualVectors(dt1_da, Vector3(0.0, 0.0, 0.0)), failedUnitTest());
            assert(approxEqualVectors(dt2_da, Vector3(0.0, 0.0, 0.0)), failedUnitTest());
            // Since, for this triangle, area = 0.5*alpha^^2, we expect darea/dalpha = alpha
            assert(isClose(darea_da, 2.0*alpha), failedUnitTest());
        }

        // Build tetrahedron with equilateral triangle (side 1.0) base on xy plane.
        p0 = Vector3(0, 0, 0);
        p1 = Vector3(cos(radians(30)), sin(radians(30)), 0.0);
        p2 = Vector3(0.0, 1.0, 0.0);
        double dx = 0.5 * tan(radians(30));
        double dL = cos(radians(30));
        double dz = sqrt(dL*dL - dx*dx);
        p3 = Vector3(dx, 0.5, dz);
        number volume;
        tetrahedron_properties(p0, p1, p2, p3, centroid, volume);
        assert(approxEqualVectors(centroid, Vector3(dx,0.5,0.25*dz)), failedUnitTest());
        assert(approxEqualNumbers(volume, to!number(cos(radians(30))*0.5*dz/3)), failedUnitTest());

        // tetrahedron properties - complex derivative test
        version(complex_numbers) {
            number volume_dash;
            // Build regular tetrahedron with equal edge lengths (edge length equal to alpha).
            p0_dash = Vector3(zero, zero, zero);
            p1_dash = Vector3(alpha+ih, zero, zero);
            p2_dash = Vector3((alpha+ih)*cos(radians(60)), (alpha+ih)*sin(radians(60)), zero);
            p3_dash = Vector3(0.5*(alpha+ih), 0.5*(alpha+ih)*tan(to!number(radians(30))), (alpha+ih)*(sqrt(6.0)/3.0));
            tetrahedron_properties(p0_dash, p1_dash, p2_dash, p3_dash, centroid_dash, volume_dash);
            double dvol_da = volume_dash.im/h;
            dcentroid_da = Vector3(centroid_dash.x.im/h, centroid_dash.y.im/h, centroid_dash.z.im/h);
            assert(isClose(dvol_da, alpha^^2/(2.0*sqrt(2.0))), failedUnitTest());
            assert(approxEqualVectors(dcentroid_da, Vector3(3.0/8.0+0.25*cos(radians(60)),
                                                            1.0/4.0*sin(radians(60))+1.0/8.0*tan(radians(30)),
                                                            sqrt(6.0)/12.0)),
                   failedUnitTest());
        }

        // Build a wedge with the same equilateral-triangle base.
        p0 = Vector3(0, 0, 0);
        p1 = Vector3(cos(radians(30)), sin(radians(30)), 0.0);
        p2 = Vector3(0.0, 1.0, 0.0);
        Vector3 incz = Vector3(0, 0, -1.0);
        p3 = p0 + incz;
        Vector3 p4 = p1 + incz;
        Vector3 p5 = p2 + incz;
        wedge_properties(p0, p1, p2, p3, p4, p5, true, centroid, volume);
        assert(approxEqualVectors(centroid, Vector3(dx,0.5,-0.5)), failedUnitTest());
        assert(approxEqualNumbers(volume, to!number(cos(radians(30))*0.5*1.0)), failedUnitTest());

        // wedge properties - complex derivative test
        version(complex_numbers) {
            // Build a wedge with an equilateral-triangle base (edge length equal to alpha).
            p0_dash = Vector3(zero, zero, zero);
            p1_dash = Vector3(alpha+ih, zero, zero);
            p2_dash = Vector3((alpha+ih)*cos(radians(60)), (alpha+ih)*sin(radians(60)), zero);
            p3_dash = p0_dash + incz;
            auto p4_dash = p1_dash + incz;
            auto p5_dash = p2_dash + incz;
            wedge_properties(p0_dash, p1_dash, p2_dash, p3_dash, p4_dash, p5_dash, true, centroid_dash, volume_dash);
            dvol_da = volume_dash.im/h;
            dcentroid_da = Vector3(centroid_dash.x.im/h, centroid_dash.y.im/h, centroid_dash.z.im/h);
            assert(isClose(dvol_da, alpha*sin(radians(60))), failedUnitTest());
            assert(approxEqualVectors(dcentroid_da, Vector3((1.0/3.0)*(1.0+cos(radians(60))), (1.0/3.0)*sin(radians(60)), 0.0)), failedUnitTest());
        }

        // Pyramid
        p0 = Vector3(0,0,0); p1 = Vector3(1,0,0);
        p2 = Vector3(1,1,0); p3 = Vector3(0,1,0);
        p4 = Vector3(0.5,0.5,1); // peak
        pyramid_properties(p0, p1, p2, p3, p4, true, centroid, volume);
        assert(approxEqualVectors(centroid, Vector3(0.5,0.5,1.0/4)), failedUnitTest());
        assert(approxEqualNumbers(volume, to!number(1.0/3)), failedUnitTest());

        // pyramid properties - complex derivative test
        version(complex_numbers) {
            p0_dash = Vector3(zero,zero,zero); p1_dash = Vector3(alpha+ih,zero,zero);
            p2_dash = Vector3(alpha+ih,alpha+ih,zero); p3_dash = Vector3(zero,alpha+ih,zero);
            p4_dash = Vector3(0.5,0.5,1); // peak
            pyramid_properties(p0_dash, p1_dash, p2_dash, p3_dash, p4_dash, true, centroid_dash, volume_dash);
            dvol_da = volume_dash.im/h;
            dcentroid_da = Vector3(centroid_dash.x.im/h, centroid_dash.y.im/h, centroid_dash.z.im/h);
            assert(isClose(dvol_da, 2.0*alpha/3.0), failedUnitTest());
            // [TODO] centroid sensitivity
        }

        // Simple cube for the hex cell.
        p0 = Vector3(0,0,0); p1 = Vector3(1,0,0);
        p2 = Vector3(1,1,0); p3 = Vector3(0,1,0);
        p4 = Vector3(0,0,1); p5 = Vector3(1,0,1);
        Vector3 p6 = Vector3(1,1,1); Vector3 p7 = Vector3(0,1,1);
        number iLen, jLen, kLen;
        hex_cell_properties(p0, p1, p2, p3, p4, p5, p6, p7, true, centroid, volume,
                            iLen, jLen, kLen);
        assert(approxEqualVectors(centroid, Vector3(0.5,0.5,0.5)), failedUnitTest());
        assert(approxEqualNumbers(volume, to!number(1.0)), failedUnitTest());

        // Simple cube for the hex cell properties - complex derivative test
        version(complex_numbers) {
            p0_dash = Vector3(0,0,0); p1_dash = Vector3(alpha+ih,zero,zero);
            p2_dash = Vector3(alpha+ih,alpha+ih,zero); p3_dash = Vector3(zero,alpha+ih,zero);
            p4_dash = Vector3(zero,zero,alpha+ih); p5_dash = Vector3(alpha+ih,zero,alpha+ih);
            auto p6_dash = Vector3(alpha+ih,alpha+ih,alpha+ih); auto p7_dash = Vector3(zero,alpha+ih,alpha+ih);
            hex_cell_properties(p0_dash, p1_dash, p2_dash, p3_dash, p4_dash, p5_dash, p6_dash, p7_dash,
                                true, centroid_dash, volume_dash, iLen, jLen, kLen);
            dvol_da = volume_dash.im/h;
            dcentroid_da = Vector3(centroid_dash.x.im/h, centroid_dash.y.im/h, centroid_dash.z.im/h);
            assert(isClose(dvol_da, 3.0*alpha^^2), failedUnitTest());
            assert(approxEqualVectors(dcentroid_da, Vector3(0.5, 0.5, 0.5)), failedUnitTest());
        }

        d = Vector3(0.65, 0.0, 0.0);
        e = Vector3(0.65, -0.1, 0.0);

        p0 = Vector3(0.0, 0.0, 0.0);
        p1 = Vector3(1.0, -0.1, 0.0);
        p2 = Vector3(1.0, 1.0, 0.0);
        p3 = Vector3(0.0, 1.0, 0.0);
        p4 = Vector3(0.0, 0.0, 1.0);
        p5 = Vector3(1.0, -0.1, 1.0);
        p6 = Vector3(1.0, 1.0, 1.0);
        p7 = Vector3(0.0, 1.0, 1.0);
        assert(inside_hexahedron(p0, p1, p2, p3, p4, p5, p6, p7, d), failedUnitTest());
        assert(!inside_hexahedron(p0, p1, p2, p3, p4, p5, p6, p7, e), failedUnitTest());

        auto f = Vector3(0.1, 0.0, 0.5);
        assert(inside_tetrahedron(p0, p1, p2, p4, f), failedUnitTest());
        f = Vector3(0.0, 0.2, 0.5);
        assert(!inside_tetrahedron(p0, p1, p2, p4, f), failedUnitTest());
        f = Vector3(0.0, -0.2, 0.5);
        assert(!inside_tetrahedron(p0, p1, p2, p4, f), failedUnitTest());
        f = Vector3(1.0, 0.0, 1.0);
        assert(!inside_tetrahedron(p0, p1, p2, p4, f), failedUnitTest());
        f = Vector3(0.0, 0.0, -0.5);
        assert(!inside_tetrahedron(p0, p1, p2, p4, f), failedUnitTest());

        return 0;
    }
} // end properties_test
