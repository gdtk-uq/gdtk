/**
 * geom.d  Geometric (vector) primitives for our 3D world.
 *
 * Author: Peter J.
 * Version: 2014-06-16 first cut.
 *          2015-02-18 Attempt to reduce the number of redundant object copies
 *          by declaring some of the "in Vector3" parameters as "ref const(Vector3)".
 *          This has been somewhat successful, however, a few "in" parameters remain
 *          so that the vector arithmetic is a lot cleaner in the code.
 *          Effectively this hides the creation of Vector3 temporaries 
 *          that would otherwise have to appear explicitly in the code.
 */
module geom;

import std.conv;
import std.stdio;
import std.math;

// Nomenclature that we try to use consistently through the 3D code.
// Symbolic names and indices for the cells' faces.
// The names of the faces of the structured-grid blocks will be the same.
enum Face {
    north = 0,
    east = 1,
    south = 2,
    west = 3,
    top = 4,
    bottom = 5
}

string[] face_name = [ "north", "east", "south", "west", "top", "bottom" ];
uint face_index(string name)
{
    switch ( name ) {
    case "north": return Face.north;
    case "east": return Face.east;
    case "south": return Face.south;
    case "west": return Face.west;
    case "top": return Face.top;
    case "bottom": return Face.bottom;
    default:
	throw new Error(text("Invalid face name:", name));
    }
} // end face_index


// Vector3 is the geometric primitive object that we use in all of our
// higher-order objects.

struct Vector3 {
    public double[3] _p;

    @nogc this(in double[] p)
    {
	switch ( p.length ) {
	case 0: _p[0] = _p[1] = _p[2] = 0.0; break;
	case 1: _p[0] = p[0]; _p[1] = _p[2] = 0.0; break;
	case 2: _p[0] = p[0]; _p[1] = p[1]; _p[2] = 0.0; break;
	default: _p[0] = p[0]; _p[1] = p[1]; _p[2] = p[2]; break;
	}
    }

    @nogc this(in double p0, in double p1=0.0, in double p2=0.0)
    {
 	_p[0] = p0;
	_p[1] = p1;
	_p[2] = p2;
    }

    @nogc this(in Vector3 other)
    {
	_p[] = other._p[];
    }

    // Postblit constructor (Alexandrescu Section 7.1.3.4) so that
    // the copy of the struct can become completely independent of 
    // its source.
    this(this)
    {
	_p = _p.dup;
    }

    // For a lot of geometric work, it will be convenient to use
    // x,y,z notation.
    @nogc @property double x() const { return _p[0]; }
    @nogc @property double y() const { return _p[1]; }
    @nogc @property double z() const { return _p[2]; }
    // Note that the following three properties hand out references
    // to the elements, so that we may change their values.
    @nogc @property ref double refx() { return _p[0]; }
    @nogc @property ref double refy() { return _p[1]; }
    @nogc @property ref double refz() { return _p[2]; }

    @property Vector3 dup() const
    {
	return Vector3(this);
    }

    string toString() const
    {
	return "Vector3(" ~ to!string(_p) ~ ")";
    }

    // Some operators, at least those that make sense.
    Vector3 opUnary(string op)()
	if (op == "+")
    {
	Vector3 result;
	result._p[] = this._p[];
	return result;
    }

    Vector3 opUnary(string op)()
	if (op == "-")
    {
	Vector3 result;
	result._p[] = - this._p[];
	return result;
    }

    Vector3 opBinary(string op)(in Vector3 rhs) const
	if (op == "+")
    {
	Vector3 result;
	result._p[] = this._p[] + rhs._p[];
	return result;
    }

    Vector3 opBinary(string op)(in Vector3 rhs) const
	if (op == "-")
    {
	Vector3 result;
	result._p[] = this._p[] - rhs._p[];
	return result;
    }

    Vector3 opBinary(string op)(in double rhs) const
	if (op == "*")
    {
	Vector3 result;
	result._p[] = this._p[] * rhs;
	return result;
    }

    Vector3 opBinaryRight(string op)(in double lhs) const
	if (op == "*")
    {
	Vector3 result;
	result._p[] = this._p[] * lhs;
	return result;
    }

    Vector3 opBinary(string op)(in double rhs) const
	if (op == "/")
    {
	Vector3 result;
	result._p[] = this._p[] / rhs;
	return result;
    }

    // Assignment operators. (Alexandrescu Section 7.1.5.1)
    @nogc ref Vector3 opAssign(ref Vector3 rhs)
    {
	_p[] = rhs._p[];
	return this;
    }

    @nogc ref Vector3 opAssign(Vector3 rhs)
    {
	_p[] = rhs._p[];
	return this;
    }

    // Combined assignment operators do change the original object.
    @nogc ref Vector3 opOpAssign(string op)(in Vector3 rhs)
	if (op == "+")
    {
	this._p[] += rhs._p[];
	return this;
    }

    @nogc ref Vector3 opOpAssign(string op)(in Vector3 rhs)
	if (op == "-")
    {
	this._p[] -= rhs._p[];
	return this;
    }

    @nogc ref Vector3 opOpAssign(string op)(in double rhs)
	if (op == "*")
    {
	this._p[] *= rhs;
	return this;
    }

    @nogc ref Vector3 opOpAssign(string op)(in double rhs)
	if (op == "/")
    {
	this._p[] /= rhs;
	return this;
    }

    // Other vector-specific operations.

    /**
     * Scales the vector to unit magnitude.
     */
    @nogc ref Vector3 normalize()
    {
        double magnitude = sqrt(this.dot(this));
	if ( magnitude > 0.0 ) {
	    this /= magnitude;
	} else {
	    // Clean up, in case dot() underflows.
	    this._p[0] = this._p[1] = this._p[2] = 0.0;
	}
	return this;
    }

    @nogc double dot(ref const(Vector3) other) const
    {
	return this._p[0] * other._p[0] + 
	    this._p[1] * other._p[1] + this._p[2] * other._p[2];
    }

    // Transform functions used to reorient vector values in the CFD codes.

    /**
     * Rotate v from the global xyz coordinate system into the local frame
     * defined by the orthogonal unit vectors n,t1,t2.
     *
     * We assume, without checking, that these vectors do nicely define 
     * such a local system.
     */
    @nogc void transform_to_local_frame(ref const(Vector3) n,
					ref const(Vector3) t1,
					ref const(Vector3) t2)
    {
	double v_x = this.dot(n); // normal component
	double v_y = this.dot(t1); // tangential component 1
	double v_z = this.dot(t2); // tangential component 2
	_p[0] = v_x;
	_p[1] = v_y;
	_p[2] = v_z;
    }

    /**
     * Rotate v back into the global (xyz) coordinate system.
     */
    @nogc void transform_to_global_frame(ref const(Vector3) n,
					 ref const(Vector3) t1,
					 ref const(Vector3) t2)
    {
	double v_x = _p[0]*n._p[0] + _p[1]*t1._p[0] + _p[2]*t2._p[0]; // global-x
	double v_y = _p[0]*n._p[1] + _p[1]*t1._p[1] + _p[2]*t2._p[1]; // global-y
	double v_z = _p[0]*n._p[2] + _p[1]*t1._p[2] + _p[2]*t2._p[2]; // global-z
	_p[0] = v_x;
	_p[1] = v_y;
	_p[2] = v_z;
    }
    // Change of coordinate system; rotation with translation.

    // Transform coordinates from global frame to local (dash) frame.
    // Local frame is defined by unit vectors (n, t1 and t2) at location c.
    @nogc void transform_to_local_frame(ref const(Vector3) n,
					ref const(Vector3) t1,
					ref const(Vector3) t2,
					ref const(Vector3) c)
    {
	_p[0] -= c._p[0]; _p[1] -= c._p[1]; _p[2] -= c._p[2]; // shift to local origin
	double v_x = this.dot(n); // normal component
	double v_y = this.dot(t1); // tangential component 1
	double v_z = this.dot(t2); // tangential component 2
	_p[0] = v_x;
	_p[1] = v_y;
	_p[2] = v_z;
    }

    /**
     * Rotate v back into the global (xyz) coordinate system.
     */
    @nogc void transform_to_global_frame(ref const(Vector3) n,
					 ref const(Vector3) t1,
					 ref const(Vector3) t2,
					 ref const(Vector3) c)
    {
	double v_x = _p[0]*n._p[0] + _p[1]*t1._p[0] + _p[2]*t2._p[0] + c._p[0]; // global-x
	double v_y = _p[0]*n._p[1] + _p[1]*t1._p[1] + _p[2]*t2._p[1] + c._p[1]; // global-y
	double v_z = _p[0]*n._p[2] + _p[1]*t1._p[2] + _p[2]*t2._p[2] + c._p[2]; // global-z
	_p[0] = v_x;
	_p[1] = v_y;
	_p[2] = v_z;
    }
} // end class Vector3


/**
 * Returns the scalar dot product of two vectors.
 */
@nogc
double dot(in Vector3 v1, in Vector3 v2)
{
    double result = 0.0;
    // Maybe we should be careful with underflow and overflow...
    foreach(i; 0 .. 3) result += v1._p[i] * v2._p[i];
    return result;
}

/**
 * Returns magnitude of the vector.
 */
@nogc
double abs(in Vector3 v)
{
    return sqrt(v.dot(v));
}

/**
 * Returns a unit vector in the same direction as v.
 */
Vector3 unit(in Vector3 v)
{
    return Vector3(v).normalize();
}

/**
 * Vector cross product.
 */
@nogc
Vector3 cross(in Vector3 v1, in Vector3 v2)
{
    Vector3 v3;
    v3._p[0] = v1._p[1] * v2._p[2] - v2._p[1] * v1._p[2];
    v3._p[1] = v2._p[0] * v1._p[2] - v1._p[0] * v2._p[2];
    v3._p[2] = v1._p[0] * v2._p[1] - v2._p[0] * v1._p[1];
    return v3;
}

/**
 * Returns true is all of the components are approximately equal.
 */
@nogc
bool approxEqualVectors(in Vector3 v1, in Vector3 v2)
{
    return (approxEqual(v1._p[0], v2._p[0]) && 
	    approxEqual(v1._p[1], v2._p[1]) &&
	    approxEqual(v1._p[2], v2._p[2]));
}

unittest {
    // Check that we have separate data with the correct values.
    Vector3 a = Vector3([1.0, 2.2, 3.0]);
    Vector3 b = Vector3(1.0);
    assert(a.x == 1.0, "a.x fail");
    assert(a.y == 2.2, "a.y fail");
    assert(a.z == 3.0, "a.z fail");
    assert(a.x == b.x, "a.x == b.x fail");
    assert(b.y == 0.0, "b.y fail");
    assert(b.z == 0.0, "b.z fail");
    // Check operators
    b = -a;
    assert(b.x == -a.x && b.y == -a.y && b.z == -a.z, "unary negation");

    b = Vector3(1.0);
    Vector3 c = a + b;
    assert(c.y == a.y+b.y, "vector addition");
    c = a - b;
    assert(c.y == a.y-b.y, "vector subtraction");
    Vector3 d = a.dup;
    a.refy = 99.0;
    assert(a.y == 99.0 && d.y == 2.2, "dup followed by vector change");
    Vector3 d2 = a;
    a.refy = 3.3;
    assert(a.y == 3.3 && d2.y == 99.0, "assignment followed by vector change");

    Vector3 e = a * 2.0;
    Vector3 f = 3.0 * d;
    assert(e.z == 6.0 && f.z == 9.0, "scalar multiply");
    Vector3 g = d / 3.0;
    assert(g.z == 1.0, "scalar division");

    g += f;
    assert(g.z == 10.0, "plus-assign");
    g /= 2.0;
    assert(g.z == 5.0, "divide-assign");

    Vector3 u = unit(g);
    assert(approxEqual(abs(u), 1.0), "unit, dot, abs");

    Vector3 x = Vector3(1.0, 0.0, 0.0);
    Vector3 y = Vector3(0.0, 1.0, 0.0);
    Vector3 z = cross(x,y);
    Vector3 zref = Vector3(0.0,0.0,1.0);
    assert(approxEqualVectors(z, zref), "cross product");

    Vector3 n = unit(Vector3(1.0,1.0,0.0));
    Vector3 t1 = unit(Vector3(-1.0,1.0,0.0));
    Vector3 t2 = cross(n, t1);
    Vector3 h = Vector3(1.0,0.0,1.0);
    Vector3 h_ref = Vector3(h);
    h.transform_to_local_frame(n, t1, t2);
    assert(approxEqualVectors(h, Vector3(sqrt(1.0/2.0), -sqrt(1.0/2.0), 1.0)),
	   "to_local_frame");
    h.transform_to_global_frame(n, t1, t2);
    assert(approxEqualVectors(h, h_ref), "to_global_frame");
}


//------------------------------------------------------------------------
// Geometry functions projection and mapping.

double radians(double degrees) { return PI*degrees/180.0; }

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
    double D = dot(a, N);

    double numer = D - dot(q, N);
    double denom = dot(qr, N);

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
ref Vector3 map_neutral_plane_to_cylinder(ref Vector3 p, double H)
{
    // The axis of the hypothetical cylinder coincides with the x-axis thus
    // H is also the distance of the neutral plane above the x-axis.
    // For Hannes Wojciak and Paul Petrie-Repar's turbomachinery grids.
    if ( H > 0.0 ) {
	double theta = p.y / H;
	double old_z = p.z;
	p.refz = old_z * cos(theta);
	p.refy = old_z * sin(theta);
	// x remains the same
    }
    return p;
}

int inside_triangle(ref const(Vector3) p, ref const(Vector3) a,
		    ref const(Vector3) b, ref const(Vector3) c)
{
    Vector3 n = unit(cross(a-c, b-c)); // normal to plane of triangle
    double area1 = 0.5 * dot(cross(p-a, p-b), n); // projected area of triangle pab
    double area2 = 0.5 * dot(cross(p-b, p-c), n);
    double area3 = 0.5 * dot(cross(p-c, p-a), n);
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

unittest {
    Vector3 a = Vector3(1.0, 0.0, 0.0); // plane through a,b,c
    Vector3 b = Vector3(1.0, 1.0, 0.0);
    Vector3 c = Vector3(0.5, 0.0, 0.0);
    Vector3 qr = Vector3(3.0, 3.0, -3.0); // direction
    Vector3 q = Vector3(0.0, 0.0, 1.0); // start point
    int flag =  project_onto_plane(q, qr, a, b, c);
    assert(approxEqualVectors(q, Vector3(1.0,1.0,0.0)), "project_onto_plane");
    Vector3 myp = Vector3(1.0, 1.0, 1.0);
    map_neutral_plane_to_cylinder(myp, 1.0);
    assert(approxEqualVectors(myp, Vector3(1.0, sin(1.0), cos(1.0))), "cylinder map");
    Vector3 d = Vector3(0.65, 0.0, 0.0);
    assert(inside_triangle(d, a, b, c) > 0, "inside triangle");
    Vector3 e = Vector3(0.65, -0.1, 0.0);
    assert(!inside_triangle(e, a, b, c), "outside triangle");
}


//------------------------------------------------------------------------
// Geometry functions projection and mapping.

/** Quadrilateral properties of centroid, associated unit normals and area.
 *   p3-----p2
 *   |      |
 *   |      |
 *   p0-----p1
 * Resultant normal vector is up, toward you.
 * Assume that all points are in the one plane.
 */
void quad_properties(ref const(Vector3) p0, ref const(Vector3) p1,
		     ref const(Vector3) p2, ref const(Vector3) p3,
		     ref Vector3 centroid,
		     ref Vector3 n, ref Vector3 t1, ref Vector3 t2,
		     ref double area,
		     double tol=1.0e-12, double area_tol=1.0e-20)
{
    centroid = 0.25 * (p0 + p1 + p2 + p3);
    // Compute areas via the cross products.
    Vector3 vector_area = 0.25 * cross(p1-p0+p2-p3, p3-p0+p2-p1);
    // unit-normal and area
    area = abs(vector_area);
    if ( area > area_tol ) {
	n = unit(vector_area);
	// Tangent unit-vectors: 
	// t1 is parallel to side01 and side32, 
	// t2 is normal to n and t1
	t1 = unit((p1-p0)+(p2-p3)); // Works even if one edge has zero length.
	t2 = unit(cross(n, t1)); // Calling unit() to tighten up the magnitude.
    } else {
	// We have nothing meaningful to put into the unit vectors.
	throw new Exception("Effectively zero area quadrilateral.");
    }
    if ( fabs(abs(n)-1.0) > tol || fabs(abs(t1)-1.0) > tol || fabs(abs(t2)-1.0) > tol ) {
	throw new Exception("Failed to produce unit vectors properly.");
    }
} // end quad_properties()

double tetrahedron_volume(ref const(Vector3) p0, ref const(Vector3) p1,
			  ref const(Vector3) p2, ref const(Vector3) p3)
{
    return dot(p3-p0, cross(p1-p0, p2-p0)) / 6.0;
} // end tetrahedron_volume()

void tetrahedron_properties(ref const(Vector3) p0, ref const(Vector3) p1,
			    ref const(Vector3) p2, ref const(Vector3) p3,
			    ref Vector3 centroid, ref double volume)
{
    volume = tetrahedron_volume(p0, p1, p2, p3);
    centroid = 0.25 * (p0 + p1 + p2 + p3);
} // end tetrahedron_properties()

// Because of the way we lose precision when reading and writing files,
// it may be that the vertices are not quite in their ideal position.
// We need a couple of finite, but small, tolerances to deal with 
// collapsed volumes.
double smallButSignificantVolume = 1.0e-12;
double verySmallVolume = 1.0e-20;

void wedge_properties(in Vector3 p0, in Vector3 p1, in Vector3 p2, 
		     const Vector3 p3, in Vector3 p4, in Vector3 p5,
		     ref Vector3 centroid, ref double volume)
{
    double v1, v2, v3;
    Vector3 c1, c2, c3;
    tetrahedron_properties(p0, p4, p5, p3, c1, v1);
    tetrahedron_properties(p0, p5, p4, p1, c2, v2);
    tetrahedron_properties(p0, p1, p2, p5, c3, v3);
    volume = v1 + v2 + v3;
    if ( (volume < 0.0 && fabs(volume) < smallButSignificantVolume) ||
	 (volume >= 0.0 && volume < verySmallVolume) ) {
	// We assume that we have a collapsed wedge; no real problem.
	volume = 0.0;
	// equally-weighted tetrahedral centroids.
	centroid = (c1 + c2 + c3) / 3.0;
	return;
    }
    if ( volume < 0.0 ) {
	// Something has gone wrong with our wedge geometry.
	centroid = (c1 + c2 + c3) / 3.0;
        throw new Exception("significant negative volume.");
    }
    // Weight the tetrahedral centroids with their volumes.
    centroid = (c1*v1 + c2*v2 + c3*v3) / volume;
} // end wedge_properties()

double tetragonal_dipyramid_volume(in Vector3 p0, in Vector3 p1, 
				   in Vector3 p2, in Vector3 p3, 
				   in Vector3 pb, in Vector3 pc)
// J. Grandy (1997) Efficient Computation of Volume of Hexahedral Cells UCRL-ID-128886.
// Base of each dipyramid is specified clockwise from the outside.
// pc is apex
// pb is barycentre of base quad.
// base quad p0->p1->p2->p3->p0 counterclockwise when looking from pc toward base.
{
    double volume = dot(pc-pb, cross(p1-p0+p2-p3, p3-p0+p2-p1)) / 12.0;
    if ( (volume < 0.0 && fabs(volume) < smallButSignificantVolume) ||
	 (volume >= 0.0 && volume < verySmallVolume) ) {
	// We assume that we have a collapsed pyramid (all points coplanar);
	// no real problem here but it may be a problem for the client code.
	// That code should test the value of volume, on return.
	volume = 0.0;
    }
    if ( volume < 0.0 ) {
	// Something has gone wrong with our geometry.
        throw new Exception("significant negative volume.");
    }
    return volume;
} // end tetragonal_dipyramid_volume()

void hex_cell_properties(in Vector3 p0, in Vector3 p1, in Vector3 p2, in Vector3 p3,
			 in Vector3 p4, in Vector3 p5, in Vector3 p6, in Vector3 p7,
			 ref Vector3 centroid, ref double volume,
			 ref double iLen, ref double jLen, ref double kLen)
{
    // PJ 10-Sep-2012
    // When computing the volume of Rolf's thin, warped cells, we have to do 
    // something better than splitting our cell into six tetrahedra.
    centroid = 0.125 * (p0+p1+p2+p3+p4+p5+p6+p7);
    // Mid-points of faces.
    Vector3 pmN = 0.25*(p3+p2+p6+p7);
    Vector3 pmE = 0.25*(p1+p2+p6+p5);
    Vector3 pmS = 0.25*(p0+p1+p5+p4);
    Vector3 pmW = 0.25*(p0+p3+p7+p4);
    Vector3 pmT = 0.25*(p4+p5+p6+p7);
    Vector3 pmB = 0.25*(p0+p1+p2+p3);
    // Lengths between mid-points of faces.
    // Note that we are assuming that the hexahedron is not very skewed
    // when we later use these values as the widths of the hex cell.
    iLen = abs(pmE - pmW);
    jLen = abs(pmN - pmS);
    kLen = abs(pmT - pmB);
    // writeln("Single hexahedron divided into six tetragonal dipyramids.");
    // J. Grandy (1997) Efficient Computation of Volume of Hexahedral Cells UCRL-ID-128886.
    // Base of each dipyramid is specified clockwise from the outside.
    volume = 0.0;
    volume += tetragonal_dipyramid_volume(p6, p7, p3, p2, pmN, centroid); // North
    volume += tetragonal_dipyramid_volume(p5, p6, p2, p1, pmE, centroid); // East
    volume += tetragonal_dipyramid_volume(p4, p5, p1, p0, pmS, centroid); // South
    volume += tetragonal_dipyramid_volume(p7, p4, p0, p3, pmW, centroid); // West
    volume += tetragonal_dipyramid_volume(p7, p6, p5, p4, pmT, centroid); // Top
    volume += tetragonal_dipyramid_volume(p0, p1, p2, p3, pmB, centroid); // Bottom
    if ( (volume < 0.0 && fabs(volume) < smallButSignificantVolume) ||
	 (volume >= 0.0 && volume < verySmallVolume) ) {
	// We assume that we have a collapsed hex cell;
	// no real problem here but it may be a problem for the client code.
	// That code should test the value of volume, on return.
	volume = 0.0;
    }
    if ( volume < 0.0 ) {
	// Something has gone wrong with our geometry.
        throw new Exception("significant negative volume.");
    }
    return; 
} // end hex_cell_properties()


unittest {
    Vector3 p0 = Vector3(0.0, 0.0, 1.0);
    Vector3 p1 = Vector3(1.0, 0.0, 1.0);
    Vector3 p2 = Vector3(1.0, 1.0, 1.0);
    Vector3 p3 = Vector3(0.0, 1.0, 1.0);
    Vector3 centroid, n, t1, t2;
    double area;
    quad_properties(p0, p1, p2, p3, centroid, n, t1, t2, area);
    assert(approxEqual(area, 1.0), "quad_properties area");
    assert(approxEqualVectors(centroid, Vector3(0.5,0.5,1.0)), "quad_properties centroid");
    assert(approxEqualVectors(n, Vector3(0.0,0.0,1.0)), "quad_properties normal");
    assert(approxEqualVectors(t1, Vector3(1.0,0.0,0.0)), "quad_properties t1");
    assert(approxEqualVectors(t2, Vector3(0.0,1.0,0.0)), "quad_properties t2");

    // Build tetrahedron with equilateral triangle (side 1.0) on xy plane.
    p0 = Vector3(0, 0, 0);
    p1 = Vector3(cos(radians(30)), sin(radians(30)), 0.0);
    p2 = Vector3(0.0, 1.0, 0.0);
    double dx = 0.5 * tan(radians(30));
    double dL = cos(radians(30));
    double dz = sqrt(dL*dL - dx*dx);
    p3 = Vector3(dx, 0.5, dz);
    double volume;
    tetrahedron_properties(p0, p1, p2, p3, centroid, volume);
    assert(approxEqualVectors(centroid, Vector3(dx,0.5,0.25*dz)), "tetrahedron centroid");
    assert(approxEqual(volume, cos(radians(30))*0.5*dz/3), "tetrahedron volume");

    // Build a wedge with the same equilateral-triangle base.
    p3 = p0 + Vector3(0, 0, 1.0);
    Vector3 p4 = p1 + Vector3(0, 0, 1.0);
    Vector3 p5 = p2 + Vector3(0, 0, 1.0);
    wedge_properties(p0, p1, p2, p3, p4, p5, centroid, volume);
    assert(approxEqualVectors(centroid, Vector3(dx,0.5,0.5)), "wedge centroid");
    assert(approxEqual(volume, cos(radians(30))*0.5*1.0), "wedge volume");

    // Simple cube for the hex cell.
    p0 = Vector3(0,0,0); p1 = Vector3(1,0,0);
    p2 = Vector3(1,1,0); p3 = Vector3(0,1,0);
    p4 = Vector3(0,0,1); p5 = Vector3(1,0,1);
    Vector3 p6 = Vector3(1,1,1); Vector3 p7 = Vector3(0,1,1);
    double iLen, jLen, kLen;
    hex_cell_properties(p0, p1, p2, p3, p4, p5, p6, p7, centroid, volume,
			iLen, jLen, kLen);
    assert(approxEqualVectors(centroid, Vector3(0.5,0.5,0.5)), "hex centroid");
    assert(approxEqual(volume, 1.0), "hex volume");
}

