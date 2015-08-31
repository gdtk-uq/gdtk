/** surface.d
 * Geometry-building elements for our 3D world -- two-parameter surfaces.
 *
 * Author: Peter J and Rowan G.
 * Version: 2015-02-19 first code
 */

module surface;

import std.math;
import std.stdio;
import std.conv;
import geom;
import gpath;

// Nomenclature for the parametric distances, bounding paths and corners.
//
//         north
//     p01-------p11 s=1
//      |         |
// west |         | east
//      |         |
//     p00-------p10 s=0
//         south
//     r=0       r=1
//
// We'll try to use this notation consistently in the classes below.

class ParametricSurface {
public:
    abstract Vector3 opCall(double r, double s) const;
    abstract ParametricSurface dup() const;
    abstract override string toString() const;
} // end class ParametricSurface


class CoonsPatch : ParametricSurface {
public:
    Path north, east, south, west; // bounding paths
    Vector3 p00, p10, p11, p01;    // corners

    this(in Vector3 p00, in Vector3 p10, in Vector3 p11, in Vector3 p01)
    {
	north = new Line(p01, p11);
	east = new Line(p10, p11);
	south = new Line(p00, p10);
	west = new Line(p00, p01);
	this(south, north, west, east);
    }

    this(in Path south, in Path north, in Path west, in Path east)
    // The particular order for the boundary surfaces goes way back
    // to the original grid generation paper, so it doesn't match the
    // default order of NESW in the rest of the flow code.
    {
	this.north = north.dup();
	this.east = east.dup();
	this.south = south.dup();
	this.west = west.dup();
	p00 = south(0.0);
	p10 = south(1.0);
	p01 = north(0.0);
	p11 = north(1.0);
	// Check alternate evaluation of corners for consistency.
	Vector3 p00_alt = west(0.0);
	Vector3 p10_alt = east(0.0);
	Vector3 p01_alt = west(1.0);
	Vector3 p11_alt = east(1.0);
	if (!approxEqualVectors(p00, p00_alt)) {
	    throw new Error(text("CoonsPatch open corner p00= ", p00,
				 " p00_alt= ", p00_alt));
	}
	if (!approxEqualVectors(p10, p10_alt)) {
	    throw new Error(text("CoonsPatch open corner p10= ", p10,
				 " p10_alt= ", p10_alt));
	}
	if (!approxEqualVectors(p11, p11_alt)) {
	    throw new Error(text("CoonsPatch open corner p11= ", p11,
				 " p11_alt= ", p11_alt));
	}
	if (!approxEqualVectors(p01, p01_alt)) {
	    throw new Error(text("CoonsPatch open corner p01= ", p01,
				 " p01_alt= ", p01_alt));
	}
    }

    this(ref const(CoonsPatch) other)
    {
	this.north = other.north.dup();
	this.east = other.east.dup();
	this.south = other.south.dup();
	this.west = other.west.dup();
	p00 = other.p00;
	p10 = other.p10;
	p01 = other.p01;
	p11 = other.p11;
    }

    override CoonsPatch dup() const
    {
	return new CoonsPatch(this.south, this.north, this.west, this.east);
    }

    override Vector3 opCall(double r, double s) const 
    {
	Vector3 south_r = south(r); 
	Vector3 north_r = north(r);
	Vector3 west_s = west(s); 
	Vector3 east_s = east(s);
	Vector3 p = (1.0-s)*south_r + s*north_r + (1.0-r)*west_s + r*east_s - 
	    ( (1.0-r)*(1.0-s)*p00 + (1.0-r)*s*p01 + r*(1.0-s)*p10 + r*s*p11 );
	return p;
    }

    override string toString() const
    {
	return "CoonsPatch(south=" ~ to!string(south) ~
	    ", north=" ~ to!string(north) ~
	    ", west=" ~ to!string(west) ~
	    ", east=" ~ to!string(east) ~
	    ")";
    }
} // end class CoonsPatch


unittest {
    auto p00 = Vector3([0.0, 0.1, 3.0]);
    auto p10 = Vector3(1.0, 0.1, 3.0);
    auto p11 = Vector3(1.0, 1.1, 3.0);
    auto p01 = Vector3(0.0, 1.1, 3.0);
    auto my_patch = new CoonsPatch(p00, p10, p11, p01);
    auto c = my_patch(0.5, 0.5);
    assert(approxEqualVectors(c, Vector3(0.5, 0.6, 3.0)), "CoonsPatch mid");
    c = my_patch(0.1, 0.1);
    assert(approxEqualVectors(c, Vector3(0.1, 0.2, 3.0)), "CoonsPatch lower-left");
}


class AOPatch : ParametricSurface {
    // Another surface defined as a blend of 4 bounding paths.
    //
    // The topology of the parametric surface is the same as that of the CoonsPatch.
    // The difference, however, is in the evaluation of points on the surface.
    // Here, points are interpolated within a background mesh that has been fitted
    // with the Area-Orthogonality (AO) elliptic grid generator 
    // described by Patrick M. Knupp "A Robust Elliptic Grid Generator"
    // J. Computational Physics Vol.100 pp409-418 (1992)
public:
    Path north, east, south, west; // bounding paths
    Vector3 p00, p10, p11, p01;    // corners

private:
    CoonsPatch tfi_surface;
    int _nx, _ny;
    Vector3[][] _bgmesh;

public:
    this(in Vector3 p00, in Vector3 p10, in Vector3 p11, in Vector3 p01,
	 int nx=10, int ny=10)
    {
	north = new Line(p01, p11);
	east = new Line(p10, p11);
	south = new Line(p00, p10);
	west = new Line(p00, p01);
	this(south, north, west, east, nx, ny);
    }

    this(in Path south, in Path north, in Path west, in Path east,
	 int nx=10, int ny=10)
    {
	this.north = north.dup();
	this.east = east.dup();
	this.south = south.dup();
	this.west = west.dup();
	_nx = nx; _ny = ny;
	// Set up internal representation.
	p00 = south(0.0);
	p10 = south(1.0);
	p01 = north(0.0);
	p11 = north(1.0);
	// Check alternate evaluation of corners for consistency.
	Vector3 p00_alt = west(0.0);
	Vector3 p10_alt = east(0.0);
	Vector3 p01_alt = west(1.0);
	Vector3 p11_alt = east(1.0);
	if (!approxEqualVectors(p00, p00_alt)) {
	    throw new Error(text("AOPatch open corner p00= ", p00,
				 " p00_alt= ", p00_alt));
	}
	if (!approxEqualVectors(p10, p10_alt)) {
	    throw new Error(text("AOPatch open corner p10= ", p10,
				 " p10_alt= ", p10_alt));
	}
	if (!approxEqualVectors(p11, p11_alt)) {
	    throw new Error(text("AOPatch open corner p11= ", p11,
				 " p11_alt= ", p11_alt));
	}
	if (!approxEqualVectors(p01, p01_alt)) {
	    throw new Error(text("AOPatch open corner p01= ", p01,
				 " p01_alt= ", p01_alt));
	}
	// The TFI grid is a fall-back way of evaluating a point on
	// the surface and will be used when the parameter values 
	// place us very close to a boundary.
	// The coarse background mesh tends to cut across 
	// curved boundaries and so should not be used as our
	// guide near those boundaries.
	tfi_surface = new CoonsPatch(south, north, west, east);
	compute_background_mesh();
    }

    this(ref const(AOPatch) other)
    {
	this.north = other.north.dup();
	this.east = other.east.dup();
	this.south = other.south.dup();
	this.west = other.west.dup();
	p00 = other.p00;
	p10 = other.p10;
	p01 = other.p01;
	p11 = other.p11;
	tfi_surface = other.tfi_surface.dup();
	compute_background_mesh();
    }

    override AOPatch dup() const
    {
	return new AOPatch(this.south, this.north, this.west, this.east,
			   this._nx, this._ny);
    }

    override Vector3 opCall(double r, double s) const 
    {
	// Use TFI close to the boundaries.  
	// The background mesh is usually pretty coarse.
	double tol = 1.0e-4;
	if ( r < tol || s < tol || r > 1.0-tol || s > 1.0-tol ) {
	    Vector3 p = tfi_surface(r, s);
	    return p;
	}
	// Interpolate within the AO background mesh.
	// This involves finding the relevant coarse cell and
	// using a bilinear interpolation within that cell.
	double dr = 1.0 / _nx;
	double ds = 1.0 / _ny;
	int ix_coarse = to!int(r / dr);
	int iy_coarse = to!int(s / ds);
	// Parametric coordinate within the coarse cell.
	double local_r = (r - dr * ix_coarse) / dr;
	double local_s = (s - ds * iy_coarse) / ds;
	// BiLinear interpolation for each component.
	Vector3 p = (1.0 - local_r) * (1.0 - local_s) * _bgmesh[ix_coarse][iy_coarse] +
	    (1.0 - local_r) * local_s * _bgmesh[ix_coarse][iy_coarse + 1] +
	    local_r * (1.0 - local_s) * _bgmesh[ix_coarse + 1][iy_coarse] +
	    local_r * local_s * _bgmesh[ix_coarse + 1][iy_coarse + 1];
 	return p;
    } // end opCall()

    override string toString() const
    {
	return "AOPatch(south=" ~ to!string(south) ~
	    ", north=" ~ to!string(north) ~
	    ", west=" ~ to!string(west) ~
	    ", east=" ~ to!string(east) ~
	    ", nx=" ~ to!string(_nx) ~ ", ny=" ~ to!string(_ny) ~
	    ")";
    } // end toString()

private:
    void compute_background_mesh()
    {
	// Inflate array for holding nicely orthogonal mesh.
	_bgmesh.length = _nx + 1;
	foreach (ix; 0 .. _nx+1) _bgmesh[ix].length = _ny + 1;
	// Initial positions of the mesh points are just TFI locations.
	double dXi = 1.0 / _nx;
	double dEta = 1.0 / _ny;
	foreach (ix; 0 .. _nx+1) {
	    foreach (iy; 0 .. _ny+1) {
		double r = dXi * ix;
		double s = dEta * iy;
		_bgmesh[ix][iy] = tfi_surface(r, s);
	    }
	}
	// Now, adjust the mesh point locations.
	// Note that the current formulation is for the (x,y)-plane only.
	// z-components remain as the TFI value.
	double dXi2 = dXi * dXi;
	double dEta2 = dEta * dEta;
	double x_tol = 1.0e-6;
	double y_tol = 1.0e-6;
	double largest_x_move;
	double largest_y_move;
	int max_count = 500;
	int count;
	for (count = 1; count <= max_count; ++count) {
	    largest_x_move = 0.0;
	    largest_y_move = 0.0;
	    // Adjust the internal points only.
	    foreach (ix; 1 .. _nx) {
		foreach (iy; 1 .. _ny) {
		    // Save the old position.
		    double x_old = _bgmesh[ix][iy].x;
		    double y_old = _bgmesh[ix][iy].y;
		    // Calculate the partial derivatives.
		    double dxdXi = (_bgmesh[ix+1][iy].x - _bgmesh[ix-1][iy].x) / (2.0 * dXi);
		    double dxdEta = (_bgmesh[ix][iy+1].x - _bgmesh[ix][iy-1].x) / (2.0 * dEta);
		    double d2xdXidEta = ((_bgmesh[ix+1][iy+1].x - _bgmesh[ix-1][iy+1].x) -
					 (_bgmesh[ix+1][iy-1].x - _bgmesh[ix-1][iy-1].x)) / 
			(4.0 * dXi * dEta);
		    double dydXi = (_bgmesh[ix+1][iy].y - _bgmesh[ix-1][iy].y) / (2.0 * dXi);
		    double dydEta = (_bgmesh[ix][iy+1].y - _bgmesh[ix][iy-1].y) / (2.0 * dEta);
		    double d2ydXidEta = ((_bgmesh[ix+1][iy+1].y - _bgmesh[ix-1][iy+1].y) -
					 (_bgmesh[ix+1][iy-1].y - _bgmesh[ix-1][iy-1].y)) / 
			(4.0 * dXi * dEta);
		    // Calculate intermediate quantities
		    double B = dxdXi * dydEta + dxdEta * dydXi;
		    double Ax = (4.0 * dxdXi * dxdEta * d2xdXidEta +
				 2.0 * B * d2ydXidEta) * (dXi2 * dEta2);
		    double Ay = (4.0 * dydXi * dydEta * d2ydXidEta +
				 2.0 * B * d2xdXidEta) * (dXi2 * dEta2);
		    double g11 = dxdXi * dxdXi + dydXi * dydXi;
		    double g22 = dxdEta * dxdEta + dydEta * dydEta;
		    // Update the node position.
		    double numer = Ax + g22 * dEta2 * (_bgmesh[ix+1][iy].x + _bgmesh[ix-1][iy].x) +
			g11 * dXi2 * (_bgmesh[ix][iy+1].x + _bgmesh[ix][iy-1].x);
		    double denom = 2.0 * (g22 * dEta2 + g11 * dXi2);
		    _bgmesh[ix][iy].refx = numer / denom;

		    numer = Ay + g22 * dEta2 * (_bgmesh[ix+1][iy].y + _bgmesh[ix-1][iy].y) +
			g11 * dXi2 * (_bgmesh[ix][iy+1].y + _bgmesh[ix][iy-1].y);
		    denom = 2.0 * (g22 * dEta2 + g11 * dXi2);
		    _bgmesh[ix][iy].refy = numer / denom;
		    double dx = fabs(_bgmesh[ix][iy].x - x_old);
		    double dy = fabs(_bgmesh[ix][iy].y - y_old);
		    if ( dx > largest_x_move ) largest_x_move = dx;
		    if ( dy > largest_y_move ) largest_y_move = dy;
		    // writefln("Iteration %d node[%d, %d] x,y(%f, %f)\n",
		    //          count, ix, iy, _bgmesh[ix][iy].x, _bgmesh[ix][iy].y);
		} //iy loop
	    } // ix loop
	    // Check for convergence
	    if (largest_x_move <= x_tol && largest_y_move <= y_tol) break;
	} // count loop
	if (count > max_count) {
	    writeln("Warning: AO iteration did not converge.");
	}
	return;
    } // end compute_background_mesh()
} // end class AOPatch

unittest {
    auto p00 = Vector3([0.0, 0.1, 3.0]);
    auto p10 = Vector3(1.0, 0.4, 3.0);
    auto p11 = Vector3(1.0, 1.1, 3.0);
    auto p01 = Vector3(0.0, 1.1, 3.0);
    auto my_patch = new AOPatch(p00, p10, p11, p01);
    auto c = my_patch(0.5, 0.5);
    assert(approxEqualVectors(c, Vector3(0.443775, 0.652291, 3.0)),
	   "AOPatch mid");
    c = my_patch(0.1, 0.1);
    assert(approxEqualVectors(c, Vector3(0.0892476, 0.215529, 3.0)),
	   "AOPatch lower-left");
}


class ChannelPatch : ParametricSurface {
public:
    Path cA; // south edge
    Path cB; // north edge
    bool ruled;
    bool pure2D;
    Vector3 p00, p10, p11, p01;
    
    this(const Path cA, const Path cB, bool ruled=false, bool pure2D=false)
    {
	this.cA = cA.dup();
	this.cB = cB.dup();
	this.ruled = ruled;
	this.pure2D = pure2D;
	p00 = cA(0.0);
	p10 = cA(1.0);
	p01 = cB(0.0);
	p11 = cB(1.0);
    }
    
    this(ref const(ChannelPatch) other)
    {
	cA = other.cA.dup();
	cB = other.cB.dup();
	ruled = other.ruled;
	pure2D = other.pure2D;
	p00 = other.p00;
	p10 = other.p10;
	p01 = other.p01;
	p11 = other.p11;
    }

    override ChannelPatch dup() const
    {
	return new ChannelPatch(this.cA, this.cB, this.ruled, this.pure2D);
    }

    override Vector3 opCall(double r, double s) const 
    {
	auto bridge_path = make_bridging_path(r);
	Vector3 p = bridge_path(s);
	if ( pure2D ) p.refz = 0.0;
	return p;
    }

    override string toString() const
    {
	return "ChannelPatch(cA=" ~ to!string(cA) ~ ", cB=" ~ to!string(cB) ~
	    ", ruled=" ~ to!string(ruled) ~ ", pure2D=" ~ to!string(pure2D) ~ ")";
    }
    
    Path make_bridging_path(double r) const
    {
	Vector3 cAr = cA(r); 
	Vector3 cBr = cB(r);
	if ( pure2D ) { cAr.refz = 0.0; cBr.refz = 0.0; }
	if ( ruled ) {
	    // Bridge with a straight line for a ruled surface.
	    return new Line(cAr, cBr);
	} else {
	    // Bridge with a Bezier3 path that is normal to both defining curves.
	    Vector3 pBminuspA = cBr - cAr;
	    double L = abs(pBminuspA);
	    Vector3 dcArdt = cA.dpdt(r);
	    Vector3 dcBrdt = cB.dpdt(r);
	    // Out-of-plane vectors
	    Vector3 oopvA = cross(dcArdt, pBminuspA);
	    Vector3 oopvB = cross(dcBrdt, pBminuspA);
	    // Inward-facing normal vectors on the surface.
	    Vector3 nA = cross(oopvA, dcArdt).normalize();
	    Vector3 nB = cross(dcBrdt, oopvB).normalize();
	    // Intermediate control points for the cubic Bezier.
	    Vector3 p1 = cAr + L/3.0*nA;
	    Vector3 p2 = cBr + L/3.0*nB;
	    if ( pure2D ) { p1.refz = 0.0; p2.refz = 0.0; }
	    return new Bezier([cAr, p1, p2, cBr]);
	}
    } // end make_bridging_path()
} // end class ChannelPatch


class SubRangedSurface : ParametricSurface {
public:
    ParametricSurface underlying_surface;
    double r0; // to subrange r, when evaluating a point on the surface
    double r1;
    double s0;
    double s1;

    this(const ParametricSurface psurf,
	 double r0=0.0, double r1=1.0, double s0=0.0, double s1=1.0)
    {
	underlying_surface = psurf.dup();
	this.r0 = r0;
	this.r1 = r1;
	this.s0 = s0;
	this.s1 = s1;
    }
    override Vector3 opCall(double r, double s) const
    {
	r = r0 + (r1-r0)*r; // subrange the parameter
	s = s0 + (s1-s0)*s;
	return underlying_surface(r,s);
    }
    override ParametricSurface dup() const
    {
	return new SubRangedSurface(underlying_surface, r0, r1, s0, s1);
    }
    override string toString() const
    {
	return "SubRangedSurface(underlying_surface=" ~
	    to!string(underlying_surface) ~
	    ", r0=" ~ to!string(r0) ~ ", r1=" ~ to!string(r1) ~
	    ", s0=" ~ to!string(s0) ~ ", s1=" ~ to!string(s1) ~
	    ")";
    }
} // end class SubRangedSurface
