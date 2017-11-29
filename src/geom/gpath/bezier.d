// bezier.d
// Peter J. 2017-11-29: Split out of the original path module.

module geom.gpath.bezier;

import std.conv;
import std.math;

import geom.elements;
import geom.gpath.path;


class Bezier : Path {
public:
    Vector3[] B; // collection of control points
    this(in Vector3[] B)
    {
	if (B.length == 0) {
	    throw new Error(text("Bezier() No control points present."));
	}
	if (B.length == 1) {
	    throw new Error(text("Bezier() Only one control point, not enough for a curve."));
	}
	this.B = B.dup();
	set_deriv_control_points();
    }
    this(ref const(Bezier) other)
    {
	this(other.B);
    }
    override Bezier dup() const
    {
	return new Bezier(B);
    }
    override Vector3 opCall(double t) const 
    {
	// Evaluate B(t)
	return deCasteljau(B, t);
    } // end opCall()
    override Vector3 dpdt(double t) const
    {
	return deCasteljau(C, t);
    }
    override Vector3 d2pdt2(double t) const
    {
	return deCasteljau(D, t);
    }
    override string toString() const
    {
	return "Bezier(B=" ~ to!string(B) ~ ")";
    }
    override string classString() const
    {
	return "Bezier";
    }

protected:
    Vector3[] C; // derivative curve
    Vector3[] D; // second derivative curve
    void set_deriv_control_points()
    {
	size_t n = B.length - 1;
	if ( n == 0) {
	    // shouldn't reach here due to check in constructor
	    throw new Error(text("Bezier() Curve is a point, derivative not defined"));
	}
	C.length = n;
	foreach (i; 0 .. n){
	    C[i] = n*(B[i+1] - B[i]);
	}
	if ( n == 1 ) {
	    D = [Vector3(0)];
	    return;
	}
	D.length = n - 1;
	foreach (i; 0 .. n-1){
	    D[i] = (n-1)*(C[i+1] - C[i]);
	}
    }
    Vector3 deCasteljau(ref const(Vector3[]) B, double t) const
    {
	if ( B.length == 1 ) return B[0];
	size_t n_order = B.length - 1;
	// Apply de Casteljau's algorithm. 
	Vector3[] Q = B.dup(); // work array will be overwritten
	foreach (k; 0 .. n_order) {
	    foreach (i; 0 .. n_order-k) {
		Q[i] = (1.0 - t) * Q[i] + t * Q[i+1];
	    }
	}
	return Q[0];
    }
} // end class Bezier


unittest {
    import geom.gpath.arc;
    auto a = Vector3([2.0, 2.0, 0.0]);
    auto b = Vector3([1.0, 2.0, 1.0]);
    auto c = Vector3([1.0, 2.0, 0.0]);
    auto abc = new Arc(a, b, c);
    auto d = abc(0.5);
    auto adb = new Bezier([a, d, b]);
    auto e = adb(0.5);
    assert(approxEqualVectors(e, Vector3(1.60355, 2, 0.603553)), "Bezier");
    auto ab = new Bezier([a, b]);
    assert(approxEqualVectors(ab.dpdt(0.5), Vector3(-1, 0, 1)), "Bezier");
    assert(approxEqualVectors(ab.d2pdt2(0.5), Vector3(0)), "Bezier");
    auto acb = new Bezier([a, c, b]);
    assert(approxEqualVectors(acb.dpdt(0.5), Vector3(-1, 0, 1)), "Bezier");
    assert(approxEqualVectors(acb.d2pdt2(0.5), Vector3(2,0,2)), "Bezier");
}
