// bezier.d
// Peter J. 2017-11-29: Split out of the original path module.

module geom.gpath.bezier;

import std.conv;
import std.math;
import ntypes.complex;
import nm.number;

import geom.geometry_exception;
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
    void elevateDegree(int newDegree)
    {
        /*
         * This algorithm is described in Section 2.5 of
         * Rogers (2001), An Introduction to NURBS: with historical perspective."
         *
         * Rogers cites Forrest for this procedure:
         * Forrest, A.R. (1972)
         * Interactive interpolation and approximation by Bezier polynomials.
         * Comp. J. 15:pp.71--79
         */
        auto currentDegree = B.length - 1;

        if (newDegree < currentDegree) {
            throw new GeometryException("Desired Bezier degree elevation is less than current degree.");
        }

        if (currentDegree == newDegree)
            return;

        auto n_elevations = newDegree - currentDegree;
        Vector3[] B_star;
        B_star.length = B.length + n_elevations;
        foreach (j; 0 .. n_elevations) {
            auto n = B.length - 1;
            B_star[0] = B[0];
            // Blend internal points
            foreach (i; 1 .. n+1) {
                double alpha = to!double(i)/(n + 1);
                B_star[i] = alpha*B[i-1] + (1.0 - alpha)*B[i];
            }
            B_star[n+1] = B[n];
            // Get B ready for use, or prepare for next iteration
            B.length = B.length + 1;
            foreach (i; 0 .. n+2) B[i] = B_star[i];
        }
        // Remember to reset control points for derivative
        // and second derivative curves.
        set_deriv_control_points();
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
                Q[i] = (1.0-t)*Q[i] + t*Q[i+1];
            }
        }
        return Q[0];
    }
} // end class Bezier


version(bezier_test) {
    import util.msg_service;
    int main() {
        import geom.gpath.arc;
        auto a = Vector3([2.0, 2.0, 0.0]);
        auto b = Vector3([1.0, 2.0, 1.0]);
        auto c = Vector3([1.0, 2.0, 0.0]);
        auto abc = new Arc(a, b, c);
        auto d = abc(0.5);
        auto adb = new Bezier([a, d, b]);
        auto e = adb(0.5);
        assert(approxEqualVectors(e, Vector3(1.60355, 2, 0.603553)), failedUnitTest());
        auto ab = new Bezier([a, b]);
        assert(approxEqualVectors(ab.dpdt(0.5), Vector3(-1, 0, 1)), failedUnitTest());
        assert(approxEqualVectors(ab.d2pdt2(0.5), Vector3(0)), failedUnitTest());
        auto acb = new Bezier([a, c, b]);
        assert(approxEqualVectors(acb.dpdt(0.5), Vector3(-1, 0, 1)), failedUnitTest());
        assert(approxEqualVectors(acb.d2pdt2(0.5), Vector3(2,0,2)), failedUnitTest());
        //
        auto adbCopy = adb.dup();
        adbCopy.elevateDegree(4);
        assert(approxEqualVectors(adb(0.5), adbCopy(0.5)), failedUnitTest());
        version(complex_numbers) {
            // Try out the complex derivative evaluation.
            double cubic_bezier_analytic_derivative(double t, size_t pt) {
                // compute the analytic derivative of cubic Bezier curve.
                double value;
                if (pt == 0) value = (1-t)^^3;
                else if (pt == 1) value = 3*t*(1-t)^^2;
                else if (pt == 2) value = 3*(1-t)*t^^2;
                else value = t^^3; // assume pt = 3
                return value;
            }
            // define cubic Bezier curve
            auto P0 = Vector3([2.0, 2.0, 0.0]);
            auto P1 = Vector3([1.0, 2.0, 1.0]);
            auto P2 = Vector3([1.0, 2.0, 0.0]);
            auto P3 = Vector3([2.0, 2.0, 1.0]);
            auto myBez = new Bezier([P0, P1, P2, P3]);
            Vector3[] P; // copy of points to be perturbed
            P ~= P0; P ~= P1; P ~= P2; P ~= P3;
            Vector3[] Po; // copy of original points
            Po ~= P0; Po ~= P1; Po ~= P2; Po ~= P3;
            Bezier myNewBez; double dPt_dP_analytic;
            double dPt_dP_complex_x;
            double dPt_dP_complex_y;
            double dPt_dP_complex_z;
            double h = 1.0e-20;
            number ih = complex(0, h); // complex step-size
            double t = 0.5;
            // we will compute the sensitivity of the midpoint of the Bezier curve
            // with respect to x-coord of each Bezier point.
            foreach ( idx; 0..P.length) {
                // compute analytical derivative
                dPt_dP_analytic = cubic_bezier_analytic_derivative(t, idx);
                // compute complex step derivative
                P[idx].x += ih;
                myNewBez = new Bezier([P[0], P[1], P[2], P[3]]); // perturbed Bezier curve
                dPt_dP_complex_x = myNewBez(t).x.im/ih.im;
                dPt_dP_complex_y = myNewBez(t).y.im/ih.im;
                dPt_dP_complex_z = myNewBez(t).z.im/ih.im;
                // compare deriatives
                assert(std.math.isClose(dPt_dP_analytic, dPt_dP_complex_x), failedUnitTest());
                assert(std.math.isClose(0.0, dPt_dP_complex_y), failedUnitTest());
                assert(std.math.isClose(0.0, dPt_dP_complex_z), failedUnitTest());
                // restore point to original position
                P[idx].x = Po[idx].x;
            }
        }
        return 0;
    }
} // end bezier_test
