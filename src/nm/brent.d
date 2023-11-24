/**
 * brent.d
 *
 * Solve a nonlinear equation f(x)=0 using the method of Brent. 
 * ref:pg_454_numerical_recipes_3rd[PRESS, 2007]
 *
 * Kyle D.
 * D version 16-Oct-2015
 * Bracketing taken from ridder.d
 */
module nm.brent;

import std.string;
import std.math;
import std.algorithm;
import ntypes.complex;

/**
 * Locate a root of f(x) known to lie between x1 and x2. The method
 * is guaranteed to converge (according to Brent) given the initial 
 * x values bound the solution. The method uses root bracketing,
 * bisection and inverse quadratic interpolation.
 *
 * Params:
 *    f: user-supplied function f(x)
 *    x1: first end of range
 *    x2: other end of range
 *    tol: minimum size for range
 *
 * Returns:
 *    b, a point near the root.
 */
@nogc
T solve(alias f, T)(T x1, T x2, double tol=1.0e-9) 
    if (is(typeof(f(0.0)) == double) ||
        is(typeof(f(0.0)) == float) ||
        is(typeof(f(Complex!double(0.0))) == Complex!double))
{
    const int ITMAX = 100;           // maximum allowed number of iterations
    const double EPS=double.epsilon; // Machine floating-point precision
    T a = x1;
    T b = x2;
    T fa = f(a);
    T fb = f(b);
    if (abs(fa) == 0.0) { return a; }
    if (abs(fb) == 0.0) { return b; }
    if (fa * fb > 0.0) {
        // Supplied bracket does not encompass zero of the function.
        string msg = "Root must be bracketed by x1 and x2";
        debug { msg ~= format("\nx1=%g f(x1)=%g x2=%g f(x2)=%g\n", x1, fa, x2, fb); } 
        throw new Exception(msg);
    }
    T c = b;
    T fc = fb;
    // Define d, e outside the loop body so that
    // they don't get initialized to nan each pass.
    T d, e;
    foreach (iter; 0 .. ITMAX) {
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            // On first pass fc==fb and we expect to enter here.
            c = a;
            fc = fa;
            e = d = b-a;
        }
        if (abs(fc) < abs(fb)) {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        // Convergence check
        T tol1 = 2.0*EPS*abs(b)+0.5*tol;
        T xm = 0.5*(c-b);
        if (abs(xm) <= tol1 || fb == 0.0) {
            // If converged, let's return the best estimate
            return b;
        }
        if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
            // Attempt inverse quadratic interpolation for new bound
            T p, q;
            T s = fb/fa;
            if ( a == c ) {
                p = 2.0*xm*s;
                q = 1.0-s;
            } else {
                q = fa/fc;
                T r = fb/fc;
                p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q = (q-1.0)*(r-1.0)*(s-1.0);
            }
            // Check whether quadratic interpolation is in bounds
            if (p > 0.0) { q = -q; }
            p = abs(p);
            T min1 = 3.0*xm*q-abs(tol1*q);
            T min2 = abs(e*q);
            if (2.0*p < (min1 < min2 ? min1 : min2) ) {
                // Accept interpolation
                e = d;
                d = p/q;
            } else {
                // Interpolation failed, use bisection
                d = xm;
                e = d;
            }
        } else {
            // Bounds decreasing too slowly, use bisection
            d = xm;
            e = d;
        }
        // Move last guess to a
        a = b;
        fa = fb;
        // Evaluate new trial root
        if (abs(d) > tol1) {
            b += d;
        } else {
            b += copysign(tol1, xm);
        }
        fb = f(b);      
    }
    throw new Exception("Maximum number of iterations exceeded in brent.solve()");
} // end solve()


version(brent_test) {
    import util.msg_service;
    import std.conv;
    import nm.number;
    int main() {
        @nogc number test_fun_1(number x) {
            return pow(x,3) + pow(x,2) - 3*x - 3;
        }
        @nogc number test_fun_2(number x, number a) {
            return a*x + sin(x) - exp(x);
        }
        assert(fabs(solve!(test_fun_1,number)(to!number(1), to!number(2)) - 1.732051) < 1.0e-5,
               failedUnitTest());
        number my_a = 3.0;
        @nogc number test_fun_3 (number x) {
            return test_fun_2(x, my_a);
        }
        assert(fabs(solve!(test_fun_3,number)(to!number(0), to!number(1)) - 0.3604217) < 1.0e-5,
               failedUnitTest());

        return 0;
    }
}

