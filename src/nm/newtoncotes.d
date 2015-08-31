/**
 * newtoncotes.d
 *
 * Adaptive quadrature using Newton-Cotes 5- and 3-point rules.
 *
 * Author: PA Jacobs, School of Mechanical and Mining Engineering, UQ
 * Version: 2014-06-20, Adapted from MECH2700 demonstration code written back in 2003.
 */

import std.math;

/**
 * Apply Newton-Cotes 5- and 3-point quadrature rules to the segment [a,b].
 *
 * Params:
 *     f: user-supplied function, f(x)
 *     a, b: range of integration
 *     tol: maximum difference between rules, above which the range is split
 *
 * Returns: 
 *     integral of f(x) from a to b.
 */
double integrate(alias f)(double a, double b, double tol=1.0e-5)
    if ( is(typeof(f(0.0)) == double) || is(typeof(f(0.0)) == float) )
{
    double dx = b - a;
    double f0 = f(a);
    double f1 = f(a + 0.25 * dx);
    double f2 = f(a + 0.5 * dx);
    double f3 = f(a + 0.75 * dx);
    double f4 = f(b);
    double I2 = dx/6.0 * (f0 + 4 * f2 + f4);
    double I4 = dx/90.0 * (7*f0 + 32*f1 + 12*f2 + 32*f3 + 7*f4);
    double I;
    if ( abs(I4 - I2) > tol ) {
        double mid = 0.5 * (a + b);
	I = integrate!f(a, mid, tol/2.0) + integrate!f(mid, b, tol/2.0);
    } else {
        I = I4;
    }
    return I;
}


unittest{
    import std.math;
    double fun1(double x) { return abs(x) < 1.0 ? sqrt(1.0 - x*x): 0.0; }
    double fun2(double x) { return 1.0 / (1.0 + x * x); }
    assert(approxEqual(PI/4, integrate!fun1(0.0, 1.0)), "integrate fun1"); 
    assert(approxEqual(PI/4, integrate!fun2(0.0, 1.0)), "integrate fun2"); 
}
