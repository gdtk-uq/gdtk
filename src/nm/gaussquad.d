/**
 * gaussquad.d
 *
 * Adaptive integration of a function using Gauss quadrature rules.
 *
 * Author: Peter J.
 * Version: 2014-06-21: adapted from newtoncotes.d
 */

import std.math;

// Rule coefficients computed by Maxima. 
// See the scripts in directory ./notes

double[] xs3 = [-7.7459666924148337704e-1,
		0.0,
		7.7459666924148337704e-1];
double[] ws3 = [5.5555555555555555556e-1,
		8.8888888888888888889e-1,
		5.5555555555555555556e-1];

double[] xs4 = [-8.6113631159405257523e-1, -3.399810435848562648e-1,
		3.399810435848562648e-1, 8.6113631159405257523e-1];
double[] ws4 = [3.4785484513745385738e-1, 6.5214515486254614263e-1,
		6.5214515486254614263e-1, 3.4785484513745385738e-1];

/**
 * Adaptively apply Gauss 4- and 3-point quadrature rules to the segment [a,b].
 *
 * Params:
 *     f: user-supplied function, f(x)
 *     a, b: range of integration
 *     tol: maximum difference between rules, above which the range is split
 *
 * Returns: 
 *     integral of f(x) from a to b.
 */
double apply_gauss_rule(alias f)(double a, double b, double[] xs, double[] ws)
    if ( is(typeof(f(0.0)) == double) || is(typeof(f(0.0)) == float) )
{
    double xmid = 0.5 * (b + a);
    double xrange2 = 0.5 * (b - a);
    size_t N = xs.length;
    double result = 0.0;
    foreach(i; 0 .. N) {
	double x = xmid + xs[i] * xrange2; // transform from range -1.0 +1.0
	result += ws[i] * f(x);
    }
    return xrange2 * result;
}

/**
 * Adaptively apply Gauss 4- and 3-point quadrature rules to the segment [a,b].
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
    double I3 = apply_gauss_rule!f(a, b, xs3, ws3);
    double I4 = apply_gauss_rule!f(a, b, xs4, ws4);
    double I;
    if ( abs(I4 - I3) > tol ) {
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

