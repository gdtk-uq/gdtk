/**
 * gaussquad.d
 *
 * Adaptive integration of a function using Gauss quadrature rules.
 *
 * Author: Peter J.
 * Version: 2014-06-21: adapted from newtoncotes.d
 */

module nm.gaussquad;
import std.math;
import ntypes.complex;

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
T apply_gauss_rule(alias f, T)(T a, T b, double[] xs, double[] ws)
    if ( is(typeof(f(0.0)) == double) ||
         is(typeof(f(0.0)) == float)  ||
         is(typeof(f(Complex!double(0.0))) == Complex!double))
{
    T xmid = 0.5 * (b + a);
    T xrange2 = 0.5 * (b - a);
    size_t N = xs.length;
    T result = 0.0;
    foreach(i; 0 .. N) {
        T x = xmid + xs[i] * xrange2; // transform from range -1.0 +1.0
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
T integrate(alias f, T)(T a, T b, double tol=1.0e-5)
    if ( is(typeof(f(0.0)) == double) ||
         is(typeof(f(0.0)) == float)  ||
         is(typeof(f(Complex!double(0.0))) == Complex!double))
{
    T I3 = apply_gauss_rule!(f,T)(a, b, xs3, ws3);
    T I4 = apply_gauss_rule!(f,T)(a, b, xs4, ws4);
    T I;
    if ( abs(I4 - I3) > tol ) {
        T mid = 0.5 * (a + b);
        I = integrate!(f,T)(a, mid, tol/2.0) + integrate!(f,T)(mid, b, tol/2.0);
    } else {
        I = I4;
    }
    return I;
}

version(gaussquad_test) {
    import util.msg_service;
    import std.conv;
    import nm.number;
    int main() {
        number fun1(number x) { return abs(x) < to!number(1.0) ? sqrt(1.0 - x*x): to!number(0.0); }
        number fun2(number x) { return 1.0 / (1.0 + x * x); }
        assert(approxEqualNumbers(to!number(PI/4), integrate!fun1(to!number(0.0), to!number(1.0)), 1.0e-6), failedUnitTest()); 
        assert(approxEqualNumbers(to!number(PI/4), integrate!fun2(to!number(0.0), to!number(1.0)), 1.0e-6), failedUnitTest());
        
        return 0;
    }
}

