/**
 * newtoncotes.d
 *
 * Adaptive quadrature using Newton-Cotes 5- and 3-point rules.
 *
 * Author: PA Jacobs, School of Mechanical and Mining Engineering, UQ
 * Version: 2014-06-20, Adapted from MECH2700 demonstration code written back in 2003.
 */
module nm.newtoncotes;
import std.math;
import ntypes.complex;

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
T integrate(alias f, T)(T a, T b, double tol=1.0e-5)
    if ( is(typeof(f(0.0)) == double) ||
         is(typeof(f(0.0)) == float)  ||
         is(typeof(f(Complex!double(0.0))) == Complex!double))
{
    T dx = b - a;
    T f0 = f(a);
    T f1 = f(a + 0.25 * dx);
    T f2 = f(a + 0.5 * dx);
    T f3 = f(a + 0.75 * dx);
    T f4 = f(b);
    T I2 = dx/6.0 * (f0 + 4 * f2 + f4);
    T I4 = dx/90.0 * (7*f0 + 32*f1 + 12*f2 + 32*f3 + 7*f4);
    T I;
    if (abs(I4 - I2) > tol) {
        T mid = 0.5 * (a + b);
        I = integrate!(f,T)(a, mid, tol/2.0) + integrate!(f,T)(mid, b, tol/2.0);
    } else {
        I = I4;
    }
    return I;
}


version(newtoncotes_test) {
    import util.msg_service;
    import std.conv;
    import nm.number;
    int main() {
        const number one = 1.0;
        const number zero = 0.0;
        number fun1(number x) { return abs(x) < one ? sqrt(one - x*x): zero; }
        number fun2(number x) { return 1.0 / (1.0 + x * x); }
        assert(approxEqualNumbers(to!number(PI/4), integrate!(fun1,number)(zero, one), 1.0e-6),
               failedUnitTest()); 
        assert(approxEqualNumbers(to!number(PI/4), integrate!(fun2,number)(zero, one), 1.0e-6),
               failedUnitTest()); 
        return 0;
    }
}
