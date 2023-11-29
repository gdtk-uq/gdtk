/**
 * linesearch.d
 *
 * Implementation of an algorithm for optimization from Gerald and Wheatley.
 *
 * A class demo for mech2700, originally in Python.
 * Author: PJ, 14-Jun-2014 
*/

module nm.linesearch;
import ntypes.complex;

/**
 * Returns the bracket xL,xR containing the minimum of the function f.
 *
 * Params:
 *     f: a user supplied function
 *     a,b: the bracket containing a minimum
 *     tol: the final size of the bracket.
 *          It should not be set too small.
 */
void minimize(alias f, T)(ref T a, ref T b, double tol=1.0e-4)
    if (is(typeof(f(0.0)) == double) ||
        is(typeof(f(0.0)) == float) ||
        is(typeof(f(Complex!double(0.0))) == Complex!double))
{
    T r = 0.618034;
    T xL = a + (1-r)*(b-a);
    T xR = a + r*(b-a);
    T FL = f(xL);
    T FR = f(xR);
    while ((xR - xL) > tol) {
        if (FR > FL) {
            b = xR;
            xR = xL;
            FR = FL;
            xL = a + (1.0-r)*(b-a);
            FL = f(xL);
        } else {
            a = xL;
            xL = xR;
            FL = FR;
            xR = a + r*(b-a);
            FR = f(xR);
        }
    }
    a = xL;
    b = xR;
 } // end minimize()


version(linesearch_test) {
    import std.math;
    import util.msg_service;
    import nm.number;
    int main() {
        number fdemo(number x) {
            return exp(x) + 2.0 - cos(x);
        }
        number a = -3;
        number b = 1;
        minimize!(fdemo,number)(a, b, 1.0e-6);
        number xminimum = -0.588534;
        assert(fabs(a - xminimum) < 1.0e-4, failedUnitTest());
        assert(fabs(b - xminimum) < 1.0e-4, failedUnitTest());

        return 0;
    }
}
