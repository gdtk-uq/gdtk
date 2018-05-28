/**
 * linesearch.d
 *
 * Implementation of an algorithm for optimization from Gerald and Wheatley.
 *
 * A class demo for mech2700, originally in Python.
 * Author: PJ, 14-Jun-2014 
*/

module nm.linesearch;
import std.stdio;
import nm.complex;
import nm.number;

/**
 * Returns the bracket xL,xR containing the minimum of the function f.
 *
 * Params:
 *     f: a user supplied function
 *     a,b: the bracket containing a minimum
 *     tol: the final size of the bracket.
 *          It should not be set too small.
 */
void minimize(alias f)(ref number a, ref number b, double tol=1.0e-4)
    if ( is(typeof(f(0.0)) == double) ||
         is(typeof(f(0.0)) == float)  ||
         is(typeof(f(Complex!double(0.0))) == Complex!double))
{
    number r = 0.618034;
    number xL = a + (1-r)*(b-a);
    number xR = a + r*(b-a);
    number FL = f(xL);
    number FR = f(xR);
    while ( (xR - xL) > tol ) {
        if ( FR > FL ) {
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
        //writeln("xL=", xL, " xR=", xR, " FL=", FL, " FR=", FR);
    }
    a = xL;
    b = xR;
 } // end minimize()


version(linesearch_test) {
    import std.math;
    import util.msg_service;
    int main() {
        number fdemo(number x) {
            return exp(x) + 2.0 - cos(x);
        }
        number a = -3;
        number b = 1;
        minimize!fdemo(a, b, 1.0e-6);
        number xminimum = -0.588534;
        assert(fabs(a - xminimum) < 1.0e-4, failedUnitTest());
        assert(fabs(b - xminimum) < 1.0e-4, failedUnitTest());

        return 0;
    }
}
