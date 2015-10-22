/**
 * Bracket a root of f(x).
 *
 * Params:
 *    f: user-supplied function f(x)
 *    x1: first end of range
 *    x2: other end of range
 *  x1min: limits the value of x1 so that it does not become lower than a certain value
 *
 * Returns:
 *    0: successfully bracketed a root
 *   -1: failed to bracket a root
 *
 * On return (x1, x2) should bracketing the root.
 */

module bracket;
import std.math;
import std.algorithm;

int bracket(alias f)(ref double x1, ref double x2,
                     int max_try=50, double factor=1.6, double x1_min = -1.0e99)
    if ( is(typeof(f(0.0)) == double) || is(typeof(f(0.0)) == float) )
{
    if ( x1 == x2 ) {
    throw new Exception("Bad initial range given to bracket.");
    }
    double f1 = f(x1);
    double f2 = f(x2);
    for ( int i = 0; i < max_try; ++i ) {
        if ( f1*f2 < 0.0 ) return 0; // we have success
        if ( abs(f1) < abs(f2) ) {
            x1 += factor * (x1 - x2);
            x1 = max(x1_min, x1);//prevent the bracket from being expanded beyond a specified domain
            f1 = f(x1);
        } else {
            x2 += factor * (x2 - x1);
            f2 = f(x2);
    }
    }
    // If we leave the loop here, we were unsuccessful.
    return -1;
 } // end bracket()
