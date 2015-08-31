/**
 * linesearch.d
 *
 * Implementation of an algorithm for optimization from Gerald and Wheatley.
 *
 * A class demo for mech2700, originally in Python.
 * Author: PJ, 14-Jun-2014 
*/

module linesearch;
import std.stdio;

/**
 * Returns the bracket xL,xR containing the minimum of the function f.
 *
 * Params:
 *     f: a user supplied function
 *     a,b: the bracket containing a minimum
 *     tol: the final size of the bracket.
 *          It should not be set too small.
 */
void minimize(alias f)(ref double a, ref double b, double tol=1.0e-4)
    if ( is(typeof(f(0.0)) == double) || is(typeof(f(0.0)) == float) )
{
    double r = 0.618034;
    double xL = a + (1-r)*(b-a);
    double xR = a + r*(b-a);
    double FL = f(xL);
    double FR = f(xR);

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


unittest {
    import std.math;
    double fdemo(double x) {
        return exp(x) + 2.0 - cos(x);
    }
    double a = -3;
    double b = 1;
    minimize!fdemo(a, b, 1.0e-6);
    double xminimum = -0.588534;
    assert( abs(a - xminimum) < 1.0e-4 );
    assert( abs(b - xminimum) < 1.0e-4 );
}
