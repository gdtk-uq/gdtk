/**
 * ridder.d
 *
 * Solve a nonlinear equation f(x)=0 using the method of Ridder.
 *
 * Peter J. 
 * mech3750 demo code 12-Mar-2014
 * added bracketing 19-Mar-2014
 * D version 13-Jun-2014
 */
module ridder;
import std.math;

/**
 * Locate a root of f(x) by subdividing the original range,
 * assuming a linear model of the underlying transformed function.
 *
 * Params:
 *    f: user-supplied function f(x)
 *    x1: first end of range
 *    x2: other end of range
 *    tol: minimum size for range
 *
 * Returns:
 *    x, a point near the root.
 */
double solve(alias f)(double x1, double x2, double tol=1.0e-9) 
    if ( is(typeof(f(0.0)) == double) || is(typeof(f(0.0)) == float) )
{
    double x3, f3;
    double x4 = x1; // So that g++ doesn't warn on maybe unitialized.
    double f4, eps;

    double f1 = f(x1); 
    double f2 = f(x2);
    if ( abs(f1) == 0.0 ) return x1;
    if ( abs(f2) == 0.0 ) return x2;
    if ( x1 == x2 ) {
	throw new Exception("Bad initial range given to bracket.");
    }
    if ( f1 * f2 > 0.0 ) {
	throw new Exception("Range does not clearly bracket a root.");
    }
    while ( abs(x2 - x1) > tol ) {
	x3 = 0.5*(x1+x2);
	f3 = f(x3);
	if ( f3 == 0.0 ) return x3;
	eps = (f3 + copysign(sqrt(f3*f3-f1*f2),f2))/f2;
	x4 = x3 - f3*eps*(x1-x3)/(f1 - eps*f3);
	f4 = f(x4);
	if ( f4 == 0.0 ) return x4;
	// Contract the bracket.
	if ( f3*f2 < 0.0 ) {
	    if ( f4*f2 < 0.0 ) {
		x1 = x4; f1 = f4;
	    } else {
		x1 = x3; f1 = f3;
		x2 = x4; f2 = f4;
	    }
	} else {
	    if ( f4*f1 < 0.0 ) {
		x2 = x4; f2 = f4;
	    } else {
		x1 = x4; f1 = f4;
		x2 = x3; f2 = f3;
	    }
	}
    } // end while
    return x4;
} // end solve()


/**
 * Bracket a root of f(x).
 *
 * Params:
 *    f: user-supplied function f(x)
 *    x1: first end of range
 *    x2: other end of range
 *
 * Returns:
 *    0: successfully bracketed a root
 *   -1: failed to bracket a root
 *
 * On return (x1, x2) should bracketing the root.
 */
int bracket(alias f)(ref double x1, ref double x2,
                     int max_try=50, double factor=1.6)
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
            f1 = f(x1);
        } else {
            x2 += factor * (x2 - x1);
            f2 = f(x2);
	}
    }
    // If we leave the loop here, we were unsuccessful.
    return -1;
} // end bracket()


unittest {
    double test_fun_1(double x) {
	return pow(x,3) + pow(x,2) - 3*x - 3;
    }
    double test_fun_2(double x, double a) {
	return a*x + sin(x) - exp(x);
    }
    assert( abs(solve!test_fun_1(1, 2) - 1.732051) < 1.0e-5 );
    double my_a = 3.0;
    auto test_fun_3 = delegate (double x) { return test_fun_2(x, my_a); };
    assert( abs(solve!test_fun_3(0, 1) - 0.3604217) < 1.0e-5 );
    double x1 = 0.4;
    double x2 = 0.5;
    assert( bracket!test_fun_3(x1, x2) == 0 );
}
