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
double solve(alias f)(double x1, double x2, double tol=1.0e-9) 
    if ( is(typeof(f(0.0)) == double) || is(typeof(f(0.0)) == float) ) {
        const int ITMAX = 100;           // maximum allowed number of iterations
        const double EPS=double.epsilon; // Machine floating-point precision
        double a = x1;
	double b = x2;
	double fa = f(a);
	double fb = f(b);
	if ( abs(fa) == 0.0 ) return a;
        if ( abs(fb) == 0.0 ) return b;
	if ( fa * fb > 0.0 ) {
	    // Supplied bracket does not encompass zero of the function.
	    string msg = "Root must be bracketed by x1 and x2\n";
	    msg ~= format("x1=%g f(x1)=%g x2=%g f(x2)=%g\n", x1, fa, x2, fb); 
     	    throw new Exception(msg);
	}
	double c = b;
	double fc = fb;
        for ( int iter=0; iter<ITMAX; iter++ ) {
	    double d, e, tol1, xm;
            if ( (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) ) {
                c = a;
		fc = fa;
		e = d = b-a;
	    }
            if ( abs(fc) < abs(fb) ) {
                a = b;
		b = c;
		c = a;
		fa = fb;
		fb = fc;
		fc = fa;
	    }
	    tol1 = 2.0*EPS*abs(b)+0.5*tol;  // Convergence check
	    xm = 0.5*(c-b);
	    if ( abs(xm) <= tol1 || fb == 0.0 ) return b;   // if converged let's return the best estimate
	    if ( abs(e) >= tol1 && abs(fa) > abs(fb) ) {
		double p, q, r, s;
	        s = fb/fa;         // Attempt inverse quadratic interpolation for new bound
	        if ( a == c ) {
	            p = 2.0*xm*s;
	            q = 1.0-s;
	        }
	        else {
	            q = fa/fc;
		    r = fb/fc;
		    p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
		    q = (q-1.0)*(r-1.0)*(s-1.0);
		}
		if ( p > 0.0 )
	            q = -q;    // Check whether quadratic interpolation is in bounds
		p = abs(p);
		double min1 = 3.0*xm*q-abs(tol1*q);
		double min2 = abs(e*q);
		if (2.0*p < (min1 < min2 ? min1 : min2) ) {
	            e = d;     // Accept interpolation
		    d = p/q;
		}
		else {
	            d = xm;  // else Interpolation failed, use bisection
	            e = d;
		}
	    }
	    else {
                d = xm;   // Bounds decreasing too slowly, use bisection
		e = d;
	    }
	    a = b;       // Move last guess to a
	    fa = fb;
	    if ( abs(d) > tol1 )  // Evaluate new trial root
	        b += d;
	    else {
	        b += copysign(tol1, xm);
	    }
	    fb = f(b);	    
	}
	throw new Exception("Maximum number of iterations exceeded in brent.d");
    }


version(brent_test) {
    import util.msg_service;
    int main() {
	double test_fun_1(double x) {
	    return pow(x,3) + pow(x,2) - 3*x - 3;
	}
	double test_fun_2(double x, double a) {
	    return a*x + sin(x) - exp(x);
	}
	assert(abs(solve!test_fun_1(1, 2) - 1.732051) < 1.0e-5, failedUnitTest());
	double my_a = 3.0;
	auto test_fun_3 = delegate (double x) { return test_fun_2(x, my_a); };
	assert(abs(solve!test_fun_3(0, 1) - 0.3604217) < 1.0e-5, failedUnitTest());

	return 0;
    }
}

