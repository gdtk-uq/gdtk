/** nelmin_demo.d
 *
 * Try out the Nelder-Mead Simplex optimizer.
 * 
 * Author: Peter J
 * Version: 2014-06-14, Adapted from nelmin_test.cxx
 */

import std.stdio;
import std.math;
import nelmin;

/** Test objective function 1.
 *  x is expected to be a list of ccordinates.
 *  Returns a single float value.
 */
double testFunction1(double[] x)
{
    double sum = 0.0;
    foreach (elem; x) sum += (elem - 1.0) * (elem - 1.0); 
    return sum;
}

// Test objective function 2.
// Example 3.3 from Olsson and Nelson.
double testFunction2(double[] x)
{
    double x1 = x[0]; double x2 = x[1];   // rename to match the paper
    if ( (x1 * x1 + x2 * x2) > 1.0 ) {
	return 1.0e38;
    } else {
	double yp = 53.69 + 7.26 * x1 - 10.33 * x2 + 7.22 * x1 * x1 
	    + 6.43 * x2 * x2 + 11.36 * x1 * x2;
	double ys = 82.17 - 1.01 * x1 - 8.61 * x2 + 1.40 * x1 * x1 
	    - 8.76 * x2 * x2 - 7.20 * x1 * x2;
	return -yp + fabs(ys - 87.8);
    }
}
 
// Test objective function 3.
// Example 3.5 from Olsson and Nelson; least-squares.
double testFunction3(double[] z)
{
    double[] x = [0.25, 0.50, 1.00, 1.70, 2.00, 4.00];
    double[] y = [0.25, 0.40, 0.60, 0.58, 0.54, 0.27];
    double a1 = z[0]; double a2 = z[1]; 
    double alpha1 = z[2]; double alpha2 = z[3];
    double sum_residuals = 0.0;
    foreach (i; 0 .. 6) {
	double t = x[i];
	double eta = a1 * exp(alpha1 * t) + a2 * exp(alpha2 * t);
	double r = y[i] - eta;
	sum_residuals += r * r;
    }
    return sum_residuals;
}
   
// -------------------------------------------------------------------

void main() {
    writeln("Begin nelmin demostration...");

    writeln("---------------------------------------------------");
    writeln("test 1: simple quadratic with zero at (1,1,...)");
    auto x = [0.0, 0.0, 0.0];
    double[] dx;
    double fx;
    int nfe, nres;
    bool conv_flag = minimize!testFunction1(x, fx, nfe, nres, dx);
    writeln("x = ", x, " fx = ", fx);
    writeln("convergence-flag = ", conv_flag);
    writeln("number-of-fn-evaluations = ", nfe);
    writeln("number-of-restarts = ", nres);

    writeln("---------------------------------------------------");
    writeln("test 2: Example 3.3 in Olsson and Nelson f(0.811,-0.585)=-67.1");
    double[] x2 = [0.0, 0.0];
    double[] dx2 = [0.5, 0.5];
    conv_flag = minimize!testFunction2(x2, fx, nfe, nres, dx2, 1.0e-4);
    writeln("x = ", x2, " fx = ", fx);
    writeln("convergence-flag = ", conv_flag);
    writeln("number-of-fn-evaluations = ", nfe);
    writeln("number-of-restarts =", nres);

    writeln("---------------------------------------------------");
    writeln("test 3: Example 3.5 in Olsson and Nelson, nonlinear least-squares");
    writeln("f(1.801, -1.842, -0.463, -1.205)=0.0009");
    double[] x3 = [1.0, 1.0, -0.5, -2.5];
    double[] dx3 = [0.1, 0.1, 0.1, 0.1];
    conv_flag = minimize!testFunction3(x3, fx, nfe, nres, dx3, 1.0e-9, 800);
    writeln("x = ", x3, " fx = ", fx);
    writeln("convergence-flag = ", conv_flag);
    writeln("number-of-fn-evaluations = ", nfe);
    writeln("number-of-restarts = ", nres);

    writeln("---------------------------------------------------");
    writeln("Done.");
} // end main()
