/**
 * rungekutta.d
 *
 * Ordinary differential equation integration.
 *
 * Author: Peter J.
 * Version: 2014-Jun-15, adapted from the mech2700 class example
 *          2014-Jul-09, preallocate work arrays and pass them in.
 *
 */

module rungekutta;
import std.math;

/**
 * Allocate workspace arrays for the ODE stepper.
 *
 * To avoid the allocation of the intermediate arrays at every update,
 * we get the user to preallocate them with this function.
 */
double[][] allocate_rk45_workspace(uint n)
{
    double[][] workspace_arrays;
    workspace_arrays.length = 7;
    foreach(ref vect; workspace_arrays) vect.length = n;
    return workspace_arrays;
} // end allocate_rk45_workspace()

/**
 * Steps the set of ODEs by the Runge-Kutta-Fehlberg method.
 *
 * Params:
 *     f: a callable function that returns the derivative of y wrt t
 *        The signature of this function is f(t, y) where
 *        t is a float value, y is an array of double values.
 *     t0: is the starting value of the independent variable
 *     h: the requested step size
 *     y0: an array of starting values for the dependent variables
 *         It is assumed that the y-elements are indexed 0 .. n-1
 *         where n = y0.length
 *     y1: an array of final values of the dependent variables
 *     err: estimates of the errors in the values of y1
 *
 * Returns:
 *     the final value of the dependent variable
 */
double rkf45_step(alias f)(double t0, double h, double[] y0,
			   ref double[] y1, ref double[] err,
			   ref double[][] work_arrays)
    if ( is(typeof(f(0.0, [0.0,0.0])) == double[]) )
{
    // Assuming a system of equations, we need arrays for the intermediate data.
    double[] k1 = work_arrays[1];
    if ( k1.length != y0.length ) {
	throw new Exception("Array lengths don't match the workspace and.");
    }
    double[] k2 = work_arrays[2];
    double[] k3 = work_arrays[3]; 
    double[] k4 = work_arrays[4];
    double[] k5 = work_arrays[5];
    double[] k6 = work_arrays[6];
    double[] ytmp = work_arrays[0];
    // Build up the sample point information as per the text book descriptions.
    // We assign the result of intermediate array expressions to ytmp
    // because that's needed for D.
    k1[] = f(t0, y0);
    ytmp[] = y0[] + 0.25*h*k1[];
    k2[] = f(t0 + h/4.0, ytmp);
    ytmp[] = y0[] + 3.0*h*k1[]/32.0 + 9.0*h*k2[]/32.0;
    k3[] = f(t0 + 3.0*h/8.0, ytmp);
    ytmp[] = y0[] + 1932.0*h*k1[]/2197.0 - 7200.0*h*k2[]/2197.0 + 
	7296.0*h*k3[]/2197.0;
    k4[] = f(t0 + 12.0*h/13.0, ytmp);
    ytmp[] = y0[] + 439.0*h*k1[]/216.0 - 8.0*h*k2[] + 
	3680.0*h*k3[]/513.0 - 845.0*h*k4[]/4104.0;
    k5[] = f(t0 + h, ytmp);
    ytmp[] = y0[] - 8.0*h*k1[]/27.0 + 2.0*h*k2[] - 
	3544.0*h*k3[]/2565.0 + 1859.0*h*k4[]/4104.0 - 11.0*h*k5[]/40.0;
    k6[] = f(t0 + h/2.0, ytmp);
    // Now, do the integration as a weighting of the sampled data.
    y1[] = y0[] + 16.0*h*k1[]/135.0 + 6656.0*h*k3[]/12825.0 + 
	28561.0*h*k4[]/56430.0 - 9.0*h*k5[]/50.0 + 2.0*h*k6[]/55.0;
    err[] = h*k1[]/360.0 - 128.0*h*k3[]/4275.0 - 2197.0*h*k4[]/75240.0 + 
	h*k5[]/50.0 + 2.0*h*k6[]/55.0;
    foreach(ref e; err) e = fabs(e);
    return t0 + h;
} // end rkf45_step()

unittest {
    double[] testSystem1(double t, double[] x)
    {
	double dx0dt =  -8.0/3.0*x[0] -  4.0/3.0*x[1] +     x[2] + 12.0;
	double dx1dt = -17.0/3.0*x[0] -  4.0/3.0*x[1] +     x[2] + 29.0;
	double dx2dt = -35.0/3.0*x[0] + 14.0/3.0*x[1] - 2.0*x[2] + 48.0;
	return [dx0dt, dx1dt, dx2dt];
    }
    double[] solution1(double t)
    {
	double x = exp(-3.0*t)/6.0*(6.0-50.0*exp(t)+10.0*exp(2.0*t)+34.0*exp(3.0*t));
	double y = exp(-3.0*t)/6.0*(12.0-125.0*exp(t)+40.0*exp(2.0*t)+73.0*exp(3.0*t));
	double z = exp(-3.0*t)/6.0*(14.0-200.0*exp(t)+70.0*exp(2.0*t)+116.0*exp(3.0*t));
	return [x, y, z];
    }
    double[] x0=[0.0, 0.0, 0.0];
    double[] x1=x0.dup;
    double[] err=x0.dup;
    auto work = allocate_rk45_workspace(3);
    double t1 = rkf45_step!(testSystem1)(0.0, 0.2, x0, x1, err, work);
    assert(approxEqual(x1, solution1(t1)), "Single step of rkf45.");
}
