/**
 * rungekutta.d
 *
 * Ordinary differential equation integration.
 *
 * Author: Peter J.
 * Version: 2014-Jun-15, adapted from the mech2700 class example
 *          2014-Jul-09, preallocate work arrays and pass them in.
 *          2018-May-26, work with double or complex numbers
 *          2018-May-30, accept the type of the dependent variables as a parameter
 */

module nm.rungekutta;

import std.math;
import ntypes.complex;

/**
 * Allocate workspace arrays for the ODE stepper.
 *
 * To avoid the allocation of the intermediate arrays at every update,
 * we get the user to preallocate them with this function.
 */
T[][] allocate_rk45_workspace(T)(uint n)
{
    T[][] workspace_arrays;
    workspace_arrays.length = 7;
    foreach(ref vect; workspace_arrays) { vect.length = n; }
    return workspace_arrays;
}

/**
 * Steps the set of ODEs by the Runge-Kutta-Fehlberg method.
 *
 * Params:
 *     f: a callable function that returns the derivative of y wrt t
 *        The signature of this function is f(t, y) where
 *        t is a float value, y is an array of number values.
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
T rkf45_step(alias f, T)(T t0, T h, T[] y0, ref T[] y1, ref T[] err,
                         ref T[][] work_arrays)
    if (is(typeof(f(0.0, [0.0,])) == double[]) ||
        is(typeof(f(0.0, [0.0,])) == float[]) ||
        is(typeof(f(Complex!double(0.0), [Complex!double(0.0),Complex!double(0.0)]))
           == Complex!double[]))
{
    // Assuming a system of equations, we need arrays for the intermediate data.
    size_t n = y0.length;
    T[] k1 = work_arrays[1];
    if ( k1.length != y0.length ) {
        throw new Exception("Array lengths don't match the workspace and.");
    }
    T[] k2 = work_arrays[2];
    T[] k3 = work_arrays[3];
    T[] k4 = work_arrays[4];
    T[] k5 = work_arrays[5];
    T[] k6 = work_arrays[6];
    T[] ytmp = work_arrays[0];
    // Build up the sample point information as per the text book descriptions.
    // We assign the result of intermediate array expressions to ytmp
    // because that's needed for D.
    k1[] = f(t0, y0);
    foreach (i; 0 .. n) { ytmp[i] = y0[i] + 0.25*h*k1[i]; }
    k2[] = f(t0 + h/4.0, ytmp);
    foreach (i; 0 .. n) { ytmp[i] = y0[i] + 3.0*h*k1[i]/32.0 + 9.0*h*k2[i]/32.0; }
    k3[] = f(t0 + 3.0*h/8.0, ytmp);
    foreach (i; 0 .. n) {
        ytmp[i] = y0[i] + 1932.0*h*k1[i]/2197.0 - 7200.0*h*k2[i]/2197.0 +
            7296.0*h*k3[i]/2197.0;
    }
    k4[] = f(t0 + 12.0*h/13.0, ytmp);
    foreach (i; 0 .. n) {
        ytmp[i] = y0[i] + 439.0*h*k1[i]/216.0 - 8.0*h*k2[i] +
            3680.0*h*k3[i]/513.0 - 845.0*h*k4[i]/4104.0;
    }
    k5[] = f(t0 + h, ytmp);
    foreach (i; 0 .. n) {
        ytmp[i] = y0[i] - 8.0*h*k1[i]/27.0 + 2.0*h*k2[i] -
            3544.0*h*k3[i]/2565.0 + 1859.0*h*k4[i]/4104.0 - 11.0*h*k5[i]/40.0;
    }
    k6[] = f(t0 + h/2.0, ytmp);
    // Now, do the integration as a weighting of the sampled data.
    foreach (i; 0 .. n) {
        y1[i] = y0[i] + 16.0*h*k1[i]/135.0 + 6656.0*h*k3[i]/12825.0 +
            28561.0*h*k4[i]/56430.0 - 9.0*h*k5[i]/50.0 + 2.0*h*k6[i]/55.0;
        err[i] = h*k1[i]/360.0 - 128.0*h*k3[i]/4275.0 - 2197.0*h*k4[i]/75240.0 +
            h*k5[i]/50.0 + 2.0*h*k6[i]/55.0;
    }
    foreach(ref e; err) { e = abs(e); }
    return t0 + h;
} // end rkf45_step()

version(rungekutta_test) {
    import util.msg_service;
    import std.conv;
    import nm.number;
    int main() {
        number[] testSystem1(number t, number[] x)
        {
            number dx0dt =  -8.0/3.0*x[0] -  4.0/3.0*x[1] +     x[2] + 12.0;
            number dx1dt = -17.0/3.0*x[0] -  4.0/3.0*x[1] +     x[2] + 29.0;
            number dx2dt = -35.0/3.0*x[0] + 14.0/3.0*x[1] - 2.0*x[2] + 48.0;
            return [dx0dt, dx1dt, dx2dt];
        }
        number[] solution1(number t)
        {
            number x = exp(-3.0*t)/6.0*(6.0-50.0*exp(t)+10.0*exp(2.0*t)+34.0*exp(3.0*t));
            number y = exp(-3.0*t)/6.0*(12.0-125.0*exp(t)+40.0*exp(2.0*t)+73.0*exp(3.0*t));
            number z = exp(-3.0*t)/6.0*(14.0-200.0*exp(t)+70.0*exp(2.0*t)+116.0*exp(3.0*t));
            return [x, y, z];
        }
        number[] x0=[to!number(0.0), to!number(0.0), to!number(0.0)];
        number[] x1=x0.dup;
        number[] err=x0.dup;
        auto work = allocate_rk45_workspace!number(3);
        number t0 = 0.0;
        number h = 0.2;
        number t1 = rkf45_step!(testSystem1, number)(t0, h, x0, x1, err, work);
        assert(approxEqualNumbers(x1, solution1(t1), 1.0e-5), failedUnitTest());

        return 0;
    }
}
