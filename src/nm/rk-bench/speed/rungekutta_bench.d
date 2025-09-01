/** rungekutta_demo.d
 *
 * Try out the Runge-Kutta ODE stepper as a benchmark program.
 *
 * Author: Peter J.
 * Version: 2014-Jun-15, adapted from the mech2700 class example
 *          2014-Jul-09, preallocate work arrays and pass them in.
 *          2018-May-26, work with double or complex numbers
 *          2018-May-30, accept the type of the dependent variables as a parameter
 *          2022-May-20, Build as a single-source-file program.
 */

import std.stdio;
import std.math;
import core.time;

/**
 * Steps the set of ODEs by the Runge-Kutta-Fehlberg method.
 *
 * Params:
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
 *
 * Types:
 *     T  may be double
 *     TA may be double[] or double[3], etc (a slice or array)
 *     f is a callable function that returns the derivative of y wrt t
 *        The signature of this function is f(t, y) where
 *        t is a float value, y is an array of number values.
 */
double rkf45_step(T, TA, alias f)(T t0, T h, TA y0,
                                  ref TA y1, ref TA err, ref TA[7] work_arrays)
{
    // Assuming a system of equations, we need arrays for the intermediate data.
    TA k1 = work_arrays[1];
    TA k2 = work_arrays[2];
    TA k3 = work_arrays[3];
    TA k4 = work_arrays[4];
    TA k5 = work_arrays[5];
    TA k6 = work_arrays[6];
    TA ytmp = work_arrays[0];
    // Build up the sample point information as per the text book descriptions.
    // We assign the result of intermediate array expressions to ytmp
    // because that's needed for D.
    k1[] = f(t0, y0);
    ytmp[] = y0[] + 0.25*h*k1[];
    k2[] = f(t0 + h/4.0, ytmp);
    ytmp[] = y0[] + 3.0*h*k1[]/32.0 + 9.0*h*k2[]/32.0;
    k3[] = f(t0 + 3.0*h/8.0, ytmp);
    ytmp[] = y0[] + 1932.0*h*k1[]/2197.0 - 7200.0*h*k2[]/2197.0 + 7296.0*h*k3[]/2197.0;
    k4[] = f(t0 + 12.0*h/13.0, ytmp);
    ytmp[] = y0[] + 439.0*h*k1[]/216.0 - 8.0*h*k2[] + 3680.0*h*k3[]/513.0 - 845.0*h*k4[]/4104.0;
    k5[] = f(t0 + h, ytmp);
    ytmp[] = y0[] - 8.0*h*k1[]/27.0 + 2.0*h*k2[] -
        3544.0*h*k3[]/2565.0 + 1859.0*h*k4[]/4104.0 - 11.0*h*k5[]/40.0;
    k6[] = f(t0 + h/2.0, ytmp);
    // Now, do the integration as a weighting of the sampled data.
    y1[] = y0[] + 16.0*h*k1[]/135.0 + 6656.0*h*k3[]/12825.0 +
        28561.0*h*k4[]/56430.0 - 9.0*h*k5[]/50.0 + 2.0*h*k6[]/55.0;
    err[] = h*k1[]/360.0 - 128.0*h*k3[]/4275.0 - 2197.0*h*k4[]/75240.0 +
        h*k5[]/50.0 + 2.0*h*k6[]/55.0;
    foreach(ref e; err) { e = abs(e); }
    return t0 + h;
} // end rkf45_step()



void main()
{
    writeln("Begin demonstration of ODE stepper...");
    /** Test system 1
     * Third-order system with a simple analytic solution.
     * Adapted from section 11.3 in Cheney and Kincaid, 6th ed.
     * Except for the zero-based indexing, the notation is
     * chosen to match that in the text.
     */
    double[] delegate(double t, double[] x) testSystem1 = delegate(double t, double[] x)
        {
         double dx0dt =  -8.0/3.0*x[0] -  4.0/3.0*x[1] +     x[2] + 12.0;
         double dx1dt = -17.0/3.0*x[0] -  4.0/3.0*x[1] +     x[2] + 29.0;
         double dx2dt = -35.0/3.0*x[0] + 14.0/3.0*x[1] - 2.0*x[2] + 48.0;
         return [dx0dt, dx1dt, dx2dt];
        };

    double[] testSystem2(double t, double[] x)
    {
        double dx0dt =  -8.0/3.0*x[0] -  4.0/3.0*x[1] +     x[2] + 12.0;
        double dx1dt = -17.0/3.0*x[0] -  4.0/3.0*x[1] +     x[2] + 29.0;
        double dx2dt = -35.0/3.0*x[0] + 14.0/3.0*x[1] - 2.0*x[2] + 48.0;
        return [dx0dt, dx1dt, dx2dt];
    }

    @nogc
    double[3] testSystem3(double t, double[3] x)
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

    double[] x1, err, errsum;
    x1.length = 3; err.length = 3; errsum.length = 3;
    double[][7] work; foreach (ref w; work) { w.length = 3; }

    writeln("Using a delegate with slices.");
    double t0 = 0.0;
    double t1 = 0.0;
    int Nstep = 10000;
    double h = 1.0/Nstep;
    double[] x0 = [0.0, 0.0, 0.0];
    auto start_time = MonoTime.currTime;
    foreach (i; 0 .. Nstep) {
        t1 = rkf45_step!(double, double[], testSystem1)(t0, h, x0, x1, err, work);
        x0[] = x1[];
        t0 = t1;
    }
    auto elapsed_time = MonoTime.currTime - start_time;
    writeln("  elapsed_time=", elapsed_time);
    writeln("  x1 = ", x1);
    writeln("  exact = ", solution1(t1));

    writeln("Using a normal function with slices.");
    t0 = 0.0;
    t1 = 0.0;
    x0 = [0.0, 0.0, 0.0];
    start_time = MonoTime.currTime;
    foreach (i; 0 .. Nstep) {
        t1 = rkf45_step!(double, double[], testSystem2)(t0, h, x0, x1, err, work);
        x0[] = x1[];
        t0 = t1;
    }
    elapsed_time = MonoTime.currTime - start_time;
    writeln("  elapsed_time=", elapsed_time);
    writeln("  x1 = ", x1);
    writeln("  exact = ", solution1(t1));

    writeln("Using a @nogc function with fixed array size.");
    t0 = 0.0;
    t1 = 0.0;
    double[3] xx0 = [0.0, 0.0, 0.0];
    double[3] xx1, eerr;
    double[3][7] wwork;
    start_time = MonoTime.currTime;
    foreach (i; 0 .. Nstep) {
        t1 = rkf45_step!(double, double[3], testSystem3)(t0, h, xx0, xx1, eerr, wwork);
        xx0[] = xx1[];
        t0 = t1;
    }
    elapsed_time = MonoTime.currTime - start_time;
    writeln("  elapsed_time=", elapsed_time);
    writeln("  x1 = ", xx1);
    writeln("  exact = ", solution1(t1));

    writeln("Done.");
}
