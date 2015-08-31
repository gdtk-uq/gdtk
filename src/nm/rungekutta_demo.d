/** rungekutta_demo.d
 *
 * Try out the Runge-Kutta ODE stepper.
 *
 * Author: Peter J.
 * Version: 2014-06-15, fresh start with Dlang
 */

import std.stdio;
import std.math;
import rungekutta;

/** Test system 1 
 * Third-order system with a simple analytic solution.
 * Adapted from section 11.3 in Cheney and Kincaid, 6th ed.
 * Except for the zero-based indexing, the notation is
 * chosen to match that in the text.
 */
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

void main()
{
    writeln("Begin demonstration of ODE stepper...");
    double[] x0=[0.0, 0.0, 0.0];
    double[] x1=x0.dup;
    double[] err=x0.dup;
    auto work = allocate_rk45_workspace(3);
    double t1 = rkf45_step!(testSystem1)(0.0, 0.2, x0, x1, err, work);
    writeln("x1 = ", x1);
    writeln("exact = ", solution1(t1));
    writeln("Done.");
}
