// power_benchmark.d
// Compare the standard function and a truncated expansion.
// PJ, 2021-10-28
//
/+
peterj@prandtl ~/work/play/dlang/benchmark $ dmd power_benchmark.d 
peterj@prandtl ~/work/play/dlang/benchmark $ ./power_benchmark 
math.pow: 905 msecs
qd_power: 253 msecs
peterj@prandtl ~/work/play/dlang/benchmark $ dmd -O -release -noboundscheck power_benchmark.d
peterj@prandtl ~/work/play/dlang/benchmark $ ./power_benchmark 
math.pow: 848 msecs
qd_power: 275 msecs
peterj@prandtl ~/work/play/dlang/benchmark $ ldmd2 -O -release -noboundscheck power_benchmark.d
peterj@prandtl ~/work/play/dlang/benchmark $ ./power_benchmark
math.pow: 756 msecs
qd_power: 113 msecs
+/

import std.datetime.stopwatch : benchmark;
import std.math : pow;
import std.stdio : writefln;

void main()
{
    immutable N = 100_000;
    auto a = new double[N];
    auto b = new double[N];
    auto c = new double[N];
    auto d = new double[N];
    foreach (i; 0 .. N) { a[i] = 1.1;  b[i] = 1.3; }
    auto bm = benchmark!({
       foreach (i; 0 .. N) c[i] = pow(a[i], b[i]);
    }, {
       foreach (i; 0 .. N) d[i] = qd_power(a[i], b[i]);
    })(100); // number of executions of each tested function
    //
    writefln("math.pow: %s msecs", bm[0].total!"msecs");
    writefln("qd_power: %s msecs", bm[1].total!"msecs");
}

// Adapted from the qd_power macro in cfcfd3/lib/nm/source/

/** \file qd_power.h
 * \ingroup nm
 * \brief Quick_and_dirty power MACRO for real arguments.
 *
 * Compute z = x**y approximately for moderate values of y and
 * values of x close to 1.
 * The routine produces results within 1% relative error for
 * the following ranges of x (base) and y (exponent).
 * y = 0.2:  0.12 <= x <= 8.4
 * y = 4.0:  0.75 <= x <= 1.55
 * y = -3.0: 0.55 <= x <= 1.48
 *
 * \param x : base
 * \param y : exponent
 * \param z : an approximation to x**y for x near 1 and y small.
 *
 * \version 1.0     17-Mar-91
 *
 * \author PA Jacobs, ICASE
 *
 */

double qd_power(double x, double y)
{
    // Compute logarithm.
    double r = (x - 1.0) / (x + 1.0);
    double temp = r * r * r;
    double lnx = r + 0.333333333 * temp;
    temp *= r * r;
    lnx += 0.2 * temp;
    temp *= r * r;
    lnx += 0.142857143 * temp;  // This term = r**7 / 7
    lnx *= 2.0;
    // Scale by the exponent.
    double t = lnx * y;
    // Compute the exponential.
    temp = 0.5 * t * t;
    double z = 1 + t + temp;
    temp *= 0.333333333 * t;
    z += temp;
    temp *= 0.25 * t;  // This term = t**4/4!
    z += temp;
    temp *= 0.2 * t;
    z += temp;
    return z;
}
