/**
 number.d
 Specifies the type of numbers to use in the numerical methods.
*/

module nm.number;

import std.math;
import ntypes.complex;

version(complex_numbers) {
    immutable alias number = Complex!double;
} else {
    immutable alias number = double;
}

/**
 * Returns true if both components are approximately equal.
 */
@nogc
bool approxEqualNumbers(in Complex!double v1, in Complex!double v2,
                        double maxRelDiff=1.0e-2, double maxAbsDiff=1.0e-5)
{
    return (isClose(v1.re, v2.re, maxRelDiff, maxAbsDiff) &&
            isClose(v1.im, v2.im, maxRelDiff, maxAbsDiff));
}
@nogc
bool approxEqualNumbers(in Complex!double[] v1, in Complex!double[] v2,
                        double maxRelDiff=1.0e-2, double maxAbsDiff=1.0e-5)
{
    if (v1.length != v2.length) return false;
    foreach (i; 0 .. v1.length) {
        if (!isClose(v1[i].re, v2[i].re, maxRelDiff, maxAbsDiff) ||
            !isClose(v1[i].im, v2[i].im, maxRelDiff, maxAbsDiff)) return false;
    }
    return true;
}

@nogc
bool approxEqualNumbers(in double v1, in double v2,
                        double maxRelDiff=1.0e-2, double maxAbsDiff=1.0e-5)
{
    return isClose(v1, v2, maxRelDiff, maxAbsDiff);
}
@nogc
bool approxEqualNumbers(in double[] v1, in double[] v2,
                        double maxRelDiff=1.0e-2, double maxAbsDiff=1.0e-5)
{
    if (v1.length != v2.length) return false;
    foreach (i; 0 .. v1.length) {
        if (!isClose(v1[i], v2[i], maxRelDiff, maxAbsDiff)) return false;
    }
    return true;
}

