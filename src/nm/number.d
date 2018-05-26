/**
 number.d
 Specifices what type of numbers to use in the numerical methods.
*/

module nm.number;

import std.math;
import nm.complex; 
 
version(complex_numbers) {
    alias number = Complex!double;
} else {
    alias number = double;
}

/**
 * Returns true if both components are approximately equal.
 */
@nogc
bool approxEqualNumbers(in number v1, in number v2,
                        double maxRelDiff=1.0e-2, double maxAbsDiff=1.0e-5)
{
    return (approxEqual(v1.re, v2.re, maxRelDiff, maxAbsDiff) && 
            approxEqual(v1.im, v2.im, maxRelDiff, maxAbsDiff));
}
@nogc
bool approxEqualNumbers(in number[] v1, in number[] v2,
                        double maxRelDiff=1.0e-2, double maxAbsDiff=1.0e-5)
{
    if (v1.length != v2.length) return false;
    foreach (i; 0 .. v1.length) {
        if (!approxEqual(v1[i].re, v2[i].re, maxRelDiff, maxAbsDiff)) return false;
    }
    return true;
}

