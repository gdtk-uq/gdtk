// number.d
module number;
import std.format;
import std.math;
import complexify; 
 
version(complex_numbers) {
    alias number = Complex!double;
} else {
    alias number = double;
}

/*
  Returns true if both components are approximately equal.
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
    if (v1.length != v2.length) { return false; }
    foreach (i; 0 .. v1.length) {
        if (!approxEqual(v1[i].re, v2[i].re, maxRelDiff, maxAbsDiff)) { return false; }
    }
    return true;
}

/*
  Commonly used message in assert statements.
*/
string failedUnitTest(size_t lineNo = __LINE__,
                      string fileName = __FILE__)
{
    return format("Unit test failure on line %d in file %s\n", lineNo, fileName);
}

