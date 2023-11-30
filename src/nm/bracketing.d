/**
 * Bracket a root of f(x).
 *
 * Params:
 *    f: user-supplied function f(x)
 *    x1: first end of range
 *    x2: other end of range
 *  x1min: limits the value of x1 so that it does not become lower than a certain value
 *
 * Returns:
 *    0: successfully bracketed a root
 *   -1: failed to bracket a root
 *
 * On return (x1, x2) should bracketing the root.
 */

module nm.bracketing;
import std.math;
import std.algorithm;
import ntypes.complex;
debug { import std.stdio; }

int bracket(alias f, T)(ref T x1, ref T x2,
                        T x1_min = -1.0e38, T x2_max = +1.0e38,
                        int max_try=50, T factor=1.6)
    if (is(typeof(f(0.0)) == double) ||
        is(typeof(f(0.0)) == float)  ||
        is(typeof(f(Complex!double(0.0))) == Complex!double))
{
    if (x1 == x2) {
        throw new Exception("Bad initial range given to bracket.");
    }
    // We assume that x1 < x2.
    if (x1 > x2) { swap(x1, x2); }
    if (x1_min > x2_max) { swap(x1_min, x2_max); }
    T f1 = f(x1);
    T f2 = f(x2);
    for (int i = 0; i < max_try; ++i) {
        if (f1*f2 < 0.0) { return 0; } // we have success
        if (fabs(f1) < fabs(f2)) {
            x1 += factor * (x1 - x2);
            //prevent the bracket from being expanded beyond a specified domain
            x1 = fmax(x1_min, x1);
            try {
                f1 = f(x1);
            } catch (Exception e) {
                // Presume that we have gone into an invalid region,
                // so reset to the boundary.
                // debug { writeln("bracket trying x1 e.msg=", e.msg); }
                x1 = x1_min;
                f1 = f(x1);
            }
        } else {
            x2 += factor * (x2 - x1);
            x2 = fmin(x2_max, x2);
            try {
                f2 = f(x2);
            } catch (Exception e) {
                // Presume that we have gone into an invalid region,
                // so reset to the boundary.
                // debug { writeln("bracket trying x2 e.msg=", e.msg); }
                x2 = x2_max;
                f2 = f(x2);
            }
        }
    }
    // If we leave the loop here, we were unsuccessful.
    return -1;
} // end bracket()

version(bracketing_test) {
    import util.msg_service;
    import std.conv;
    import nm.number;
    int main() {
        number test_fun_2(number x, number a) {
            return a*x + sin(x) - exp(x);
        }
        number my_a = 3.0;
        auto test_fun_3 = delegate (number x) { return test_fun_2(x, my_a); };
        number x1 = 0.4;
        number x2 = 0.5;
        assert(bracket!(test_fun_3,number)(x1, x2) == to!number(0), failedUnitTest());

        return 0;
    }
}
