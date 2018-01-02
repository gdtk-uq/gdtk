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
import std.stdio; // for debugging writes

int bracket(alias f)(ref double x1, ref double x2,
                     double x1_min = -1.0e99, double x2_max = +1.0e99,
                     int max_try=50, double factor=1.6)
    if (is(typeof(f(0.0)) == double) || is(typeof(f(0.0)) == float))
{
    if (x1 == x2) {
        throw new Exception("Bad initial range given to bracket.");
    }
    // We assume that x1 < x2.
    if (x1 > x2) { swap(x1, x2); }
    if (x1_min > x2_max) { swap(x1_min, x2_max); }
    double f1 = f(x1);
    double f2 = f(x2);
    for (int i = 0; i < max_try; ++i) {
        if (f1*f2 < 0.0) return 0; // we have success
        if (abs(f1) < abs(f2)) {
            x1 += factor * (x1 - x2);
            //prevent the bracket from being expanded beyond a specified domain
            x1 = fmax(x1_min, x1);
            f1 = f(x1);
        } else {
            x2 += factor * (x2 - x1);
            x2 = fmin(x2_max, x2);
            f2 = f(x2);
        }
    }
    // If we leave the loop here, we were unsuccessful.
    return -1;
} // end bracket()

version(bracketing_test) {
    import util.msg_service;
    int main() {
        double test_fun_2(double x, double a) {
            return a*x + sin(x) - exp(x);
        }
        double my_a = 3.0;
        auto test_fun_3 = delegate (double x) { return test_fun_2(x, my_a); };
        double x1 = 0.4;
        double x2 = 0.5;
        assert(bracket!test_fun_3(x1, x2) == 0, failedUnitTest());

        return 0;
    }
}
