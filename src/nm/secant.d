/**
 * secant.d
 *
 * Solve a nonlinear equation f(x)=0 using the secant method.
 *
 * Peter J.
 * 2020-09-27: D version adapted from Python3 version
 */
module nm.secant;
import std.math;
import std.algorithm;
import ntypes.complex;

/**
 * The iterative secant method for zero-finding in one-dimension.
 *
 * Params:
 *   f: user-supplied function f(x)
 *   x0: first guess
 *   x1: second guess, presumably close to x0
 *   tol: stopping tolerance for f(x)=0
 *   max_iterations: to stop the iterations running forever, just in case...
 *
 * Returns:
 *   x such that f(x)=0
 */
T solve(alias f, T)(T x0, T x1, double tol=1.0e-11,
                    double limit0=0.0, double limit1=0.0,
                    int max_iterations=1000)
    if ( is(typeof(f(0.0)) == double) ||
         is(typeof(f(0.0)) == float)  ||
         is(typeof(f(Complex!double(0.0))) == Complex!double))
{
    // We're going to arrange x0 as the oldest (furtherest) point
    // and x1 and the closer-to-the-solution point.
    // x2, when we compute it, will be the newest sample point.
    T f0 = f(x0);
    T f1 = f(x1);
    if (abs(f0) < abs(f1)) { swap(x0, x1); swap(f0, f1); }
    bool have_limits = !((limit0 == 0.0) && (limit1 == 0.0));
    foreach (i; 0..max_iterations) {
        T df = f0 - f1;
        if (df.re == 0.0) { throw new Exception("Cannot proceed with zero slope."); }
        T x2 = x1 - f1 * (x0 - x1) / df;
        if (have_limits) {
            x2.re = max(limit0, x2.re);
            x2.re = min(limit1, x2.re);
        }
        T f2 = f(x2);
        x0 = x1; f0 = f1; x1 = x2; f1 = f2;
        if (abs(f2) < tol) { return x2; }
    }
    throw new Exception("Did not converge after max iterations.");
} // end solve()

version(secant_test) {
    import util.msg_service;
    import std.conv;
    import nm.number;
    int main() {
        number test_fun_1(number x) {
            return pow(x,3) + pow(x,2) - 3*x - 3;
        }
        number test_fun_2(number x, number a) {
            return a*x + sin(x) - exp(x);
        }
        assert(fabs(solve!(test_fun_1,number)(1.5, 1) - 1.732051) < 1.0e-5, failedUnitTest());
        number my_a = 3.0;
        auto test_fun_3 = delegate (number x) { return test_fun_2(x, my_a); };
        assert(fabs(solve!(test_fun_3,number)(0, 0.1) - 0.3604217) < 1.0e-5, failedUnitTest());

        return 0;
    }
}
