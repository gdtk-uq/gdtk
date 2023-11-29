/**
3 * Solve a nonlinear equation f(x)=0 using Newton's method.
 *
 * Caller needs to supply f(x) and f'(x).
 *
 * Authors: Rowan G. and Peter J.
 * Version: 2019-12-05, first cut
 */

module nm.newton;

import std.math;
import ntypes.complex;
import std.format;

/**
 * Locate a root of f(x) using Newton's method.
 *
 * This algorithm is a variation on the classic Newton
 * method in that it has some safeguards if the method
 * wants to wander into territory beyond some bounds.
 * Here we denote those bounds as xMin and xMax.
 *
 * This is an implementation of the function 'rtsafe'
 * that appears in Section 9.4 of Press et al.
 * Note I have added a modification such that the user
 * supplies an initial guess, whereas the original
 * suggestion of Press et al. is to use:
 *    guess = 0.5*(xMin + xMax);
 *
 * Reference:
 * Press, Teukolsky, Vetterling, Flannery (1992)
 * Numerical Recipes in C, The Art of Scientific Computing, 2nd ed.
 * Cambridge University Press
 */

import nm.nm_exception : NumericalMethodException;


version(complex_numbers) {
// For the complex numbers version of the code we need
// a Newton's method with a fixed number of iterations.
// An explanation can be found in:
//     Efficient Construction of Discrete Adjoint Operators on Unstructured Grids
//     by Using Complex Variables, pg. 10, Nielsen et al., AIAA Journal, 2006.
// TODO: add in safeguards from the real number version.
@nogc
T solve(alias f, alias dfdx, T)(T x0, T xMin, T xMax, double tol=1.0e-9)
if ( (is(typeof(f(0.0)) == double) && is(typeof(dfdx(0.0)) == double)) ||
     (is(typeof(f(0.0)) == float) && is(typeof(dfdx(0.0)) == float)) ||
     (is(typeof(f(Complex!double(0.0))) == Complex!double) &&
      is(typeof(dfdx(Complex!double(0.0))) == Complex!double)) )
    {
        const int MAXIT = 10; // maximum number of iterations                
        T rts = x0;
        T dx = 0.0;
        T f0 = f(rts);
        T df0 = dfdx(rts);
        
        foreach (j; 0 .. MAXIT) {
            dx = f0/df0;
            rts -= dx;
            f0 = f(rts);
            df0 = dfdx(rts);
        }
        
        return rts;
    }
} else {
@nogc
T solve(alias f, alias dfdx, T)(T x0, T xMin, T xMax, double tol=1.0e-9)
if ( (is(typeof(f(0.0)) == double) && is(typeof(dfdx(0.0)) == double)) ||
     (is(typeof(f(0.0)) == float) && is(typeof(dfdx(0.0)) == float)) ||
     (is(typeof(f(Complex!double(0.0))) == Complex!double) &&
      is(typeof(dfdx(Complex!double(0.0))) == Complex!double)) )
    {
        const int MAXIT = 30; // maximum number of iterations
        T xL = xMin;
        T xH = xMax;
        T fL = f(xL);
        T fH = f(xH);
        if ( (fL > 0.0 && fH > 0.0) || (fL < 0.0 && fH < 0.0) ) {
            string msg = "Root must be bracketed in call to nm.newton.solve";
            throw new NumericalMethodException(msg);
        }
        
        if (fL == 0.0) return xMin;
        if (fH == 0.0) return xMax;
        
        // Orient the search so that that f(xL) < 0.0
        if (fL < 0.0) {
            xL = xMin;
            xH = xMax;
        }
        else {
            xH = xMin;
            xL = xMax;
        }
        
        T rts = x0;
        T dxold = (xMax - xMin);
        T dx = dxold;
        T f0 = f(rts);
        T df0 = dfdx(rts);
        foreach (j; 0 .. MAXIT) {
            // Bisect if Newton prediction is out of range.
            if ( (((rts-xH)*df0 - f0)*((rts-xL)*df0 - f0) > 0.0) ||
                 (fabs(2.0*f0) > fabs(dxold*df0)) ) {
                dxold = dx;
                dx = 0.5*(xH - xL);
                rts = xL + dx;
                // Check if change in root is negligible.
                // Accept it, if it is.
                if (xL == rts) return rts; 
            }
            else { // A Newton step is acceptable
                dxold = dx;
                dx = f0/df0;
                T tmp = rts;
                rts -= dx;
                if (tmp == rts) return rts;
            }
            if (fabs(dx) < tol) return rts;
            // Otherwise, re-evaluate for next iteration
            f0 = f(rts);
            df0 = dfdx(rts);
            if (f0 < 0.0) {
                xL = rts;
            }
            else {
                xH = rts;
            }
        }
        // When successful, we should never reach here.
        string msg = "Newton method failed to converge";
        debug { msg ~= format(" in %d iterations.", MAXIT); }
        throw new NumericalMethodException(msg);
    }
 }

version(newton_test) {
    import std.conv;
    import util.msg_service;
    import nm.number;
    int main() {
        @nogc number test_fun_1(number x) {
            return pow(x,3) + pow(x,2) - 3*x - 3;
        }
        @nogc number test_dfun_1(number x) {
            return 3.0*pow(x,2) + 2.0*x - 3.0;
        }
        assert(fabs(solve!(test_fun_1, test_dfun_1, number)(to!number(1.0), to!number(-5.0), to!number(5.0)) - 1.732051) < 1.0e-5,
               failedUnitTest());


        @nogc number test_fun_2(number x, number a) {
            return a*x + sin(x) - exp(x);
        }
        @nogc number test_dfun_2(number x, number a) {
            return a + cos(x) - exp(x);
        }
        number my_a = 3.0;
        @nogc number test_fun_3 (number x) {
            return test_fun_2(x, my_a);
        }
        @nogc number test_dfun_3 (number x) {
            return test_dfun_2(x, my_a);
        }
        assert(fabs(solve!(test_fun_3, test_dfun_3, number)(to!number(0.2), to!number(0.0), to!number(1.0)) - 0.3604217) < 1.0e-5,
               failedUnitTest());
        return 0;
    }
}
