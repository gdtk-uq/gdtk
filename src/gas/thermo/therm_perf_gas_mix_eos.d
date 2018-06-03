/**
 * therm_perf_gas_mix_eos.d 
 * This implements a caloric equation of state
 * relating temperature and internal energy
 * for a mixture of thermally perfect gases.
 * In this model, the specific heat at constant
 * pressure (Cp) is a function of temperature.
 * This gas model is commonly used in reacting gas
 * flow calculations at moderate-to-high temperatures.
 *
 * The energy is a function of temperature and composition
 * only in this model. The energy has no dependence on density
 * or pressure. In thermodynamics jargon, the energy is
 * evaluated at the zero-pressure limit.
 *
 * Author: Rowan G. and Peter J.
 */

module gas.thermo.therm_perf_gas_mix_eos;

import std.math;
import std.stdio;
import std.string;
import std.conv;
import core.stdc.stdlib : exit;
import nm.complex;
import nm.number;
import util.lua;
import util.lua_service;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;
import gas.thermo.evt_eos;
import gas.thermo.cea_thermo_curves;

// These bracket limits are conservative limits that
// are used when we fail to bracket the temperature using
// a nearby guess based on the temperature at the previous
// timestep.
immutable double T_bracket_low = 20.0; // K
immutable double T_bracket_high = 100000.0; // K

/++
  ThermallyPerfectGasMixEOS is a caloric equation of state.

  In this model, Cp is a function of temperature. Consequently,
  the other thermodynamic properties are also functions of temperature.
  The "perfect" nature of the gas means that the energy has no dependence
  on pressure or density.
+/
class ThermallyPerfectGasMixEOS : EVT_EOS {
public:
    this(double[] R, CEAThermo[] curves)
    {
        _R = R.dup;
        _curves = curves.dup;
        _energy.length = R.length;
        
    }
    override void update_energy(ref GasState Q)
    {
        foreach ( isp, ref e; _energy ) {
            e = _curves[isp].eval_h(Q.T) - _R[isp]*Q.T;
        }
        Q.u = mass_average(Q, _energy);
    }
    override void update_temperature(ref GasState Q)
    {
        number Tsave = Q.T; // Keep a copy for diagnostics purpose.
        // We'll adjust the temperature estimate until the energy is within TOL Joules.
        // Surely 1/1000 of a Joule is sufficient precision when we are talking of megaJoules.
        double TOL = 1.0e-3;
        // The "target" energy is the value we will iterate to find.
        // We search (via a numerical method) for the temperature
        // value that gives this target energy.
        number e_tgt = Q.u;
        // delT is the initial guess for a bracket size.
        // We set this quite agressivley at 100 K hoping to
        // keep the number of iterations required to a small
        // value. In general, we are taking small timesteps
        // so the value of temperature from the previous step
        // should not change too much. If it does, there should
        // be enough robustness in the bracketing and
        // the function-solving method to handle this.
        number delT = 100.0;
        version(complex_numbers) {
            number T1 = nm.complex.fmax(Q.T - 0.5*delT, to!number(T_MIN));
        } else {
            double T1 = std.math.fmax(Q.T - 0.5*delT, T_MIN);
        }
        number T2 = T1 + delT;

        if (bracket(T1, T2, e_tgt, Q, to!number(T_MIN)) == -1) {
            // We have a fall back if our aggressive search failed.
            // We apply a very conservative range:
            T1 = T_bracket_low;
            T2 = T_bracket_high;
            if (bracket(T1, T2, e_tgt, Q, to!number(T_MIN)) == -1) {
                string msg = "The 'bracket' function failed to find temperature values\n";
                msg ~= "that bracketed the zero function in ThermallyPerfectGasMixEOS.eval_temperature().\n";
                msg ~= format("The final values are: T1 = %12.6f and T2 = %12.6f\n", T1, T2);
                msg ~= format("The initial temperature guess was: %12.6f\n", Tsave);
                msg ~= format("The target energy value was: %12.6f\n", e_tgt);
                msg ~= format("The bracket limits were: low=%12.6f  high=%12.6f\n", T_bracket_low, T_bracket_high);
                Q.T = 20.0;
                update_energy(Q);
                msg ~= format("Energy at 20.0: %12.6f\n", Q.u);
                Q.T = 100000.0;
                update_energy(Q);
                msg ~= format("Energy at 100000.0: %12.6f\n", Q.u);
                msg ~= format("The GasState is currently:\n");
                msg ~= Q.toString() ~ "\n";
                throw new GasModelException(msg);
            }
        }
        try {
            Q.T = solve(T1, T2, TOL, e_tgt, Q);
        }
        catch ( Exception e ) {
            string msg = "There was a problem iterating to find temperature\n";
            msg ~= "in function ThermallyPerfectGasMix.eval_temperature().\n";
            msg ~= format("The initial temperature guess was: %12.6f\n", Tsave);
            msg ~= format("The target energy value was: %12.6f\n", e_tgt);
            msg ~= format("The GasState is currently:\n");
            msg ~= Q.toString() ~ "\n";
            msg ~= "The message from the solve function is:\n";
            msg ~= e.msg;
            throw new GasModelException(msg);
        }
    } // end update_temperature()
private:
    double[] _R;
    CEAThermo[] _curves;
    // Private working arrays
    number[] _energy;

    //--------------------------------------------------------------------------------
    // Bracketing and solve functions copied in from nm/bracketing.d and nm/brent.d
    // and explicitly specialized for the temperature update.
    //
    // Originally, we would call the solving functions, passing to them a delegate f,
    // and let them determine T such that f(T)=0, however, creating a delegate involves
    // memory allocation and we suspect that it hinders parallel calculations.
    //
    // PJ, 10-Sep-2016

    import std.math;
    import std.algorithm;

    number zeroFun(number T, number e_tgt, ref GasState Q)
    // Helper function for update_temperature.
    {
        Q.T = T;
        update_energy(Q);
        return e_tgt - Q.u;
    }

    int bracket(ref number x1, ref number x2, number e_tgt, ref GasState Q,
                number x1_min = -1.0e99, number x2_max = +1.0e99,
                int max_try=50, double factor=1.6)
    {
        if ( x1 == x2 ) {
            throw new GasModelException("Bad initial range given to bracket.");
        }
        number f1 = zeroFun(x1, e_tgt, Q);
        number f2 = zeroFun(x2, e_tgt, Q);
        for ( int i = 0; i < max_try; ++i ) {
            if ( f1*f2 < 0.0 ) return 0; // we have success
            if ( abs(f1) < abs(f2) ) {
                x1 += factor * (x1 - x2);
                //prevent the bracket from being expanded beyond a specified domain
                version(complex_numbers) {
                    x1 = nm.complex.fmax(x1_min, x1);
                } else {
                    x1 = std.math.fmax(x1_min, x1);
                }
                f1 = zeroFun(x1, e_tgt, Q);
            } else {
                x2 += factor * (x2 - x1);
                version(complex_numbers) {
                    x2 = nm.complex.fmin(x2_max, x2);
                } else {
                    x2 = std.math.fmin(x2_max, x2);
                }
                f2 = zeroFun(x2, e_tgt, Q);
            }
        }
        // If we leave the loop here, we were unsuccessful.
        return -1;
    } // end bracket()

    /**
     * Locate a root of f(x) known to lie between x1 and x2. The method
     * is guaranteed to converge (according to Brent) given the initial 
     * x values bound the solution. The method uses root bracketing,
     * bisection and inverse quadratic interpolation.
     *
     * Params:
     *    f: user-supplied function f(x)
     *    x1: first end of range
     *    x2: other end of range
     *    tol: minimum size for range
     *
     * Returns:
     *    b, a point near the root.
     */
    number solve(number x1, number x2, double tol, number e_tgt, ref GasState Q) 
    {
        const int ITMAX = 100;           // Maximum allowed number of iterations
        const double EPS=double.epsilon; // Machine floating-point precision
        number a = x1;
        number b = x2;
        number fa = zeroFun(a, e_tgt, Q);
        number fb = zeroFun(b, e_tgt, Q);
        if (abs(fa) == 0.0) { return a; }
        if (abs(fb) == 0.0) { return b; }
        if (fa * fb > 0.0) {
            // Supplied bracket does not encompass zero of the function.
            string msg = "Root must be bracketed by x1 and x2\n";
            msg ~= format("x1=%g f(x1)=%g x2=%g f(x2)=%g\n", x1, fa, x2, fb); 
            throw new GasModelException(msg);
        }
        number c = b;
        number fc = fb;
        // Define d, e outside the loop body so that
        // they don't get initialized to nan each pass.
        number d, e;
        foreach (iter; 0 .. ITMAX) {
            if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
                // On first pass fc==fb and we expect to enter here.
                c = a;
                fc = fa;
                e = d = b-a;
            }
            if (abs(fc) < abs(fb)) {
                a = b;
                b = c;
                c = a;
                fa = fb;
                fb = fc;
                fc = fa;
            }
            // Convergence check
            number tol1 = 2.0*EPS*abs(b)+0.5*tol;
            number xm = 0.5*(c-b);
            if (abs(xm) <= tol1 || fb == 0.0) {
                // If converged, let's return the best estimate
                return b;
            }
            if (abs(e) >= tol1 && abs(fa) > abs(fb)) {
                // Attempt inverse quadratic interpolation for new bound
                number p, q;
                number s = fb/fa;
                if (a == c) {
                    p = 2.0*xm*s;
                    q = 1.0-s;
                } else {
                    q = fa/fc;
                    number r = fb/fc;
                    p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                    q = (q-1.0)*(r-1.0)*(s-1.0);
                }
                // Check whether quadratic interpolation is in bounds
                if (p > 0.0) { q = -q; }
                p = abs(p);
                number min1 = 3.0*xm*q-abs(tol1*q);
                number min2 = abs(e*q);
                if (2.0*p < (min1 < min2 ? min1 : min2)) {
                    // Accept interpolation
                    e = d;
                    d = p/q;
                } else {
                    // Interpolation failed, use bisection
                    d = xm;
                    e = d;
                }
            } else {
                // Bounds decreasing too slowly, use bisection
                d = xm;
                e = d;
            }
            // Move last guess to a
            a = b;
            fa = fb;
            // Evaluate new trial root
            if (abs(d) > tol1) {
                b += d;
            } else {
                version(complex_numbers) {
                    b += nm.complex.copysign(tol1, xm);
                } else {
                    b += std.math.copysign(tol1, xm);
                }
            }
            fb = zeroFun(b, e_tgt, Q);      
        }
        throw new GasModelException("Maximum number of iterations exceeded in solve.");
    } // end solve()
} // end class ThermallyPerfectGasMixEOS


ThermallyPerfectGasMixEOS createThermallyPerfectGasMixEOS(string[] species, lua_State* L)
{
    double[] R = new double[species.length];
    CEAThermo[] curves;
    foreach ( isp, s; species ) {
        lua_getfield(L, LUA_GLOBALSINDEX, s.toStringz);
        double M = getDouble(L, -1, "M");
        R[isp] = R_universal/M;
        lua_getfield(L, -1, "cea_thermo");
        curves ~= createCEAThermo(L, R[isp]);
        lua_pop(L, 1);
        lua_pop(L, 1);
    }
    return new ThermallyPerfectGasMixEOS(R, curves);
}

version(therm_perf_gas_mix_eos_test) {
    import std.math;
    import util.msg_service;
    int main() {
        auto L = init_lua_State();
        doLuaFile(L, "sample-data/O2-N2-H2.lua");
        string[] species;
        getArrayOfStrings(L, LUA_GLOBALSINDEX, "species", species);
        ThermallyPerfectGasMixEOS tpgm = createThermallyPerfectGasMixEOS(species, L);
        auto Q = new GasState(3, 1);
        Q.massf[0] = 0.2; Q.massf[1] = 0.7; Q.massf[2] = 0.1;
        Q.T = 1000.0;
        tpgm.update_energy(Q);
        assert(approxEqual(1031849.875, Q.u, 1.0e-6), failedUnitTest());
        // Now set T a little off, say 1500.0.
        // Using Newton iterations, finding a temperature near the
        // CEA polynomial breaks was problematic. Brent's method
        // should do better.
        Q.T = 1500.0;
        tpgm.update_temperature(Q);
        assert(approxEqual(1000.0, Q.T, 1.0e-6), failedUnitTest());

        return 0;
    }
}
