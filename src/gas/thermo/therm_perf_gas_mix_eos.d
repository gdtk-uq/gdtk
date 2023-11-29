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
import nm.nm_exception : NumericalMethodException;
import ntypes.complex;
import nm.number;
import nm.newton : solve;
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
    this(double[] R, CEAThermoCurve[] curves)
    {
        _R = R.dup;
        _curves = curves.dup;
        _vals.length = R.length;
    }
    override void update_energy(ref GasState Q)
    {
        foreach ( isp, ref e; _vals ) {
            e = _curves[isp].eval_h(Q.T) - _R[isp]*Q.T;
        }
        Q.u = mass_average(Q, _vals);
    }
    override void update_temperature(ref GasState Q)
    {
        number Tsave = Q.T; // Keep a copy for diagnostics purpose.
        double TOL = 1.0e-6;
        //
        // PJ 2020-01-01
        // Hack to cover over low energy problems with large ionization fractions.
        // Should not be needed if things are being solved well.
        //   number u_original = Q.u;
        //   Q.T = to!number(T_MIN+1.0);
        //   update_energy(Q);
        //   if (u_original > Q.u) { Q.u = u_original; }
        // End of hack.
        //
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
        number delT = 1000.0;
        version(complex_numbers) {
            number T1 = fmax(Q.T - 0.5*delT, to!number(T_MIN));
        } else {
            double T1 = std.math.fmax(Q.T - 0.5*delT, T_MIN);
        }
        number T2 = T1 + delT;

        number zeroFn(number T)
        {
            Q.T = T;
            update_energy(Q);
            /*
            debug {
                writefln("T= %.3f  u= %.3e  Fn= %.6e ", T, Q.u, e_tgt - Q.u);
            }
            */
            return e_tgt - Q.u;
        }

        number dzdT(number T)
        {
            // We evaluate Cv.
            foreach ( isp, ref Cv; _vals ) {
                Cv = _curves[isp].eval_Cp(T) - _R[isp];
            }
            return -1.0*mass_average(Q, _vals);
        }

        try {
            Q.T = solve!(zeroFn, dzdT)(Tsave, T1, T2, TOL);
        }
        catch (NumericalMethodException e) {
            // If we fail, we'll have a second attempt with an extended temperature range.
            try {
                Q.T = solve!(zeroFn, dzdT)(Tsave, to!number(T_MIN), to!number(T_MAX), TOL);
            }
            catch (NumericalMethodException e) {
                string msg = "There was a problem iterating to find temperature";
                debug {
                    msg ~= "\nin function ThermallyPerfectGasMix.update_temperature().\n";
                    msg ~= format("The initial temperature guess was: %12.6f\n", Tsave);
                    msg ~= format("The target energy value was: %12.6f\n", e_tgt);
                    msg ~= format("zeroFn(Tsave)=%g\n", zeroFn(Tsave));
                    msg ~= format("zeroFn(T_MIN)=%g\n", zeroFn(to!number(T_MIN)));
                    msg ~= format("zeroFn(T_MAX)=%g\n", zeroFn(to!number(T_MAX)));
                    msg ~= format("The GasState is currently:\n");
                    msg ~= Q.toString() ~ "\n";
                    msg ~= "The message from the solve function is:\n";
                    msg ~= e.msg;
                }
                // If we fail, leave temperature at value upon entry to method.
                Q.T = Tsave;
                update_energy(Q);
                throw new GasModelException(msg);
            }
        }
    } // end update_temperature()
private:
    double[] _R;
    CEAThermoCurve[] _curves;
    // Private working array
    number[] _vals;

} // end class ThermallyPerfectGasMixEOS


ThermallyPerfectGasMixEOS createThermallyPerfectGasMixEOS(string[] species, lua_State* L)
{
    double[] R = new double[species.length];
    CEAThermoCurve[] curves;
    foreach ( isp, s; species ) {
        lua_getglobal(L, s.toStringz);
        double M = getDouble(L, -1, "M");
        R[isp] = R_universal/M;
        lua_getfield(L, -1, "cea_thermo");
        curves ~= new CEAThermoCurve(L, R[isp]);
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
        getArrayOfStrings(L, "species", species);
        ThermallyPerfectGasMixEOS tpgm = createThermallyPerfectGasMixEOS(species, L);
        auto Q = GasState(3, 1);
        Q.massf[0] = 0.2; Q.massf[1] = 0.7; Q.massf[2] = 0.1;
        Q.T = 1000.0;
        tpgm.update_energy(Q);
        assert(isClose(1031849.875, Q.u, 1.0e-6), failedUnitTest());
        // Now set T a little off, say 1500.0.
        Q.T = 1500.0;
        tpgm.update_temperature(Q);
        assert(isClose(1000.0, Q.T, 1.0e-6), failedUnitTest());

        return 0;
    }
}
