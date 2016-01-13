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
import core.stdc.stdlib : exit;
import gas.gas_model;
import gas.physical_constants;
import gas.thermo.evt_eos;
import gas.thermo.cea_thermo_curves;
// For testing, use solve from either brent or ridder.
import nm.brent;
// import nm.ridder;
import nm.bracketing;
import util.lua;
import util.lua_service;

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
	    e = _curves[isp].eval_h(Q.T[0]) - _R[isp]*Q.T[0];
	}
	Q.e[0] = mass_average(Q, _energy);
    }
    override void update_temperature(ref GasState Q)
    {
	double Tsave = Q.T[0]; // Keep a copy for diagnostics purpose.
	// We'll adjust the temperature estimate until the energy is within TOL Joules.
	// Surely 1/100 of a Joule is sufficient precision when we are talking of megaJoules.
	double TOL = 1.0e-2;
	// The "target" energy is the value we will iterate to find.
	// We search (via a numerical method) for the temperature
	// value that gives this target energy.
	double e_tgt = Q.e[0];
	// delT is the initial guess for a bracket size.
	// We set this quite agressivley at 10 K hoping to
	// keep the number of iterations required to a small
	// value. In general, we are taking small timesteps
	// so the value of temperature from the previous step
	// should not change too much. If it does, there should
	// be enough robustness in the bracketing and
	// the function-solving method to handle this.
	double delT = 10.0;
	double T1 = fmax(Q.T[0] - delT/2, T_MIN);
	double T2 = T1 + delT;

	auto zeroFun = delegate (double T) {
	    Q.T[0] = T;
	    update_energy(Q);
	    return e_tgt - Q.e[0];
	};

	if ( bracket!zeroFun(T1, T2, T_MIN) == -1 ) {
	    string msg = "The 'bracket' function failed to find temperature values\n";
	    msg ~= "that bracketed the zero function in ThermallyPerfectGasMixEOS.eval_temperature().\n";
	    msg ~= format("The final values are: T1 = %12.6f and T2 = %12.6f\n", T1, T2);
	    msg ~= format("The initial temperature guess was: %12.6f\n", Tsave);
	    msg ~= format("The target energy value was: %12.6f\n", e_tgt);
	    msg ~= format("The GasState is currently:\n");
	    msg ~= Q.toString() ~ "\n";
	    throw new Exception(msg);
	}
	try {
	    Q.T[0] = solve!zeroFun(T1, T2, TOL);
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
	    throw new Exception(msg);
	}
    }
private:
    double[] _R;
    CEAThermo[] _curves;
    // Private working arrays
    double[] _energy;
}

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
	auto L = init_lua_State("sample-data/O2-N2-H2.lua");
	string[] species;
	getArrayOfStrings(L, LUA_GLOBALSINDEX, "species", species);
	ThermallyPerfectGasMixEOS tpgm = createThermallyPerfectGasMixEOS(species, L);
	auto Q = new GasState(3, 1);
	Q.massf[0] = 0.2; Q.massf[1] = 0.7; Q.massf[2] = 0.1;
	Q.T[0] = 1000.0;
	tpgm.update_energy(Q);
	assert(approxEqual(1031849.875, Q.e[0], 1.0e-6), failedUnitTest());
	// Now set T[0] a little off, say 1500.0.
	// Using Newton iterations, finding a temperature near the
	// CEA polynomial breaks was problematic. Brent's method
	// should do better.
	Q.T[0] = 1500.0;
	tpgm.update_temperature(Q);
	assert(approxEqual(1000.0, Q.T[0], 1.0e-6), failedUnitTest());

	return 0;
    }
}


