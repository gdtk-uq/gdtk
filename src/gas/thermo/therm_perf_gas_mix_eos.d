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
	    e = _curves[isp].eval_h(Q.Ttr) - _R[isp]*Q.Ttr;
	}
	Q.u = mass_average(Q, _energy);
    }
    override void update_temperature(ref GasState Q)
    {
	double Tsave = Q.Ttr; // Keep a copy for diagnostics purpose.
	// We'll adjust the temperature estimate until the energy is within TOL Joules.
	// Surely 1/100 of a Joule is sufficient precision when we are talking of megaJoules.
	double TOL = 1.0e-2;
	// The "target" energy is the value we will iterate to find.
	// We search (via a numerical method) for the temperature
	// value that gives this target energy.
	double e_tgt = Q.u;
	// delT is the initial guess for a bracket size.
	// We set this quite agressivley at 10 K hoping to
	// keep the number of iterations required to a small
	// value. In general, we are taking small timesteps
	// so the value of temperature from the previous step
	// should not change too much. If it does, there should
	// be enough robustness in the bracketing and
	// the function-solving method to handle this.
	double delT = 10.0;
	double T1 = fmax(Q.Ttr - delT/2, T_MIN);
	double T2 = T1 + delT;

	if ( bracket(T1, T2, e_tgt, Q, T_MIN) == -1 ) {
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
	    Q.Ttr = solve(T1, T2, TOL, e_tgt, Q);
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

    double zeroFun(double T, double e_tgt, ref GasState Q)
    // Helper function for update_temperature.
    {
	Q.Ttr = T;
	update_energy(Q);
	return e_tgt - Q.u;
    }

    int bracket(ref double x1, ref double x2, double e_tgt, ref GasState Q,
		double x1_min = -1.0e99, double x2_max = +1.0e99,
		int max_try=50, double factor=1.6)
    {
	if ( x1 == x2 ) {
	    throw new Exception("Bad initial range given to bracket.");
	}
	double f1 = zeroFun(x1, e_tgt, Q);
	double f2 = zeroFun(x2, e_tgt, Q);
	for ( int i = 0; i < max_try; ++i ) {
	    if ( f1*f2 < 0.0 ) return 0; // we have success
	    if ( abs(f1) < abs(f2) ) {
		x1 += factor * (x1 - x2);
		//prevent the bracket from being expanded beyond a specified domain
		x1 = fmax(x1_min, x1);
		f1 = zeroFun(x1, e_tgt, Q);
	    } else {
		x2 += factor * (x2 - x1);
		x2 = fmin(x2_max, x2);
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
    double solve(double x1, double x2, double tol, double e_tgt, ref GasState Q) 
    {
	const int ITMAX = 100;           // maximum allowed number of iterations
	const double EPS=double.epsilon; // Machine floating-point precision
	double a = x1;
	double b = x2;
	double fa = zeroFun(a, e_tgt, Q);
	double fb = zeroFun(b, e_tgt, Q);
	if ( abs(fa) == 0.0 ) return a;
	if ( abs(fb) == 0.0 ) return b;
	if ( fa * fb > 0.0 ) {
	    // Supplied bracket does not encompass zero of the function.
	    string msg = "Root must be bracketed by x1 and x2\n";
	    msg ~= format("x1=%g f(x1)=%g x2=%g f(x2)=%g\n", x1, fa, x2, fb); 
	    throw new Exception(msg);
	}
	double c = b;
	double fc = fb;
	for ( int iter=0; iter<ITMAX; iter++ ) {
	    double d, e, tol1, xm;
	    if ( (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) ) {
		c = a;
		fc = fa;
		e = d = b-a;
	    }
	    if ( abs(fc) < abs(fb) ) {
		a = b;
		b = c;
		c = a;
		fa = fb;
		fb = fc;
		fc = fa;
	    }
	    tol1 = 2.0*EPS*abs(b)+0.5*tol;  // Convergence check
	    xm = 0.5*(c-b);
	    if ( abs(xm) <= tol1 || fb == 0.0 ) return b;   // if converged let's return the best estimate
	    if ( abs(e) >= tol1 && abs(fa) > abs(fb) ) {
		double p, q, r, s;
		s = fb/fa;         // Attempt inverse quadratic interpolation for new bound
		if ( a == c ) {
		    p = 2.0*xm*s;
		    q = 1.0-s;
		}
		else {
		    q = fa/fc;
		    r = fb/fc;
		    p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
		    q = (q-1.0)*(r-1.0)*(s-1.0);
		}
		if ( p > 0.0 )
		    q = -q;    // Check whether quadratic interpolation is in bounds
		p = abs(p);
		double min1 = 3.0*xm*q-abs(tol1*q);
		double min2 = abs(e*q);
		if (2.0*p < (min1 < min2 ? min1 : min2) ) {
		    e = d;     // Accept interpolation
		    d = p/q;
		}
		else {
		    d = xm;  // else Interpolation failed, use bisection
		    e = d;
		}
	    }
	    else {
		d = xm;   // Bounds decreasing too slowly, use bisection
		e = d;
	    }
	    a = b;       // Move last guess to a
	    fa = fb;
	    if ( abs(d) > tol1 )  // Evaluate new trial root
		b += d;
	    else {
		b += copysign(tol1, xm);
	    }
	    fb = zeroFun(b, e_tgt, Q);	    
	}
	throw new Exception("Maximum number of iterations exceeded in solve (therm_perf_gas_mix)");
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
	Q.Ttr = 1000.0;
	tpgm.update_energy(Q);
	assert(approxEqual(1031849.875, Q.u, 1.0e-6), failedUnitTest());
	// Now set Ttr a little off, say 1500.0.
	// Using Newton iterations, finding a temperature near the
	// CEA polynomial breaks was problematic. Brent's method
	// should do better.
	Q.Ttr = 1500.0;
	tpgm.update_temperature(Q);
	assert(approxEqual(1000.0, Q.Ttr, 1.0e-6), failedUnitTest());

	return 0;
    }
}
