/**
 * Authors: Rowan G. and Peter J.
 * Date: 2017-07-13
 *
 * A model for for vibrational relaxation kinetics in
 * single-species diatomic nitrogen.
 */

module kinetics.two_temperature_nitrogen_kinetics;

import std.math : exp, pow;

import gas;
import kinetics.thermochemical_reactor;

final class VibRelaxNitrogen : ThermochemicalReactor {
    this(string fname, GasModel gmodel)
    {
	super(gmodel);
	_Q_eq = new GasState(gmodel);
    }

    override void opCall(GasState Q, double tInterval, ref double dtSuggest, ref double[] params)
    {
	// Hard-code Blackman relaxation time for present.
	double pAtm = Q.p/P_atm;
	double tau = (7.12e-9/pAtm)*exp(124.07/pow(Q.T, 1./3.));
	// Find the total internal energy in the gas
	double uTotal = Q.u + Q.u_modes[0];
	// Find the vib energy at equilibrium with T
	_Q_eq.T = Q.T;
	_Q_eq.T_modes[0] = Q.T;
	_Q_eq.p = Q.p;
	_gmodel.update_thermo_from_pT(_Q_eq);
	double u_v_eq = _Q_eq.u_modes[0];
	double u_v = Q.u_modes[0];
	Q.u_modes[0] = u_v_eq + (u_v - u_v_eq)*exp(-tInterval/tau);
	Q.u = uTotal - Q.u_modes[0];
	_gmodel.update_thermo_from_rhou(Q);
	_gmodel.update_sound_speed(Q);
    }

private:
    GasState _Q_eq;
}
