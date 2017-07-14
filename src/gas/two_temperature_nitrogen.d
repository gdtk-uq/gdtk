/**
 * Authors: Rowan G. and Peter J.
 * Date: 2017-07-13
 *
 * A model for single-species diatomic nitrogen with
 * two temperatures: a transrotational temperature and
 * a vibrational temperature.
 *
 */

module gas.two_temperature_nitrogen;

import std.math;
import std.conv;
import std.stdio;

import gas.gas_model;
import gas.physical_constants;

class TwoTemperatureNitrogen : GasModel {
public:
    this()
    {
	_n_species = 1;
	_n_modes = 1;
	_species_names.length = 1;
	_species_names[0] = "N2";
	create_species_reverse_lookup();
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "TwoTemperatureNitrogen";
	return to!string(repr);
    }

    override void update_thermo_from_pT(GasState Q) const
    {
	Q.rho = Q.p/(Q.Ttr*_R_N2);
	Q.u = (5./2.)*Q.Ttr*_R_N2;
	Q.e_modes[0] = (_R_N2*_theta_N2)/(exp(_theta_N2/Q.T_modes[0]) - 1.0);
    }

    override void update_thermo_from_rhou(GasState Q) const
    {
	Q.Ttr = Q.u/((5./2.)*_R_N2);
	Q.T_modes[0] = _theta_N2/log((_R_N2*_theta_N2/Q.e_modes[0]) + 1.0);
	Q.p = Q.rho*_R_N2*Q.Ttr;
    }

    override void update_thermo_from_rhoT(GasState Q) const
    {
	Q.p = Q.rho*_R_N2*Q.Ttr;
	Q.u = (5./2.)*Q.Ttr*_R_N2;
	Q.e_modes[0] = (_R_N2*_theta_N2)/(exp(_theta_N2/Q.T_modes[0]) - 1.0);
    }
    
    override void update_thermo_from_rhop(GasState Q) const
    {
	Q.Ttr = Q.p/(Q.rho*_R_N2);
	// Assume Q.T_modes[0] is set independently, and correct.
	Q.u = (5./2.)*Q.Ttr*_R_N2;
	Q.e_modes[0] = (_R_N2*_theta_N2)/(exp(_theta_N2/Q.T_modes[0]) - 1.0);
    }

    override void update_thermo_from_ps(GasState Q, double s) const
    {
	throw new GasModelException("update_thermo_from_ps not implemented in TwoTemperatureNitrogen.");
    }

    override void update_thermo_from_hs(GasState Q, double h, double s) const
    {
	throw new GasModelException("update_thermo_from_hs not implemented in TwoTemperatureNitrogen.");
    }

    override void update_sound_speed(GasState Q) const
    {
	Q.a = sqrt(_gamma*_R_N2*Q.Ttr);
    }

    override void update_trans_coeffs(GasState Q) const
    {
	// Not done correctly as we are only running inviscid flows presently.
	Q.mu = 0.0;
	Q.k = 0.0;
	Q.k_modes[0] = 0.0;
    }

    override double dudT_const_v(in GasState Q) const
    {
	return (5./2.)*_R_N2;
    }
    override double dhdT_const_p(in GasState Q) const
    {
	return (7./2.)*_R_N2;
    }
    override double dpdrho_const_T(in GasState Q) const
    {
	return _R_N2*Q.Ttr;
    }
    override double gas_constant(in GasState Q) const
    {
	return _R_N2;
    }
    override double internal_energy(in GasState Q) const
    {
	return Q.u;
    }
    override double enthalpy(in GasState Q) const
    {
	return Q.u + Q.p/Q.rho;
    }
    override double entropy(in GasState Q) const
    {
	throw new GasModelException("entropy not implemented in TwoTemperatureNitrogen.");
    }

private:
    double _R_N2 = 296.805; // gas constant for N2
    double _theta_N2 = 3354.0; // K, characteristic vib temp for N2
    double _gamma = 7./5.; // ratio of specific heats.
}

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
	double tau = (7.12e-9/pAtm)*exp(124.07/pow(Q.Ttr, 1./3.));
	// Find the total internal energy in the gas
	double uTotal = Q.u + Q.e_modes[0];
	// Find the vib energy at equilibrium with Ttr
	_Q_eq.Ttr = Q.Ttr;
	_Q_eq.T_modes[0] = Q.Ttr;
	_Q_eq.p = Q.p;
	_gmodel.update_thermo_from_pT(_Q_eq);
	double u_v_eq = _Q_eq.e_modes[0];
	double u_v = Q.e_modes[0];
	Q.e_modes[0] = u_v_eq + (u_v - u_v_eq)*exp(-tInterval/tau);
	Q.u = uTotal - Q.e_modes[0];
	_gmodel.update_thermo_from_rhou(Q);
	_gmodel.update_sound_speed(Q);
    }

private:
    GasState _Q_eq;
}
