/**
 * powers_aslam_gas.d
 *
 * Two-component reacting gas as described in.
 * JM Powers and TD Aslam (2006)
 * Exact solution for multidimensional compressible reactive flow
 * for verifying numerical algorithms.
 * AIAA Journal Vol. 44 No. 2 pages 337-344
 *
 * This gas model is useful as a demonstration of building a custom
 * reacting gas model and as a flow-solver verification tool.
 *
 * Authors: Peter J. and Rowan G.
 * Version: 2016-01-07: initial cut.
 */

module gas.powers_aslam_gas;

import gas.gas_model;
import gas.physical_constants;
import gas.diffusion.viscosity;
import gas.diffusion.therm_cond;
import std.math;
import std.stdio;
import std.string;
import std.file;
import std.json;
import std.conv;
import util.lua;
import util.lua_service;
import core.stdc.stdlib : exit;

// First, the basic gas model.

class PowersAslamGas: GasModel {
public:

    this(lua_State *L) {
	// Some parameters are fixed and some come from the gas model file.
	_n_species = 2;
	_n_modes = 0;
	_species_names.length = 2;
	_species_names[0] = "A";
	_species_names[1] = "B";
	// Bring table to TOS
	lua_getglobal(L, "PowersAslamGas");
	// [TODO] test that we actually have the table as item -1
	// Now, pull out the remaining numeric value parameters.
	_Rgas = getDouble(L, -1, "R");
	_mol_masses.length = 2;
	_mol_masses[0] = R_universal / _Rgas;
	_mol_masses[1] = _mol_masses[0];
	_gamma = getDouble(L, -1, "gamma");
	// Heat of reaction
	_q = getDouble(L, -1, "q");
	_alpha = getDouble(L, -1, "alpha");
	_Ti = getDouble(L, -1, "Ti");
	lua_pop(L, 1); // dispose of the table
	// Entropy reference, same as for IdealAir
	_s1 = 0.0;
	_T1 = 298.15;
	_p1 = 101.325e3;
	// Compute derived parameters
	_Cv = _Rgas / (_gamma - 1.0);
	_Cp = _Rgas*_gamma/(_gamma - 1.0);
	create_species_reverse_lookup();
    } // end constructor

    override string toString() const
    {
	char[] repr;
	repr ~= "PowersAslamGas =(";
	repr ~= "species=[\"A\", \"B\"]";
	repr ~= ", Mmass=[" ~ to!string(_mol_masses[0]);
	repr ~= "," ~ to!string(_mol_masses[1]) ~ "]";
	repr ~= ", gamma=" ~ to!string(_gamma);
	repr ~= ", q=" ~ to!string(_q);
	repr ~= ", alpha=" ~ to!string(_alpha);
	repr ~= ", Ti=" ~ to!string(_Ti);
	repr ~= ")";
	return to!string(repr);
    }

    override void update_thermo_from_pT(GasState Q) const 
    {
	Q.rho = Q.p/(Q.Ttr*_Rgas);
	Q.u = _Cv*Q.Ttr - Q.massf[1]*_q;
    }
    override void update_thermo_from_rhoe(GasState Q) const
    {
	Q.Ttr = (Q.u + Q.massf[1]*_q)/_Cv;
	Q.p = Q.rho*_Rgas*Q.Ttr;
    }
    override void update_thermo_from_rhoT(GasState Q) const
    {
	Q.p = Q.rho*_Rgas*Q.Ttr;
	Q.u = _Cv*Q.Ttr - Q.massf[1]*_q;
    }
    override void update_thermo_from_rhop(GasState Q) const
    {
	Q.Ttr = Q.p/(Q.rho*_Rgas);
	Q.u = _Cv*Q.Ttr - Q.massf[1]*_q;
    }
    
    override void update_thermo_from_ps(GasState Q, double s) const
    {
	// [FIX-ME] to account for the reaction
	Q.Ttr = _T1 * exp((1.0/_Cp)*((s - _s1) + _Rgas * log(Q.p/_p1)));
	update_thermo_from_pT(Q);
    }
    override void update_thermo_from_hs(GasState Q, double h, double s) const
    {
	// [FIX-ME] to account for the reaction
	Q.Ttr = h / _Cp;
	Q.p = _p1 * exp((1.0/_Rgas)*(_s1 - s + _Cp*log(Q.Ttr/_T1)));
	update_thermo_from_pT(Q);
    }
    override void update_sound_speed(GasState Q) const
    {
	Q.a = sqrt(_gamma*_Rgas*Q.Ttr);
    }
    override void update_trans_coeffs(GasState Q)
    {
	// The gas is inviscid for the test cases described in the AIAA paper.
	Q.mu = 0.0;
	Q.k = 0.0;
    }
    override double dudT_const_v(in GasState Q) const
    {
	return _Cv;
    }
    override double dhdT_const_p(in GasState Q) const
    {
	return _Cp;
    }
    override double dpdrho_const_T(in GasState Q) const
    {
	double R = gas_constant(Q);
	return R*Q.Ttr;
    }
    override double gas_constant(in GasState Q) const
    {
	return _Rgas;
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
	// [FIX-ME] to account for the reaction
	return _s1 + _Cp * log(Q.Ttr/_T1) - _Rgas * log(Q.p/_p1);
    }

private:
    // Thermodynamic constants
    double _Rgas; // J/kg/K
    double _gamma;   // ratio of specific heats
    double _Cv; // J/kg/K
    double _Cp; // J/kg/K
    // Reference values for entropy
    double _s1;  // J/kg/K
    double _T1;  // K
    double _p1;  // Pa
    // Molecular transport coefficents are zero.
    // Heat of reaction.
    double _q; // J/kg
    // Reaction rate constant
    double _alpha; // 1/s
    // Ignition temperature
    double _Ti; // degrees K
} // end class PowersAslamGas


// Now, for the reaction update...
//
// It is included here because it is a small amount of code and
// is specific to the gas model.

final class UpdateAB : ThermochemicalReactor {
    
    this(string fname, GasModel gmodel)
    {
	super(gmodel); // hang on to a reference to the gas model
	// We need to pick a number of pieces out of the gas-model file, again.
	// Although they exist in the GasModel object, they are private.
	auto L = init_lua_State(fname);
	lua_getglobal(L, "PowersAslamGas");
	// Now, pull out the numeric value parameters.
	_alpha = getDouble(L, -1, "alpha");
	_Ti = getDouble(L, -1, "Ti");
	lua_pop(L, 1); // dispose of the table
	lua_close(L);
    }
    
    override void opCall(GasState Q, double tInterval, ref double dtSuggest)
    {
	if (Q.Ttr > _Ti) {
	    // We are above the ignition point, proceed with reaction.
	    double massfA = Q.massf[0];
	    double massfB = Q.massf[1];
	    // This gas has a very simple reaction scheme that can be integrated explicitly.
	    massfA = massfA*exp(-_alpha*tInterval);
	    massfB = 1.0 - massfA;
	    Q.massf[0] = massfA; Q.massf[1] = massfB;
	} else {
	    // do nothing, since we are below the ignition temperature
	}
	// Since the internal energy and density in the (isolated) reactor is fixed,
	// we need to evaluate the new temperature, pressure, etc.
	_gmodel.update_thermo_from_rhoe(Q);
	_gmodel.update_sound_speed(Q);
    }

private:
    // Reaction rate constant
    double _alpha; // 1/s
    // Ignition temperature
    double _Ti; // degrees K
} // end class UpdateAB


// Unit test of the basic gas model...

version(powers_aslam_gas_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
	lua_State* L = init_lua_State("sample-data/powers-aslam-gas-model.lua");
	auto gm = new PowersAslamGas(L);
	lua_close(L);
	auto gd = new GasState(2, 0);
	gd.p = 1.0e5;
	gd.Ttr = 300.0;
	gd.massf[0] = 0.75; gd.massf[1] = 0.25;
	assert(approxEqual(gm.R(gd), 287.0, 1.0e-4), failedUnitTest());
	assert(gm.n_modes == 0, failedUnitTest());
	assert(gm.n_species == 2, failedUnitTest());
	assert(approxEqual(gd.p, 1.0e5, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.Ttr, 300.0, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.massf[0], 0.75, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.massf[1], 0.25, 1.0e-6), failedUnitTest());

	gm.update_thermo_from_pT(gd);
	gm.update_sound_speed(gd);
	double my_rho = 1.0e5 / (287.0 * 300.0);
	assert(approxEqual(gd.rho, my_rho, 1.0e-4), failedUnitTest());
	double my_Cv = gm.dudT_const_v(gd);
	double my_u = my_Cv*300.0 - 0.25*300000.0; 
	assert(approxEqual(gd.u, my_u, 1.0e-3), failedUnitTest());
	double my_Cp = gm.dhdT_const_p(gd);
	double my_a = sqrt(my_Cp/my_Cv*287.0*300.0);
	assert(approxEqual(gd.a, my_a, 1.0e-3), failedUnitTest());
	gm.update_trans_coeffs(gd);
	assert(approxEqual(gd.mu, 0.0, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.k, 0.0, 1.0e-6), failedUnitTest());

	return 0;
    }
}
