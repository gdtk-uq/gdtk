/**
 * vib_specific_nitrogen.d
 * Authors: Rowan G., Katrina Sklavos and Peter J.
 *
 * This is a 10-level vibrationally-specific model for nitrogen
 * as descrbied in:
 *
 * Giordano, et al. (1997)
 * Vibrationally Relaxing Flow of N2 past an Infinite Cylinder
 * Journal of Thermophysics and Heat Transfer, vol 11, no 1, pp 27 - 35
 *
 */

module gas.vib_specific_nitrogen;

import std.algorithm.iteration;

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

// [TODO:KS]
// Presently the number of vibrational levels
// is hard-coded. We'd like to test at small number
// then increase this to 10.
// Eventually, we might remove the hard-coded number.
immutable int N_VIB_LEVELS = 3;

class VibSpecificNitrogen: GasModel {
public:
    this()
    {
	_n_species = 1;
	_n_modes = N_VIB_LEVELS;
	_species_names.length = 1;
	_species_names[0] = "N2";
	create_species_reverse_lookup();
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "VibSpecificNitrogen=()";
	return to!string(repr);
    }

    override void update_thermo_from_pT(GasState Q) const 
    {
	// [TODO:KS]
	// Fill in Q.rho and Q.u and Q.u_modes
	// Q.rho = ....;
	// Q.u = ....;
	// foreach (imode; 0 .. _n_modes) {
	//     Q.u_modes = ....;
	// }
    }
    override void update_thermo_from_rhou(GasState Q) const
    {
	// [TODO:KS]
	// Fill in Q.p, Q.Ttr and Q.T_modes
    }
    override void update_thermo_from_rhoT(GasState Q) const
    {
	// [TODO:KS]
	// Fill in Q.p, Q.u and Q.u_modes 

    }
    override void update_thermo_from_rhop(GasState Q) const
    {
	// [TODO:KS]
	// Fill in Q.Ttr, Q.u
    }
    
    override void update_thermo_from_ps(GasState Q, double s) const
    {
	throw new Error("ERROR: VibSpecificNitrogen:update_thermo_from_ps NOT IMPLEMENTED.");

    }

    override void update_thermo_from_hs(GasState Q, double h, double s) const
    {
	throw new Error("ERROR: VibSpecificNitrogen:update_thermo_from_hs NOT IMPLEMENTED.");

    }
    override void update_sound_speed(GasState Q) const
    {
	// [TODO:KS]
	// Fill in Q.a
    }
    override void update_trans_coeffs(GasState Q)
    {
	// The gas is inviscid.
	Q.mu = 0.0;
	Q.k = 0.0;
    }
    override double dudT_const_v(in GasState Q) const
    {
	throw new Error("ERROR: VibSpecificNitrogen:dudT_const_v NOT IMPLEMENTED.");
    }
    override double dhdT_const_p(in GasState Q) const
    {
	throw new Error("ERROR: VibSpecificNitrogen:dhdT_const_p NOT IMPLEMENTED.");
    }
    override double dpdrho_const_T(in GasState Q) const
    {
	throw new Error("ERROR: VibSpecificNitrogen:dpdrho_const_T NOT IMPLEMENTED.");
    }
    override double gas_constant(in GasState Q) const
    {
	return _R_N2;
    }
    override double internal_energy(in GasState Q) const
    {
	return Q.u + sum(Q.u_modes);
    }
    override double enthalpy(in GasState Q) const
    {
	return Q.u + sum(Q.u_modes) + Q.p/Q.rho;
    }
    override double entropy(in GasState Q) const
    {
	throw new Error("ERROR: VibSpecificNitrogen:entropy NOT IMPLEMENTED.");
    }

private:
    double _R_N2 = 296.805; // gas constant for N2
    // [TODO:KS]
    // You may declare any extra data you need in this area.

} // end class VibSpecificNitrogen

version(vib_specific_nitrogen_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
	auto gm = new VibSpecificNitrogen();
	auto Q = new GasState(1, N_VIB_LEVELS);
	Q.p = 1.0e5;
	Q.Ttr = 300.0;
	foreach (imode; 0 .. gm.n_modes()) {
	    Q.T_modes[imode] = 300.0;
	}
	Q.massf[0] = 1.0;

	// [TODO:KS]
	// Test update_thermo_from_pT
	gm.update_thermo_from_pT(Q);
	// assert(approxEqual(Q.rho, ...., 1.0e-6), failedUnitTest());
	
	// and other tests.
	return 0;
    }
}
