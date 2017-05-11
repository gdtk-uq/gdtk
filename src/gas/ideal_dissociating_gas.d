/**
 * ideal_dissociating_gas.d
 *
 * This is the Lighthill-Freeman ideal dissociating gas model as described in
 * M. N. Macrossan (1990)
 * Hypervelocity flow of dissociating nitrogen downstream of a blunt nose.
 * Journal of Fluid Mechanics Vol. 217, pp 167-202.
 *
 * Authors: Peter J., Katrina Sklavos and Rowan G.
 * Version: 2017-04-22: initial shell copied from powers-aslam-gas module.
 *          2017-04-11: Filled in IDG thermo details as noted in 
 *                      PJ's workbooks 2017-04-22 through 2017-05-11
 */

module gas.ideal_dissociating_gas;

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

class IdealDissociatingGas: GasModel {
public:

    this(lua_State *L) {
	// Some parameters are fixed and some come from the gas model file.
	_n_species = 2;
	_n_modes = 0;
	_species_names.length = 2;
	// Bring table to TOS
	lua_getglobal(L, "IdealDissociatingGas");
	// [TODO] test that we actually have the table as item -1
	// Now, pull out the remaining value parameters.
	_species_names[0] = getString(L, -1, "molecule");
	_species_names[1] = getString(L, -1, "atom");
	_W = getDouble(L, -1, "W"); // molecular weight, g/mole
	_Rnn = R_universal / _W * 1000; // gas constant for molecule, J/kg/K
	_mol_masses.length = 2;
	_mol_masses[0] = _W/1000; // kg/mole
	_mol_masses[1] = _W/2000; // kg/mole
	_T_d = getDouble(L, -1, "T_d"); // characteristic dissociation temperature, K
	_rho_d = getDouble(L, -1, "rho_d"); // characteristic density, g/cm^3
	// Rate constants follow.
	_C1 = getDouble(L, -1, "C1"); 
	_n1 = getDouble(L, -1, "n1");
	_C2 = getDouble(L, -1, "C2");
	_n2 = getDouble(L, -1, "n2");
	lua_pop(L, 1); // dispose of the table
	create_species_reverse_lookup();
	// Reference conditions for entropy.
	_Tref = 298.15; _pref = 101.325e3;
    } // end constructor

    override string toString() const
    {
	char[] repr;
	repr ~= "IdealDissociatingGas =(";
	repr ~= "species=[\""~_species_names[0];
	repr ~= "\", \""~_species_names[1]~"\"]";
	repr ~= ", W=", to!string(_W);
	repr ~= ", Mmass=[" ~ to!string(_mol_masses[0]);
	repr ~= "," ~ to!string(_mol_masses[1]) ~ "]";
	repr ~= ", T_d=" ~ to!string(_T_d);
	repr ~= ", rho_d=" ~ to!string(_rho_d);
	repr ~= ", C1=" ~ to!string(_C1);
	repr ~= ", n1=" ~ to!string(_n1);
	repr ~= ", C2=" ~ to!string(_C2);
	repr ~= ", n2=" ~ to!string(_n2);
	repr ~= ")";
	return to!string(repr);
    }

    override void update_thermo_from_pT(GasState Q) const 
    {
	double alpha = Q.massf[1];
	Q.rho = Q.p/(Q.Ttr*_Rnn*(1+alpha));
	Q.u = _Rnn*alpha*_T_d + _Rnn*3*Q.Ttr;
    }
    override void update_thermo_from_rhoe(GasState Q) const
    {
	double alpha = Q.massf[1];
	Q.Ttr = (Q.u - _Rnn*alpha*_T_d)/(_Rnn*3);
	Q.p = Q.rho*(1+alpha)*_Rnn*Q.Ttr;
    }
    override void update_thermo_from_rhoT(GasState Q) const
    {
	double alpha = Q.massf[1];
	Q.p = Q.rho*(1+alpha)*_Rnn*Q.Ttr;
	Q.u = _Rnn*alpha*_T_d + _Rnn*3*Q.Ttr;
    }
    override void update_thermo_from_rhop(GasState Q) const
    {
	double alpha = Q.massf[1];
	Q.Ttr = Q.p/(Q.rho*(1+alpha)*_Rnn*Q.Ttr);
	Q.u = _Rnn*alpha*_T_d + _Rnn*3*Q.Ttr;
    }
    
    override void update_thermo_from_ps(GasState Q, double s) const
    {
	// For frozen composition.
	double alpha = Q.massf[1];
	Q.Ttr = _Tref * exp((1.0/((4+alpha)*_Rnn))*(s+(1+alpha)*_Rnn*log(Q.p/_pref)));
	update_thermo_from_pT(Q);
    }
    override void update_thermo_from_hs(GasState Q, double h, double s) const
    {
	// For frozen composition.
	double alpha = Q.massf[1];
	Q.Ttr = (h - _Rnn*alpha*_T_d) / ((4+alpha)*_Rnn);
	Q.p = _pref * exp(((4+alpha)*_Rnn*log(Q.Ttr/_Tref)-s)/((4+alpha)*_Rnn));
	update_thermo_from_pT(Q);
    }
    override void update_sound_speed(GasState Q) const
    {
	// For frozen composition.
	double alpha = Q.massf[1];
	double gamma = (4+alpha)/3.0;
	double Rgas = (1+alpha)*_Rnn;
	Q.a = sqrt(gamma*Rgas*Q.Ttr);
    }
    override void update_trans_coeffs(GasState Q)
    {
	// The gas is inviscid.
	Q.mu = 0.0;
	Q.k = 0.0;
    }
    override double dudT_const_v(in GasState Q) const
    {
	double alpha = Q.massf[1];
	return _Rnn*3; // frozen alpha
    }
    override double dhdT_const_p(in GasState Q) const
    {
	double alpha = Q.massf[1];
	return _Rnn*(4+alpha); // frozen alpha
    }
    override double dpdrho_const_T(in GasState Q) const
    {
	double alpha = Q.massf[1];
	return _Rnn*(1+alpha)*Q.Ttr; // frozen alpha
    }
    override double gas_constant(in GasState Q) const
    {
	double alpha = Q.massf[1];
	return _Rnn*(1+alpha);
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
	// Presume that we have a fixed composition mixture and that
	// Entropy for each species is zero at reference condition.
	double alpha = Q.massf[1];
	return (4+alpha)*log(Q.Ttr/_Tref) - (1+alpha)*_Rnn*log(Q.p/_pref);
    }

private:
    // Thermodynamic constants
    double _W; // Molecular mass in g/mole
    double _Rnn; // gas constant for molecule in J/kg/K
    double _T_d; // characteristic dissociation temperature, K
    double _rho_d; // characteristic density, g/cm^3
    double _Tref, _pref; // Reference conditions for entropy
    // Rate constants
    double _C1, _n1, _C2, _n2;
   // Molecular transport coefficents are zero.
} // end class IdealDissociatingGas


// Now, for the reaction update...
//
// It is included here because it is a small amount of code and
// it is specific to this particular gas model.

final class UpdateIDG : ThermochemicalReactor {
    
    this(string fname, GasModel gmodel)
    {
	super(gmodel); // hang on to a reference to the gas model
	// We need to pick a number of pieces out of the gas-model file, again.
	// Although they exist in the GasModel object, they are private.
	auto L = init_lua_State();
	doLuaFile(L, fname);
	lua_getglobal(L, "IdealDissociatingGas");
	_W = getDouble(L, -1, "W"); // molecular weight, g/mole
	_T_d = getDouble(L, -1, "T_d"); // characteristic dissociation temperature, K
	_rho_d = getDouble(L, -1, "rho_d"); // characteristic density, g/cm^3
	// Rate constants follow.
	_C1 = getDouble(L, -1, "C1"); 
	_n1 = getDouble(L, -1, "n1");
	_C2 = getDouble(L, -1, "C2");
	_n2 = getDouble(L, -1, "n2");
	lua_pop(L, 1); // dispose of the table
	lua_close(L);
    }
    
    override void opCall(GasState Q, double tInterval, ref double dtSuggest,
			 ref double[] params)
    {
	double alpha = Q.massf[1];
	//
	// Since the internal energy and density in the (isolated) reactor is fixed,
	// we need to evaluate the new temperature, pressure, etc.
	_gmodel.update_thermo_from_rhoe(Q);
	_gmodel.update_sound_speed(Q);
    }

private:
    // Thermodynamic constants
    double _W; // Molecular mass in g/mole
    double _T_d; // characteristic dissociation temperature, K
    double _rho_d; // characteristic density, g/cm^3
    // Rate constants
    double _C1, _n1, _C2, _n2;
} // end class UpdateIDG


// Unit test of the basic gas model...

version(ideal_dissociating_gas_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
	lua_State* L = init_lua_State();
	doLuaFile(L, "sample-data/idg-nitrogen.lua");
	auto gm = new IdealDissociatingGas(L);
	lua_close(L);
	auto gd = new GasState(2, 0);
	gd.p = 1.0e5;
	gd.Ttr = 300.0;
	gd.massf[0] = 1.0; gd.massf[1] = 0.0;
	double Rgas = R_universal / 0.028;
	assert(approxEqual(gm.R(gd), Rgas, 1.0e-4), failedUnitTest());
	assert(gm.n_modes == 0, failedUnitTest());
	assert(gm.n_species == 2, failedUnitTest());
	assert(approxEqual(gd.p, 1.0e5, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.Ttr, 300.0, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.massf[0], 1.0, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.massf[1], 0.0, 1.0e-6), failedUnitTest());

	gm.update_thermo_from_pT(gd);
	gm.update_sound_speed(gd);
	double my_rho = 1.0e5 / (Rgas * 300.0);
	assert(approxEqual(gd.rho, my_rho, 1.0e-4), failedUnitTest());
	double my_Cv = gm.dudT_const_v(gd);
	double my_u = my_Cv*300.0; 
	assert(approxEqual(gd.u, my_u, 1.0e-3), failedUnitTest());
	double my_Cp = gm.dhdT_const_p(gd);
	double my_a = sqrt(my_Cp/my_Cv*Rgas*300.0);
	assert(approxEqual(gd.a, my_a, 1.0e-3), failedUnitTest());
	gm.update_trans_coeffs(gd);
	assert(approxEqual(gd.mu, 0.0, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.k, 0.0, 1.0e-6), failedUnitTest());

	return 0;
    }
}
