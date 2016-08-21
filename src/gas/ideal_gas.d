/**
 * idealgas.d
 * Ideal gas model for use in the CFD codes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-06-22: initial cut, to explore options.
 */

module gas.ideal_gas;

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

class IdealGas: GasModel {
public:

    this(lua_State *L) {
	_n_species = 1;
	_n_modes = 1;
	// Bring table to TOS
	lua_getglobal(L, "IdealGas");
	// Let's overwrite that here.
	_species_names.length = 1;
	_species_names[0] = getString(L, -1, "speciesName");
	// Now, pull out the remaining numeric value parameters.
	_mol_masses.length = 1;
	_mol_masses[0] = getDouble(L, -1, "mMass"); // NOTE: changed from append to assignment, anything broken by this??
	_gamma = getDouble(L, -1, "gamma");
	// Reference values for entropy
	lua_getfield(L, -1, "entropyRefValues");
	_s1 = getDouble(L, -1, "s1");
	_T1 = getDouble(L, -1, "T1");
	_p1 = getDouble(L, -1, "p1");
	lua_pop(L, 1);
	// Molecular transport coefficent models.
	lua_getfield(L, -1, "viscosity");
	_viscModel = createViscosityModel(L);
	lua_pop(L, 1);
	
	lua_getfield(L, -1, "thermCondModel");
	auto model = getString(L, -1, "model");
	if ( model == "constPrandtl" ) {
	    _Prandtl = getDouble(L, -1, "Prandtl");
	    _constPrandtl = true;
	}
	else {
	    _thermCondModel = createThermalConductivityModel(L);
	    _constPrandtl = false;
	}
	lua_pop(L, 1);
	// Compute derived parameters
	_Rgas = R_universal/_mol_masses[0];
	_Cv = _Rgas / (_gamma - 1.0);
	_Cp = _Rgas*_gamma/(_gamma - 1.0);
	create_species_reverse_lookup();
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "IdealGas =(";
	repr ~= "name=\"" ~ _species_names[0] ~"\"";
	repr ~= ", Mmass=" ~ to!string(_mol_masses[0]);
	repr ~= ", gamma=" ~ to!string(_gamma);
	repr ~= ", s1=" ~ to!string(_s1);
	repr ~= ", T1=" ~ to!string(_T1);
	repr ~= ", p1=" ~ to!string(_p1);
	repr ~= ", constPrandtl=" ~ to!string(_constPrandtl);
	repr ~= ", Prandtl=" ~ to!string(_Prandtl);
	repr ~= ")";
	return to!string(repr);
    }

    override void update_thermo_from_pT(GasState Q) const 
    {
	assert(Q.T.length == 1, "incorrect length of temperature array");
	Q.rho = Q.p/(Q.T[0]*_Rgas);
	Q.e[0] = _Cv*Q.T[0];
    }
    override void update_thermo_from_rhoe(GasState Q) const
    {
	assert(Q.e.length == 1, "incorrect length of energy array");
	Q.T[0] = Q.e[0]/_Cv;
	Q.p = Q.rho*_Rgas*Q.T[0];
    }
    override void update_thermo_from_rhoT(GasState Q) const
    {
	assert(Q.T.length == 1, "incorrect length of temperature array");
	Q.p = Q.rho*_Rgas*Q.T[0];
	Q.e[0] = _Cv*Q.T[0];
    }
    override void update_thermo_from_rhop(GasState Q) const
    {
	assert(Q.T.length == 1, "incorrect length of temperature array");
	Q.T[0] = Q.p/(Q.rho*_Rgas);
	Q.e[0] = _Cv*Q.T[0];
	
    }
    
    override void update_thermo_from_ps(GasState Q, double s) const
    {
	Q.T[0] = _T1 * exp((1.0/_Cp)*((s - _s1) + _Rgas * log(Q.p/_p1)));
	update_thermo_from_pT(Q);
    }
    override void update_thermo_from_hs(GasState Q, double h, double s) const
    {
	Q.T[0] = h / _Cp;
	Q.p = _p1 * exp((1.0/_Rgas)*(_s1 - s + _Cp*log(Q.T[0]/_T1)));
	update_thermo_from_pT(Q);
    }
    override void update_sound_speed(GasState Q) const
    {
	Q.a = sqrt(_gamma*_Rgas*Q.T[0]);
    }
    override void update_trans_coeffs(GasState Q)
    {
	_viscModel.update_viscosity(Q);
	if ( _constPrandtl ) {
	    Q.k[0] = _Cp*Q.mu/_Prandtl;
	}
	else {
	    _thermCondModel.update_thermal_conductivity(Q);
	}
    }
    /*
    override void eval_diffusion_coefficients(ref GasState Q) {
	throw new Exception("not implemented");
    }
    */
    override double dedT_const_v(in GasState Q) const
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
	return R*Q.T[0];
    }
    override double gas_constant(in GasState Q) const
    {
	return R_universal/_mol_masses[0];
    }
    override double internal_energy(in GasState Q) const
    {
	return Q.e[0];
    }
    override double enthalpy(in GasState Q) const
    {
	return Q.e[0] + Q.p/Q.rho;
    }
    override double entropy(in GasState Q) const
    {
	return _s1 + _Cp * log(Q.T[0]/_T1) - _Rgas * log(Q.p/_p1);
    }

private:
    // Thermodynamic constants
    double _Rgas = R_universal/0.02896; // J/kg/K
    double _gamma = 1.4;   // ratio of specific heats
    double _Cv = R_universal/0.02896 / 0.4; // J/kg/K
    double _Cp = R_universal/0.02896 * 1.4/0.4; // J/kg/K
    // Reference values for entropy
    double _s1 = 0.0; // J/kg/K
    double _T1 = 298.15; // K
    double _p1 = 101.325e3; // Pa
    // Molecular transport coefficent constants.
    Viscosity _viscModel;
    // We compute thermal conductivity in one of two ways:
    // 1. based on constant Prandtl number; OR
    // 2. ThermalConductivity model
    // therefore we have places for both data
    bool _constPrandtl = false;
    double _Prandtl = 1.0;
    ThermalConductivity _thermCondModel;

} // end class Ideal_gas

version(ideal_gas_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
	lua_State* L = init_lua_State("sample-data/ideal-air-gas-model.lua");
	auto gm = new IdealGas(L);
	lua_close(L);
	auto gd = new GasState(1, 1);
	gd.p = 1.0e5;
	gd.T[0] = 300.0;
	gd.massf[0] = 1.0;
	assert(approxEqual(gm.R(gd), 287.086, 1.0e-4), failedUnitTest());
	assert(gm.n_modes == 1, failedUnitTest());
	assert(gm.n_species == 1, failedUnitTest());
	assert(approxEqual(gd.p, 1.0e5, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.T[0], 300.0, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.massf[0], 1.0, 1.0e-6), failedUnitTest());

	gm.update_thermo_from_pT(gd);
	gm.update_sound_speed(gd);
	assert(approxEqual(gd.rho, 1.16109, 1.0e-4), failedUnitTest());
	assert(approxEqual(gd.e[0], 215314.0, 1.0e-4), failedUnitTest());
	assert(approxEqual(gd.a, 347.241, 1.0e-4), failedUnitTest());
	gm.update_trans_coeffs(gd);
	assert(approxEqual(gd.mu, 1.84691e-05, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.k[0], 0.0262449, 1.0e-6), failedUnitTest());

	return 0;
    }
}
