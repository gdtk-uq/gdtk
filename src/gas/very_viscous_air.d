/**
 * very_viscous_air.d
 * A special ideal air model used for the Method of Manufactured Solutions
 * viscous case.
 *
 * Author: Rowan G.
 * Version: 2015-05-05
 */

module gas.very_viscous_air;

import gas.gas_model;
import gas.physical_constants;
import std.math;
import std.stdio;
import std.string;
import std.file;
import std.json;
import std.conv;
import std.c.stdlib : exit;
import util.msg_service;
import util.lua;

class VeryViscousAir: GasModel {
public:
    this() {
	// Default model is mostly initialized in the private data below.
	_n_species = 1;
	_n_modes = 1;
	_species_names ~= "very viscous";
	_Rgas = 287.0;
	_mol_masses ~= R_universal/_Rgas;
	_gamma = 1.4;
	_Cv = _Rgas / (_gamma - 1.0);
	_Cp = _Rgas*_gamma/(_gamma - 1.0);
	_mu = 10.0;
	double Pr = 1.0;
	_k = _mu * _Cp / Pr;
    }

    this(lua_State* L)
    {
	this();
	lua_getglobal(L, "VeryViscousAir");
	// Possibly override k and mu
	lua_getfield(L, -1, "mu");
	if ( !lua_isnil(L, -1) ) {
	    _mu = to!double(lua_tonumber(L, -1));
	}
	lua_pop(L, 1);
	lua_getfield(L, -1, "k");
	if ( !lua_isnil(L, -1) ) {
	    _k = to!double(lua_tonumber(L, -1));
	}
	lua_pop(L, 1);
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "VeryViscousAir =()";
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
	throw new Exception(format("Not implemented: line=%d, file=%s\n", __LINE__, __FILE__));
    }
    override void update_thermo_from_hs(GasState Q, double h, double s) const
    {
	throw new Exception(format("Not implemented: line=%d, file=%s\n", __LINE__, __FILE__));
    }
    override void update_sound_speed(GasState Q) const
    {
	Q.a = sqrt(_gamma*_Rgas*Q.T[0]);
    }
    override void update_trans_coeffs(GasState Q) const
    {
	Q.mu = _mu;
	Q.k[0] = _k;
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
	return _Rgas;
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
    double _Rgas = 287.0; // J/kg/K
    double _gamma = 1.4;   // ratio of specific heats
    double _Cv = R_universal/0.02896 / 0.4; // J/kg/K
    double _Cp = R_universal/0.02896 * 1.4/0.4; // J/kg/K
    // Reference values for entropy
    double _s1 = 0.0; // J/kg/K
    double _T1 = 298.15; // K
    double _p1 = 101.325e3; // Pa
    // Molecular transport coefficent constants.
    double _mu;
    double _k;

} // end class Very_viscous_air

unittest {
    auto gm = new VeryViscousAir();
    auto gs = new GasState(gm, 100.0e3, 300.0);
    assert(approxEqual(gm.R(gs), 287.0), failedUnitTest());

    gm.update_thermo_from_pT(gs);
    gm.update_sound_speed(gs);
    assert(approxEqual(gs.rho, 1.16144), failedUnitTest());
    assert(approxEqual(gs.e[0], 215250), failedUnitTest());
    assert(approxEqual(gs.a, 347.189), failedUnitTest());
    gm.update_trans_coeffs(gs);
    assert(approxEqual(gs.mu, 10.0), failedUnitTest());
    assert(approxEqual(gs.k[0], 10045), failedUnitTest());
}
