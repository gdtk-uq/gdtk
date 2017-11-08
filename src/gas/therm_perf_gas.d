/**
 * therm_perf_gas_mix.d
 * Thermally perfect gas mix model for use in the CFD codes.
 *
 * Author: Rowan G. and Peter J.
 * First code: 27-Jan-2015
 */

module gas.therm_perf_gas;

import std.math;
import std.stdio;
import std.string;
import std.conv : to;
import util.lua;
import util.lua_service;
import nm.brent; 
import nm.bracketing;
import core.stdc.stdlib : exit;

import gas.gas_model;
import gas.physical_constants;
import gas.thermo.cea_thermo_curves;
import gas.thermo.perf_gas_mix_eos;
import gas.thermo.therm_perf_gas_mix_eos;
import gas.diffusion.viscosity;
import gas.diffusion.therm_cond;
import gas.diffusion.cea_viscosity;
import gas.diffusion.cea_therm_cond;
import gas.diffusion.sutherland_viscosity;
import gas.diffusion.sutherland_therm_cond;
import gas.diffusion.wilke_mixing_viscosity;
import gas.diffusion.wilke_mixing_therm_cond;

class ThermallyPerfectGas: GasModel {
public:
    this(lua_State* L)
    {
	getArrayOfStrings(L, LUA_GLOBALSINDEX, "species", _species_names);
	_n_species = cast(uint) _species_names.length;
	_n_modes = 0;
	
	// 0. Initialise private work arrays
	_Cp.length = _n_species;
	_Cv.length = _n_species;
	_h.length = _n_species;
	_s.length = _n_species;
	// 1. Initialise gas constants from molecular mass
	_R.length = _n_species;
	_mol_masses.length = _n_species;
	foreach (isp; 0 .. _n_species) {
	    lua_getglobal(L, _species_names[isp].toStringz);
	    _mol_masses[isp] = getDouble(L, -1, "M");
	    _R[isp] = R_universal/_mol_masses[isp];
	    lua_pop(L, 1);
	}
	// 1b. Gather L-J parameters
	_LJ_sigmas.length = _n_species;
	_LJ_epsilons.length = _n_species;
	foreach (isp; 0 .. _n_species) {
	    lua_getglobal(L, _species_names[isp].toStringz);
	    _LJ_sigmas[isp] = getDouble(L, -1, "sigma");
	    _LJ_epsilons[isp] = getDouble(L, -1, "epsilon");
	    lua_pop(L, 1);
	}
	// 2. Set the p-v-T EOS
	_pgMixEOS = new PerfectGasMixEOS(_R);
	// 3. Set the e-v-T EOS (which in this case is an e-T relationship only)
	foreach ( isp; 0.._n_species ) {
	    lua_getglobal(L, _species_names[isp].toStringz);
	    lua_getfield(L, -1, "thermoCoeffs");
	    _curves ~= createCEAThermo(L, _R[isp]);
	    lua_pop(L, 1);
	    lua_pop(L, 1);
	}
	_tpgMixEOS = new ThermallyPerfectGasMixEOS(_R, _curves);
	// 4. Set the viscosity model
	Viscosity[] vms;
	foreach ( isp; 0.._n_species ) {
	    lua_getglobal(L, _species_names[isp].toStringz);
	    lua_getfield(L, -1, "viscosity");
	    lua_getfield(L, -1, "model");
	    string model = to!string(luaL_checkstring(L, -1));
	    lua_pop(L, 1);
	    switch (model) {
	    case "CEA":
		vms ~= createCEAViscosity(L);
		break;
	    case "Sutherland":
		vms ~= createSutherlandViscosity(L);
		break;
	    default:
		string errMsg = format("The viscosity model '%s' is not available.\n", model);
		errMsg ~= format("This error occurred for species %d when constructing the thermally perfect gas mix.\n");
		throw new Error(errMsg);
	    }
	    lua_pop(L, 1);
	    lua_pop(L, 1);
	}
	_viscModel = new WilkeMixingViscosity(vms, _mol_masses);
	// 5. Set the thermal conductivity model
	ThermalConductivity[] tcms;
	foreach (isp; 0.._n_species ) {
	    lua_getglobal(L, _species_names[isp].toStringz);
	    lua_getfield(L, -1, "thermal_conductivity");
	    lua_getfield(L, -1, "model");
	    string model = to!string(luaL_checkstring(L, -1));
	    lua_pop(L, 1);
	    switch (model) {
	    case "CEA":
		tcms ~= createCEAThermalConductivity(L);
		break;
	    case "Sutherland":
		tcms ~=  createSutherlandThermalConductivity(L);
		break;
	    default:
		string errMsg = format("The thermal conductivity model '%s' is not available.\n", model);
		errMsg ~= format("This error occurred for species %d when constructing the thermally perfect gas mix.\n");
		throw new Error(errMsg);
	    }
	    lua_pop(L, 1);
	    lua_pop(L, 1);
	}
	_thermCondModel = new WilkeMixingThermCond(tcms, _mol_masses);
	create_species_reverse_lookup();
    }

    this(in string fname)
    {
	auto L = init_lua_State();
	doLuaFile(L, fname);
	this(L);
	lua_close(L);
	create_species_reverse_lookup();
    }

    override string toString() const
    {
	return "";
    }

    override void update_thermo_from_pT(GasState Q) 
    {
	_pgMixEOS.update_density(Q);
	_tpgMixEOS.update_energy(Q);
    }
    override void update_thermo_from_rhou(GasState Q)
    {
	_tpgMixEOS.update_temperature(Q);
	_pgMixEOS.update_pressure(Q);
    }
    override void update_thermo_from_rhoT(GasState Q)
    {
	_tpgMixEOS.update_energy(Q);
	_pgMixEOS.update_pressure(Q);
    }
    override void update_thermo_from_rhop(GasState Q)
    {
	_pgMixEOS.update_temperature(Q);
	_tpgMixEOS.update_energy(Q);
    }
    override void update_thermo_from_ps(GasState Q, double s)
    {
	double TOL = 1.0e-6;
	double delT = 100.0;
	double T1 = Q.Ttr;
	// It's possible that Ttr is 'nan' if it's never been set, so:
	if ( isNaN(T1) )
	    T1 = 300.0; // Just set at some value to get started.
	double Tsave = T1;
	double T2 = T1 + delT;

	auto zeroFun = delegate (double T) {
	    Q.Ttr = T;
	    double s_guess = entropy(Q);
	    return s - s_guess;
	};

	if ( bracket!zeroFun(T1, T2) == -1 ) {
	    string msg = "The 'bracket' function failed to find temperature values\n";
	    msg ~= "that bracketed the zero function in ThermallyPerfectGas.update_thermo_from_ps().\n";
	    msg ~= format("The final values are: T1 = %12.6f and T2 = %12.6f\n", T1, T2);
	    throw new Exception(msg);
	}

	if ( T1 < T_MIN )
	    T1 = T_MIN;
	
	try {
	    Q.Ttr = solve!zeroFun(T1, T2, TOL);
	}
	catch ( Exception e ) {
	    string msg = "There was a problem iterating to find temperature\n";
	    msg ~= "in function ThermallyPerfectGas.update_thermo_from_ps().\n";
	    msg ~= format("The initial temperature guess was: %12.6f\n", Tsave);
	    msg ~= format("The target entropy value was: %12.6f\n", s);
	    msg ~= format("The GasState is currently:\n");
	    msg ~= Q.toString();
	    msg ~= "The message from the ridder.solve function is:\n";
	    msg ~= e.msg;
	    throw new Exception(msg);
	}
	_tpgMixEOS.update_energy(Q);
	_pgMixEOS.update_density(Q);
    }
    override void update_thermo_from_hs(GasState Q, double h, double s)
    {
	// We do this in two stages.
	// First, from enthalpy we compute temperature.
	double TOL = 1.0e-6;
	double delT = 100.0;
	double T1 = Q.Ttr;
	// It's possible that Ttr is 'nan' if it's never been set, so:
	if ( isNaN(T1) )
	    T1 = 300.0; // Just set at some value to get started.
	double Tsave = T1;
	double T2 = T1 + delT;

	auto zeroFun = delegate (double T) {
	    Q.Ttr = T;
	    double h_guess = enthalpy(Q);
	    return h - h_guess;
	};

	if ( bracket!zeroFun(T1, T2) == -1 ) {
	    string msg = "The 'bracket' function failed to find temperature values\n";
	    msg ~= "that bracketed the zero function in ThermallyPerfectGas.update_thermo_from_hs().\n";
	    msg ~= format("The final values are: T1 = %12.6f and T2 = %12.6f\n", T1, T2);
	    throw new Exception(msg);
	}

	if ( T1 < T_MIN )
	    T1 = T_MIN;
	
	try {
	    Q.Ttr = solve!zeroFun(T1, T2, TOL);
	}
	catch ( Exception e ) {
	    string msg = "There was a problem iterating to find temperature\n";
	    msg ~= "in function ThermallyPerfectGas.update_thermo_from_hs().\n";
	    msg ~= format("The initial temperature guess was: %12.6f\n", Tsave);
	    msg ~= format("The target enthalpy value was: %12.6f\n", h);
	    msg ~= format("The GasState is currently:\n");
	    msg ~= Q.toString();
	    msg ~= "The message from the ridder.solve function is:\n";
	    msg ~= e.msg;
	    throw new Exception(msg);
	}

	// Second, we can iterate to find the pressure that gives
	// correct entropy.
	TOL = 1.0e-3;
	double delp = 1000.0;
	double p1 = Q.p;
	double psave = p1;
	double p2 = p1 + delp;

	auto zeroFun2 = delegate (double p) {
	    Q.p = p;
	    double s_guess = entropy(Q);
	    return s - s_guess;
	};

	if ( bracket!zeroFun2(p1, p2) == -1 ) {
	    string msg = "The 'bracket' function failed to find pressure values\n";
	    msg ~= "that bracketed the zero function in ThermallyPerfectGas.update_thermo_from_hs().\n";
	    msg ~= format("The final values are: p1 = %12.6f and p2 = %12.6f\n", p1, p2);
	    throw new Exception(msg);
	}

	if ( p1 < 0.0 )
	    p1 = 1.0;
	
	try {
	    Q.p = solve!zeroFun2(p1, p2, TOL);
	}
	catch ( Exception e ) {
	    string msg = "There was a problem iterating to find pressure\n";
	    msg ~= "in function ThermallyPerfectGas.update_thermo_from_hs().\n";
	    msg ~= format("The initial pressure guess was: %12.6f\n", psave);
	    msg ~= format("The target entropy value was: %12.6f\n", s);
	    msg ~= format("The GasState is currently:\n");
	    msg ~= Q.toString();
	    msg ~= "The message from the ridder.solve function is:\n";
	    msg ~= e.msg;
	    throw new Exception(msg);
	}
	_tpgMixEOS.update_energy(Q);
	_pgMixEOS.update_density(Q);
    }
    override void update_sound_speed(GasState Q)
    {
	// Reference:
	// Cengel and Boles (1998)
	// Thermodynamics: an Engineering Approach, 3rd edition
	// McGraw Hill
	// Equation 16-10 on p. 849
    
	// "frozen" sound speed
	Q.a = sqrt(gamma(Q)*dpdrho_const_T(Q));
    }
    override void update_trans_coeffs(GasState Q)
    {
	_viscModel.update_viscosity(Q);
	_thermCondModel.update_thermal_conductivity(Q);
    }
    /*
    override void eval_diffusion_coefficients(ref GasState Q) {
	throw new Exception("not implemented");
    }
    */
    override double dudT_const_v(in GasState Q)
    {
	// Noting that Cv = Cp - R
	foreach ( i; 0.._n_species ) _Cv[i] = _curves[i].eval_Cp(Q.Ttr) - _R[i];
	return mass_average(Q, _Cv);
    }
    override double dhdT_const_p(in GasState Q)
    {
	foreach ( i; 0.._n_species ) _Cp[i] = _curves[i].eval_Cp(Q.Ttr);
	return mass_average(Q, _Cp);
    }
    override double dpdrho_const_T(in GasState Q)
    {
	double R = gas_constant(Q);
	return R*Q.Ttr;
    }
    override double gas_constant(in GasState Q)
    {
	return mass_average(Q, _R);
    }
    override double internal_energy(in GasState Q)
    {
	return Q.u;
    }
    override double enthalpy(in GasState Q)
    {
	foreach ( isp; 0.._n_species) {
	    _h[isp] = _curves[isp].eval_h(Q.Ttr);
	}
	return mass_average(Q, _h);
    }
    override double enthalpy(in GasState Q, int isp)
    {
	return _curves[isp].eval_h(Q.Ttr);
    }
    override double entropy(in GasState Q)
    {
	foreach ( isp; 0.._n_species ) {
	    _s[isp] = _curves[isp].eval_s(Q.Ttr) - _R[isp]*log(Q.p/P_atm);
	}
	return mass_average(Q, _s);
    }
    override double entropy(in GasState Q, int isp)
    {
	return _curves[isp].eval_s(Q.Ttr);
    }

private:
    double[] _R;
    PerfectGasMixEOS _pgMixEOS;
    ThermallyPerfectGasMixEOS _tpgMixEOS;
    CEAThermo[] _curves;
    WilkeMixingViscosity _viscModel;
    WilkeMixingThermCond _thermCondModel;
    // Working array space
    double[] _Cp, _Cv, _h, _s;
} // end class Ideal_gas


version(therm_perf_gas_test) {
    int main() {
	import util.msg_service;

	auto gm = new ThermallyPerfectGas("sample-data/therm-perf-5-species-air.lua");
	auto gd = new GasState(5, 0);
	assert(approxEqual(3.621, gm.LJ_sigmas[0]), failedUnitTest());
	assert(approxEqual(97.530, gm.LJ_epsilons[0]), failedUnitTest());
	       
	gd.p = 1.0e6;
	gd.Ttr = 2000.0;
	gd.massf = [0.2, 0.2, 0.2, 0.2, 0.2];
	gm.update_thermo_from_pT(gd);
	assert(approxEqual(11801825.6, gd.u, 1.0e-6), failedUnitTest());
	assert(approxEqual(1.2840117, gd.rho, 1.0e-6), failedUnitTest());

	gd.rho = 2.0;
	gd.u = 14.0e6;
	gm.update_thermo_from_rhou(gd);
	assert(approxEqual(3373757.4, gd.p, 1.0e-6), failedUnitTest());
	assert(approxEqual(4331.944, gd.Ttr, 1.0e-6), failedUnitTest());
    
	gd.Ttr = 10000.0;
	gd.rho = 1.5;
	gm.update_thermo_from_rhoT(gd);
	assert(approxEqual(5841068.3, gd.p, 1.0e-6), failedUnitTest());
	assert(approxEqual(20340105.9, gd.u, 1.0e-6), failedUnitTest());

	gd.rho = 10.0;
	gd.p = 5.0e6;
	gm.update_thermo_from_rhop(gd);
	assert(approxEqual(11164648.5, gd.u, 1.0e-6), failedUnitTest());
	assert(approxEqual(1284.012, gd.Ttr, 1.0e-6), failedUnitTest());

	gd.p = 1.0e6;
	double s = 10000.0;
	gm.update_thermo_from_ps(gd, s);
	assert(approxEqual(2560.118, gd.Ttr, 1.0e-6), failedUnitTest());
	assert(approxEqual(12313952.52, gd.u, 1.0e-6), failedUnitTest());
	assert(approxEqual(1.00309, gd.rho, 1.0e-6), failedUnitTest());

	s = 11000.0;
	double h = 17.0e6;
	gm.update_thermo_from_hs(gd, h, s);
	assert(approxEqual(5273.103, gd.Ttr, 1.0e-6), failedUnitTest());
	assert(approxEqual(14946629.7, gd.u, 1.0e-6), failedUnitTest());
	assert(approxEqual(0.4603513, gd.rho, 1.0e-6), failedUnitTest());
	assert(approxEqual(945271.84, gd.p, 1.0e-4), failedUnitTest());

	gd.Ttr = 4000.0;
	gm.update_trans_coeffs(gd);
	assert(approxEqual(0.00012591, gd.mu, 1.0e-6), failedUnitTest());
	assert(approxEqual(0.2448263, gd.k, 1.0e-6), failedUnitTest());
	
	// [TODO]
	// entropy, enthalpy and sound speed tests
	return 0;
    }
}
