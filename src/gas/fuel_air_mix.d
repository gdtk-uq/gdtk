/**
 * fuel_air_mix.d
 *
 * Fuel+Air->Products reacting gas as used by JJ.
 *
 * Authors: JJ Hoste, Peter J. and Rowan G.
 * Version: 2017-04-22: initial shell copied from powers-aslam-gas.
 */

module gas.fuel_air_mix;

import gas.gas_model;
import gas.physical_constants;
import gas.therm_perf_gas;
import gas.thermo.cea_thermo_curves;
import gas.thermo.perf_gas_mix_eos;
import gas.thermo.therm_perf_gas_mix_eos;
import gas.diffusion.viscosity;
import gas.diffusion.therm_cond;
import gas.diffusion.cea_viscosity;
import gas.diffusion.cea_therm_cond;
import gas.diffusion.wilke_mixing_viscosity;
import gas.diffusion.wilke_mixing_therm_cond;
import std.math;
import std.stdio;
import std.string;
import std.file;
import std.json;
import std.conv;
import std.algorithm;
import util.lua;
import util.lua_service;
import core.stdc.stdlib : exit;
import nm.bracketing;
import nm.brent; 
import kinetics.reaction_mechanism;
import kinetics.reaction;

// First, the basic gas model.

class FuelAirMix: ThermallyPerfectGas {
public:
    this(lua_State *L) {
	super(L); 
	lua_getglobal(L, "FuelAirMix"); 
	// Now, pull out the numeric value parameters.
	_A_edm = getDouble(L, -1, "Aedm");
	_B_edm = getDouble(L, -1, "Bedm");
	_laminar_limit = getBool(L, -1, "lamLimit");
	lua_pop(L, 1); // dispose of the table
    } // end constructor

    override string toString() const
    {
	char[] repr;
	repr ~= "FuelAirMix =(";
	repr ~= " Aedm=" ~ to!string(_A_edm);
	repr ~= ", Bedm=" ~ to!string(_B_edm);
	repr ~= ")";
	return to!string(repr);
    }

private:
    //settings specific to EDM model
    double _A_edm;
    double _B_edm;
    bool _laminar_limit;	
} // end class FuelAirMix


// Now, for the reaction update...
// It is included here because it is a small amount of code and
// it is specific to this particular gas model.

final class MixingLimitedUpdate : ThermochemicalReactor {
    ReactionMechanism rmech; //only used if we desire to limit the reaction rate
    this(string fname, GasModel gmodel)
    {
	super(gmodel); // hang on to a reference to the gas model
	// We need to pick a number of pieces out of the gas-model file, again.
	// Although they exist in the GasModel object, they are private.
	//lua_State *L;
	auto L = init_lua_State();
	
	// get details about type of species: fuel, ox, prod or inert
	doLuaFile(L, fname); //reading the reaction file here
	_n_species = gmodel.n_species;
	lua_getglobal(L, "speciesType");
	getArrayOfStrings(L, LUA_GLOBALSINDEX, "speciesType", _species_type);
	lua_pop(L, 1); // dispose of the table
	foreach(isp;0.._n_species){
	    if(_species_type[isp]=="fuel"){_n_fuel=_n_fuel+1;}
	    if(_species_type[isp]=="ox"){_n_ox=_n_ox+1;}
	    if(_species_type[isp]=="prod"){_n_prod=_n_prod+1;}
	}
	_n_reacting=_n_fuel+_n_ox+_n_prod;
	writeln("nb of fuel species = ",_n_fuel);
	writeln("nb of ox species = ",_n_ox);
	writeln("nb of prod species = ",_n_prod);

	// get settings for EDM model
	lua_getglobal(L, "FuelAirMix");
	// Now, pull out the numeric value parameters.
	_A_edm = getDouble(L, -1, "Aedm");
	_B_edm = getDouble(L, -1, "Bedm");
	_laminar_limit = getBool(L, -1, "lamLimit");
	lua_pop(L, 1); // dispose of the table
	//writeln(_laminar_limit);

	// if true we want to limit the reaction rate with the "no-model" one
	// in that case we initialize the ReactionMechanism object
	if (_laminar_limit) {
	    lua_getglobal(L, "config");
	    lua_getfield(L, -1, "tempLimits");
	    lua_getfield(L, -1, "lower");
	    double T_lower_limit = lua_tonumber(L, -1);
	    lua_pop(L, 1);
	    lua_getfield(L, -1, "upper");
	    double T_upper_limit = lua_tonumber(L, -1);
	    lua_pop(L, 1);
	    lua_pop(L, 1);
	    lua_pop(L, 1);
	    lua_getglobal(L, "reaction");	
	    rmech = createReactionMechanism(L, gmodel, T_lower_limit, T_upper_limit);
	    lua_rawgeti(L, -1, 1); //only 1 reaction called reaction[1]
	    lua_pop(L, 1);
	}

	lua_getglobal(L, "reaction");	
	lua_rawgeti(L, -1, 1); //only 1 reaction called reaction[1]
	int[] reacCoeffs; // required to set stochiometric mass ratio
	getArrayOfInts(L, -1, "reacCoeffs", reacCoeffs);
	writeln(reacCoeffs,reacCoeffs[0]);
	if (_n_reacting > 3) {
	    int[] prodCoeffs;
	    _nu_W.length=_n_reacting;
	    getArrayOfInts(L, -1, "prodCoeffs", prodCoeffs);
	    foreach(isp; 0 .. _n_fuel+_n_ox) {
		writeln(isp);
		_nu_W[isp]= 1.0 * reacCoeffs[isp]*gmodel.mol_masses[isp]; 
	    }
	    foreach(isp; 0 .. _n_prod) {
		_nu_W[_n_fuel+_n_ox+isp]= -1.0 * prodCoeffs[isp]*gmodel.mol_masses[isp+_n_fuel+_n_ox]; 
	    }
	    lua_pop(L, 1);
	} // else we don't need _nu_W in further computations
	lua_pop(L, 1); // I believe needed to get out 1 step back to all reactions
	lua_close(L);
	writeln("stochiometric mass ratio ",_sRatio);
	set_stochiometricRatio(reacCoeffs,gmodel.mol_masses);
	writeln("stochiometric mass ratio ",_sRatio);
    } // end of constructor


    void set_stochiometricRatio(int[] stochCoeff, double[] mol_masses)
    {	
	double num=0.0,denom=0.0;	
	foreach (isp; 0 .. _n_fuel) {
	    denom = denom + mol_masses[isp]*stochCoeff[isp];
	}
	_nuF_WF = denom; //only needed when more than 3 reacting species/prod

	foreach (isp; 0 .. _n_ox) {
	    num = num + mol_masses[_n_fuel+isp]*stochCoeff[_n_fuel+isp];
	}
	_sRatio=num/denom;
    }

    void get_massf(GasState gas, ref double[3] Ys)
    {
	// pass Ys by ref required to fill local variable
	// initialization required in case of cumulative sum below
	Ys=[0.0,0.0,0.0];
	int i=0;
	if (_n_fuel > 1) {
	    foreach (isp; 0 .. _n_fuel) {
		Ys[0]=Ys[0]+gas.massf[isp];
		i=_n_fuel-1;
	    }
	} else {
	    Ys[0]=gas.massf[0];
	    writeln(" massf sp0 ",Ys[0]);
	}
	if (_n_ox>1) {
	    foreach(isp;0.._n_ox){
		Ys[1]=Ys[1]+gas.massf[_n_fuel+isp];
	    }
	} else {
	    Ys[1]=gas.massf[_n_fuel];
	    //writeln(" massf sp1 ",Ys[1]);
	}
	if(_n_prod>1) {
	    foreach(isp;0.._n_prod){
		Ys[2]=Ys[2]+gas.massf[_n_fuel+_n_ox+isp];
		//writeln(" massf prod ",Ys[2]);
	    }
	} else {
	    Ys[2]=gas.massf[_n_fuel+_n_ox];
	}
    }

    double get_omegaDotF(double omega, double[3] Ys)
    {
	return -1.0 * _A_edm / _tau_edm * min(Ys[0],Ys[1]/_sRatio,_B_edm*Ys[2]/(1.0+_sRatio));
    }

    override void opCall(GasState Q, double tInterval, ref double dtSuggest,ref double[] params)
    {
	double local_dt_global=tInterval;
	double omega = params[0];
	_tau_edm=1.0/(_beta_star*omega);
	double _omega_dot_F = 0.0;
	double omega_dot_O = 0.0;
	double omega_dot_P = 0.0;
	double[] rates; // for laminar= "no-model" combustion
	rates.length = _n_species;
	//writeln(rmech.) 
	double[3] _Y_s;
	if(_n_reacting > 3){
	    get_massf(Q,_Y_s);
	    //writeln("species massf ",_Y_s);
	    _omega_dot_F = get_omegaDotF(omega,_Y_s);
	    //writeln("fuel RR ",_omega_dot_F);
	    //writeln("nu F x W F = ",_nuF_WF);
	    //writeln("nu  x W  = ",_nu_W);
	    if(_laminar_limit){
		double[] conc;
		conc.length=_n_species; 
		rmech.eval_rate_constants(Q);  
		_gmodel.massf2conc(Q,conc);
		writeln("concentration ",conc);
		rmech.eval_rates(conc,rates); 
		writeln("laminar RR (1/s)= ",rates);
		double[] lamRR; //in kg/(m^3 . s)
		lamRR.length=rates.length;
		double sumLamRR=0.0;
		foreach (isp; 0 .. _n_species) {
		    lamRR[isp]= _gmodel.mol_masses[isp]/Q.rho* rates[isp];
		    sumLamRR=sumLamRR+lamRR[isp];
		}
		writeln("laminar RR (1/s) = ",lamRR);
		writeln("sum of laminar RR (1/s) = ",sumLamRR);
		writeln("Fuel RR (1/s) before lam eval= ",_omega_dot_F);
		_omega_dot_F=-1.0*(min(abs(_omega_dot_F),abs(rates[0]))); 
	    }
	    writeln("Fuel RR (1/s)= ",_omega_dot_F);
	    foreach (isp; 0 .. _n_reacting) { //don't account for inert species
		Q.massf[isp]=Q.massf[isp] + local_dt_global * (_nu_W[isp])/ _nuF_WF* _omega_dot_F;
	    }
	    //double sum=0.0;
	    //foreach(isp;0.._n_species){
		//sum=sum+Q.massf[isp];					
	    //}
	    //writeln("sum of species = ",sum);
	} else {
	    _Y_s[0]=Q.massf[0];
	    _Y_s[1]=Q.massf[1];
	    _Y_s[2]=Q.massf[2];
	    _omega_dot_F = get_omegaDotF(omega,_Y_s);	// units 1/s in here 
	    if(_laminar_limit) {
		double[] conc;
		conc.length=_n_species; 
		rmech.eval_rate_constants(Q);  
		_gmodel.massf2conc(Q,conc);
		writeln("concentration ",conc);
		rmech.eval_rates(conc,rates); 
		writeln("laminar RR (1/s)= ",rates);
		double[] lamRR; //in kg/(m^3 . s)
		lamRR.length=rates.length;
		double sumLamRR=0.0;
		foreach (isp; 0 .. _n_species) {
		    lamRR[isp]= _gmodel.mol_masses[isp]/Q.rho* rates[isp];
		    sumLamRR=sumLamRR+lamRR[isp];
		}
		writeln("laminar RR (1/s) = ",lamRR);
		writeln("sum of laminar RR (1/s) = ",sumLamRR);
		writeln("Fuel RR (1/s) before lam eval= ",_omega_dot_F);
		_omega_dot_F=-1.0*(min(abs(_omega_dot_F),abs(rates[0]))); // units are same 
		// NOTE: when fuel/ox almost zero the laminar RR will give positive values for
		//       these reaction rates BUT the EDM one is much smaller and dominates	
		//		 only reason for laminar is to avoid reaction where T too low and not enough
		//		 mixing occurs.			
	    } // end if /else
	    omega_dot_O = _sRatio * _omega_dot_F;
	    omega_dot_P = - (1.0+_sRatio) *_omega_dot_F;
	    //writeln("Fuel RR (1/s)= ",_omega_dot_F);
	    //writeln(omega_dot_O);
	    //writeln(omega_dot_P);
	    Q.massf[0]=_Y_s[0] + local_dt_global * _omega_dot_F;
	    // following is actuall still correct for laminar RR
	    Q.massf[1]=_Y_s[1] + local_dt_global * _sRatio*_omega_dot_F;
	    Q.massf[2]=_Y_s[2] + local_dt_global * - (1.0+_sRatio) *_omega_dot_F;		
	} //end if/ else
    } //end function opCall

private:
    // EDM model constants
    double _A_edm=4.0; 
    double _B_edm=0.5;
    double _tau_edm;
    // boolean for use of laminar reaction rate limit
    bool _laminar_limit=false;
    static immutable _beta_star=0.09;
    string[] _species_type; // fuel, ox, prod, inert

    int _n_fuel, _n_ox, _n_prod, _n_species, _n_reacting;
    double _sRatio; // mass stochiometric ratio
    // parameters used if nb of species in reaction (_n_reacting) > 3
    double _nuF_WF; // stochiometric coeff fuel x molar mass fuel
    double[] _nu_W; // the above for every species
} // end class MixingLimitedUpdate


// Unit test of the basic gas model...

version(fuel_air_mix_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
	lua_State* L = init_lua_State();
	doLuaFile(L, "sample-data/fuel-air-mix-model.lua");
	auto gm = new FuelAirMix(L);
	lua_close(L);
	auto gd = new GasState(2, 0);
	gd.p = 1.0e5;
	gd.Ttr = 300.0;
	gd.massf[0] = 0.75; gd.massf[1] = 0.25;
	/+
	 [FIX-ME]
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
	+/
	return 0;
    }
}
