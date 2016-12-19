/**
 * SF6Virial.d
 * SF6 Virial Gas model for use in the CFD codes.
 *
 * Author: Jonathan H.
 * Version: 2015-08-25: initial cut, to explore options.
 */

module gas.sf6virial;

import gas.gas_model;
import gas.physical_constants;
import gas.diffusion.sutherland_viscosity;
import gas.diffusion.sutherland_therm_cond;
import std.math;
import std.stdio;
import std.string;
import std.file;
import std.json;
import std.conv;
import util.lua;
import util.lua_service;
import util.msg_service;
import core.stdc.stdlib : exit;
import nm.ridder;

class SF6Virial:GasModel {
public:
    this() {
	// Default model is mostly initialized in the private data below.
	_n_species = 1;
	_n_modes = 0;
	_species_names ~= "SF6";
	_mol_masses ~= 0.146055;// value for sea-level air
	create_species_reverse_lookup();
    }
	//NOT EXACTLY SURE HOW TO DEAL WITH THIS PART YET -- I am pretty sure it is used for selecting a certain gas model
    this(lua_State *L) {
	this();
	// Bring table to TOS
	lua_getglobal(L, "SF6Virial");
	// Let's overwrite that here.
	_species_names[0] = getString(L, -1, "speciesName");
	// Now, pull out the remaining numeric value parameters.
	_mol_masses ~= getDouble(L, -1, "mMass");
	//_gamma = getDouble(L, -1, "gamma");//not used
	// Reference values for entropy
	lua_getfield(L, -1, "entropyRefValues");
	//_s1 = getDouble(L, -1, "s1");
	//_T1 = getDouble(L, -1, "T1");
	//_p1 = getDouble(L, -1, "p1");
	lua_pop(L, 1);
	// Molecular transport coefficent constants.
	lua_getfield(L, -1, "sutherlandVisc");
	_mu_ref = getDouble(L, -1, "mu_ref");
	_T_mu = getDouble(L, -1, "T_ref");
	_S_mu = getDouble(L, -1, "S");
	lua_pop(L, 1);
	lua_getfield(L, -1, "sutherlandThermCond");
	_k_ref = getDouble(L, -1, "k_ref");
	_T_k = getDouble(L, -1, "T_ref");
	_S_k = getDouble(L, -1, "S");
	lua_pop(L, 1);
	// Compute derived parameters
	_Rgas = R_universal/_mol_masses[0];
	//_Cv = _Rgas / (_gamma - 1.0);
	//_Cp = _Rgas*_gamma/(_gamma - 1.0);
	create_species_reverse_lookup();
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "SF6 =(";
	repr ~= "name=\"" ~ _species_names[0] ~"\"";
	repr ~= ", Mmass=" ~ to!string(_mol_masses[0]);
	//repr ~= ", gamma=" ~ to!string(_gamma);
	//repr ~= ", s1=" ~ to!string(_s1);
	//repr ~= ", T1=" ~ to!string(_T1);
	//repr ~= ", p1=" ~ to!string(_p1);
	repr ~= ", mu_ref=" ~ to!string(_mu_ref);
	repr ~= ", T_mu=" ~ to!string(_T_mu);
	repr ~= ", S_mu=" ~ to!string(_S_mu);
	repr ~= ", k_ref=" ~ to!string(_k_ref);
	repr ~= ", T_mu= " ~ to!string(_T_k);
	repr ~= ", S_k=" ~ to!string(_S_k);
	repr ~= ")";
	return to!string(repr);
    }

    override void update_thermo_from_pT(GasState Q) const 
    {
	Q.rho = updateRho_PT(Q.p, Q.Ttr);
	Q.u = updateEnergy_rhoT(Q.rho, Q.Ttr);
    }
    override void update_thermo_from_rhoe(GasState Q) const
    {
	Q.Ttr = updateTemperature_rhoe(Q.rho, Q.u);
	Q.p = updatePressure_rhoT(Q.rho,Q.Ttr);
    }
    override void update_thermo_from_rhoT(GasState Q) const//DONE
    {
	Q.p = updatePressure_rhoT(Q.rho, Q.Ttr);
	Q.u = updateEnergy_rhoT(Q.rho, Q.u);
    }
    override void update_thermo_from_rhop(GasState Q) const
    {
	Q.u = updateEnergy_Prho(Q.p, Q.rho);//might want to fix the order that this solves in
	Q.Ttr = updateTemperature_rhoe(Q.rho, Q.u);
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
	Q.a = updateSoundSpeed_rhoT(Q.rho, Q.Ttr);
    }
    override void update_trans_coeffs(GasState Q) const
    {
	Q.mu = sutherland_viscosity(Q.Ttr, _T_mu, _mu_ref, _S_mu);
	Q.k = sutherland_thermal_conductivity(Q.Ttr, _T_k, _k_ref, _S_k);
    }
    /*
    override void eval_diffusion_coefficients(ref GasState Q) {
	throw new Exception("not implemented");
    }
    */
    override double dedT_const_v(in GasState Q) const
    {
	return get_de_dT(Q.rho,Q.Ttr);
    }
    override double dhdT_const_p(in GasState Q) const
    {
	throw new Exception(format("Not implemented: line=%d, file=%s\n", __LINE__, __FILE__));
    }
    override double dpdrho_const_T(in GasState Q) const
    {
	double R = gas_constant(Q);
	return R*Q.Ttr;
    }
    override double gas_constant(in GasState Q) const
    {
	return R_universal/_mol_masses[0];
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
	throw new Exception(format("Not implemented: line=%d, file=%s\n", __LINE__, __FILE__));
    }

private:
    // Thermodynamic constants
	double _Rgas = 8.3144621/0.146055;
	double[6] _a = [0, 0, -49.9051433, 4.124606e-2, -1.612953e-5, -4.899779e-11];
	double[6] _b = [0, 0, 5.485495e-2, -3.340088e-5, 0, 1.094195e-11];
	double[6] _c = [0, 0, -2.375924505e3, 2.819595, 0, -3.08231e-7];
	double _k = 6.883022;
	double _Tc = 318.8;
	double _d = 3.27367367e-4;
	double[6] _G = [0.0, -107.9122479, 3.94226447, -5.128665e-3, 2.422895e-6, -9.6020764e5];
	double _u0 = 50000;//dummy values
	double _T0 = 300;//dummy values
    // Molecular transport coefficent constants.
    double _mu_ref = 1.716e-5; // Pa.s
    double _T_mu = 273.0; // degrees K
    double _S_mu = 111.0; // degrees K
    double _k_ref = 0.0241; // W/(m.K) 
    double _T_k = 273.0; // degrees K
    double _S_k = 194.0; // degrees K
    
    //define some functions that will be available to the functions in the private area    
    const double updatePressure_rhoT(double rho, double T){
	double v = 1/rho; //equation is in specific volume
	double sum = 0;
	for(int i = 2; i != 6; i++) sum += (_a[i] + _b[i]*T + _c[i]*exp(-_k*T/_Tc))/(v - _d)^^i;
	return _Rgas*T/(v - _d) + sum;
	
   }

   const double updateEnergy_rhoT(double rho, double T){
	//From 1995 paper Kyle Anderson & Dimitri Mavriplis
	double v = 1/rho;
	double integralcv0 = (_G[1]-_Rgas)*(T - _T0) + 0.5*_G[2]*(T^^2 - _T0^^2) + 1.0/3.0*_G[3]*(T^^3 - _T0^^3)
						+ 0.25*_G[4]*(T^^4 - _T0^^4) - _G[5]*(1.0/T - 1.0/_T0);
	double integralDensity = 0;
	for(int i = 2; i != 6; i++) integralDensity += (_a[i] + (1 + _k*T/_Tc)*_c[i]*exp(-_k*T/_Tc))/(i - 1)/(v - _d)^^(i - 1);
	return integralcv0 + integralDensity + _u0;
   }
   const double updateTemperature_rhoe(double rho, double e, int maxIterations = 100, double Ttol = 0.1){
	double T = 400; // first approximation using totally ideal gas not possible because we don't know pressure;
	for(int i = 0; i != maxIterations; i++){
		double deltaT = (updateEnergy_rhoT(rho,T)-e)/get_de_dT(rho,T);
		/* writeln("deltaT: ", deltaT);
		writeln("Energy(T): ", getSpecificEnergy(rho,T,gas)-e);
		writeln("dE/dT: ", get_de_dT(rho,T,gas)); */
		if (abs(deltaT) < Ttol){
			//writeln("tolerance for T met");
			break;
		}
		T -= deltaT;
		if (i == maxIterations) throw new Exception(format("Newton-Cotes reached max iterations when calculating T from rho = %s, e = %s", rho, e));
		}
	assert(!isNaN(T), format("Newton-Cotes failed when calculating T from rho = %s, e = %s", rho, e));
	return T;
   }
   const double get_de_dT(double rho, double T) {
	//Gets derivative of specific energy w.r.t. Temperature for evaluating T as a function of e, rho based on Newton's method
	//conveniently de_dT is c_V
	//From 1995 paper Kyle Anderson & Dimitri Mavriplis
	double v = 1.0/rho;
	double cv0 = _G[1] - _Rgas + _G[2]*T + _G[3]*T^^2 + _G[4]*T^^3 + _G[5]/T^^2;
	double departureFunction = 0;
	for(int i = 2; i != 6; i++) departureFunction += _c[i]/(i-1)/(v-_d)^^(i-1)*-(_k/_Tc)^^2*T*exp(-_k*T/_Tc);
	return cv0 + departureFunction;
	
   }
   const double updateEnergy_Prho(double P, double rho, double[2] Tbracket = [300, 1300], double tol = 0.001){
	//double e = min(max(getSpecificEnergy(rho,P/rho/_Rgas,Gas),_u0*0.6),2.5e6);//estimate energy based on temperature, don't go too high or too low
	//writeln("Temperature guess", P/rho/_Rgas);
	//writeln("first guess at e: ", e);
	auto getPressure_T = delegate (double T){return updatePressure_rhoT(rho, T) - P;};
	double[2] e_bracket = [updateEnergy_rhoT(rho, Tbracket[0]), updateEnergy_rhoT(rho, Tbracket[1])];
	double T = solve!getPressure_T(Tbracket[0], Tbracket[1],tol);
	double e = updateEnergy_rhoT(rho, T);
	string errorString = "P: " ~ to!string(P) ~ ", rho: " ~ to!string(rho);
	assert(!isNaN(e), errorString);
	
	return e;}
   const double updateRho_PT(double P, double T, double[2] bracket = [2, 1000], double tol = 0.1){
	//when temperature is close to 	critical temperature be very careful
	//assert(T > 305, "Temperature too low and close to critical point (305K)");
	auto getPressure_rho = delegate (double rho){return updatePressure_rhoT(rho, T) - P;};
	double rho = solve!getPressure_rho(bracket[0], bracket[1], tol);
	string errorString = "P: " ~ to!string(P) ~ ", T: " ~ to!string(T);
	assert(!isNaN(rho), errorString);
	return rho;
   }
   const double updatePressure_rhoe(double rho, double e){
   	   double T = updateTemperature_rhoe(rho, e);
  	 return updatePressure_rhoT(rho, T);}
   
   const double updatePressure_Trho(double T, double rho){
   return updatePressure_rhoT(rho, T);}//just a re-casting of the same function
      
   const double updateSoundSpeed_rhoT(double rho, double T) {
	double v = 1.0/rho;
	double dP_dv = -_Rgas*T/(v-_d)^^2;
	double dP_dT = _Rgas/(v - _d);
	for(int i = 2; i != 6; i++) {
		dP_dv += -i*(v-_d)^^(-i - 1)*(_a[i] + _b[i]*T + _c[i]*exp(-_k*T/_Tc));
		dP_dT += (_b[i] + _c[i]*_k/_Tc*exp(-_k*T/_Tc))/(v-_d)^^i;
		}
	double dP_drho = -dP_dv/rho/rho;
	double cV = get_de_dT(rho, T);
	//writefln("rho: %f, T: %f, dP_drho: %f, dP_dT: %f", rho, T, dP_drho, dP_dT);
	assert(dP_drho + dP_dT ^^2 * T/rho^^2/cV > 0, 
		format("Tried to square root negative number while calculating sound speed, rho = %s, T = %s",rho, T));
	return  sqrt(dP_drho + dP_dT ^^2 * T/rho^^2/cV);
	}


} // end class SF6

unittest {//need to write these properly
    import std.stdio;
    auto gm = new IdealGas();
    assert(gm.species_name(0) == "ideal air", "species name list");
    auto gd = new GasState(gm, 100.0e3, 300.0);
    assert(approxEqual(gm.R(gd), 287.086), failedUnitTest());
    assert(gm.n_modes == 0, failedUnitTest());
    assert(gm.n_species == 1, failedUnitTest());
    assert(approxEqual(gd.p, 1.0e5), failedUnitTest());
    assert(approxEqual(gd.Ttr, 300.0), failedUnitTest());
    assert(approxEqual(gd.massf[0], 1.0), failedUnitTest());

    gm.update_thermo_from_pT(gd);
    gm.update_sound_speed(gd);
    assert(approxEqual(gd.rho, 1.16109), failedUnitTest());
    assert(approxEqual(gd.u, 215314.0), failedUnitTest());
    assert(approxEqual(gd.a, 347.241), failedUnitTest());
    gm.update_trans_coeffs(gd);
    assert(approxEqual(gd.mu, 1.84691e-05), failedUnitTest());
    assert(approxEqual(gd.k, 0.0262449), failedUnitTest());

    lua_State* L = init_lua_State("sample-data/ideal-air-gas-model.lua");
    gm = new IdealGas(L);
    lua_close(L);
    assert(approxEqual(gm.R(gd), 287.086), failedUnitTest());
    assert(gm.n_modes == 0, failedUnitTest());
    assert(gm.n_species == 1, failedUnitTest());
    assert(approxEqual(gd.p, 1.0e5), failedUnitTest());
    assert(approxEqual(gd.Ttr, 300.0), failedUnitTest());
    assert(approxEqual(gd.massf[0], 1.0), failedUnitTest());

    gm.update_thermo_from_pT(gd);
    gm.update_sound_speed(gd);
    assert(approxEqual(gd.rho, 1.16109), failedUnitTest());
    assert(approxEqual(gd.u, 215314.0), failedUnitTest());
    assert(approxEqual(gd.a, 347.241), failedUnitTest());
    gm.update_trans_coeffs(gd);
    assert(approxEqual(gd.mu, 1.84691e-05), failedUnitTest());
    assert(approxEqual(gd.k, 0.0262449), failedUnitTest());
}
