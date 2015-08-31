/**
 * CO2Gas.d
 * Bender gas model for use in the CFD codes.
 *
 * Author: Jonathan H.
 * Version: 2015-07-27: initial cut, to explore options.
 */

module gas.CO2Gas;

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
import std.c.stdlib : exit;

class CO2Gas: GasModel {
public:
    this() {
	// Default model is mostly initialized in the private data below.
	_n_species = 1;
	_n_modes = 1;
	_species_names ~= "CO2";
	_mol_masses ~= 0.04401121121333065;// value for sea-level air
    }
	//NOT EXACTLY SURE HOW TO DEAL WITH THIS PART YET -- I am pretty sure it is used for selecting a certain gas model
    this(lua_State *L) {
	this();
	// Bring table to TOS
	lua_getglobal(L, "CO2Gas");
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
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "CO2 =(";
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
	assert(Q.T.length == 1, "incorrect length of temperature array");
	Q.rho = updateRho_PT(Q.p, Q.T[0]);
	Q.e[0] = updateEnergy_rhoT(Q.rho, Q.T[0]);
    }
    override void update_thermo_from_rhoe(GasState Q) const
    {
	assert(Q.e.length == 1, "incorrect length of energy array");
	Q.T[0] = updateTemperature_rhoe(Q.rho, Q.e[0]);
	Q.p = updatePressure_rhoT(Q.rho,Q.T[0]);
    }
    override void update_thermo_from_rhoT(GasState Q) const//DONE
    {
	assert(Q.T.length == 1, "incorrect length of temperature array");
	//Calculating Pressure	
	Q.p = updatePressure_rhoT(Q.rho, Q.T[0]);
	//Calculate Energy
	Q.e[0] = updateEnergy_rhoT(Q.rho, Q.e[0]);

    }


    override void update_thermo_from_rhop(GasState Q) const
    {
	assert(Q.T.length == 1, "incorrect length of temperature array");
	Q.e[0] = updateEnergy_Prho(Q.p, Q.rho);//might want to fix the order that this solves in
	Q.T[0] = updateTemperature_rhoe(Q.rho, Q.e[0]);
	
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
	Q.a = updateSoundSpeed_rhoT(Q.rho, Q.T[0]);
    }
    override void update_trans_coeffs(GasState Q) const
    {
	Q.mu = sutherland_viscosity(Q.T[0], _T_mu, _mu_ref, _S_mu);
	Q.k[0] = sutherland_thermal_conductivity(Q.T[0], _T_k, _k_ref, _S_k);
    }
    /*
    override void eval_diffusion_coefficients(ref GasState Q) {
	throw new Exception("not implemented");
    }
    */
    override double dedT_const_v(in GasState Q) const
    {
	return get_de_dT(Q.rho,Q.T[0]);
    }
    override double dhdT_const_p(in GasState Q) const
    {
	throw new Exception(format("Not implemented: line=%d, file=%s\n", __LINE__, __FILE__));
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
	throw new Exception(format("Not implemented: line=%d, file=%s\n", __LINE__, __FILE__));
    }

private:
    // Thermodynamic constants
	double _Rgas = 188.918;
	double _u0 = 3.2174105e5;
	double _s0 = 2.1396056e2;
	double _M = 44.01;//kg/kmol
	double _Tc = 304.21;//K
	double _Pc = 7.3834e6;//Pc
	double _rhoc = 464.00;//kg/m^304
	double _T0 = 216.54;//
	double[] _As = [0.0, 2.2488558e-1, -1.3717965e2, -1.4430214e4, -2.9630491e6, -2.0606039e8, 4.5554393e-5, 7.7042840e-2, 4.0602371e1, 
		4.0029509e-7, -3.9436077e-4, 1.2115286e-10, 1.0783386e-7, 4.3962336e-11, -3.6505545e4, 1.9490511e7, -2.9186718e9,
		2.4358627e-2, -3.7546530e1, 1.1898141e4, 5.0e-6];
	double[] _G = [0.0, 8.726361e3, 1.840040e2, 1.914025, -1.667825e-3, 7.305950e-7, -1.255290e-10];
    // Molecular transport coefficent constants.
    double _mu_ref = 1.716e-5; // Pa.s
    double _T_mu = 273.0; // degrees K
    double _S_mu = 111.0; // degrees K
    double _k_ref = 0.0241; // W/(m.K) 
    double _T_k = 273.0; // degrees K
    double _S_k = 194.0; // degrees K
    
    //define some functions that will be available to the functions in the private area    
    const double updatePressure_rhoT(double rho, double T){
	double gamma = _As[20];
	double sum1_5 = 0;
	for (int i = 1; i!=6; i++){sum1_5+= _As[i]*T^^(2-i);}
	double sum6_8 = 0;
	for (int i = 6; i!=9; i++){sum6_8+= _As[i]*T^^(7-i);}
	double sum14_16 = 0;
	for (int i = 14; i!=17; i++){sum14_16+= _As[i]*T^^(12-i);}
	double sum17_19 = 0;
	for (int i = 17; i!=20; i++){sum17_19+= _As[i]*T^^(15-i);}
	
	return rho*_Rgas*T + rho^^2*sum1_5 + rho^^3 * sum6_8 +
				rho^^4*(_As[9]*T + _As[10]) + rho ^^5 *(_As[11]*T + _As[12]) + rho^^6*_As[13] + 
				(rho^^3*sum14_16 + rho^^5*sum17_19)*exp(-gamma*rho^^2);
   }

   const double updateEnergy_rhoT(double rho, double T){
	double gamma = _As[20];
   	double integral_cV_dT = _G[1]*log(T/_T0) + _G[2]*(T-_T0) + _G[3]*(T^^2 - _T0^^2)/2.0 
			+ _G[4]*(T^^3 - _T0^^3)/3.0 + _G[5]*(T^^4 - _T0^^4)/4 + _G[6]*(T^^5 - _T0^^5)/5;
	double integralDensityChange = rho*(_As[2] + 2*_As[3]*T^^-1 + 3*_As[4]*T^^-2 + 4*_As[5]*T^^-3)
		+ rho^^2/2*(_As[7] + 2*_As[8]*T^^-1)
		+ rho^^3/3*_As[10]
		+ rho^^4/4*_As[12]
		+ rho^^5/5*_As[13]
		+ exp(-gamma*rho^^2)*(-(2*gamma)^^-1*(3*_As[14]*T^^-2 + 4*_As[15]*T^^-3 + 5*_As[16]*T^^-4)
							  + (-rho^^2/(2*gamma) - (2*gamma^^2)^^-1)*(3*_As[17]*T^^-2+4*_As[18]*T^^-3 + 5*_As[19]*T^^-4));
	return _u0 + integral_cV_dT + integralDensityChange;
   }
   const double updateTemperature_rhoe(double rho, double e, int maxIterations = 100, double Ttol = 0.1){
	assert(e >= _u0*0.3, "energy below 10% of reference");//make sure you are above the reference point
	double T = 1000*(e - 441482)/(1.42035e6-411482)+400; // first approximation using totally ideal gas not possible because we don't know pressure;
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
		}
	return T;
   }
   const double get_de_dT(double rho, double T) {
	//Gets derivative of specific energy w.r.t. Temperature for evaluating T as a function of e, rho based on Newton's method
	//Polynomials were derived by hand, not sure if accurate yet
	//coincidenty de_dt = c_V
	double gamma = _As[20];
	return rho*(-2*_As[3]*T^^-2 -6*_As[4]*T^^-3 - 12*_As[5]*T^^-4)
			+ rho^^2/2.0*(-2*_As[8]*T^^-2)
			+ exp(-gamma*rho^^2)*(-(2*gamma)^^-1*(-6*_As[14]*T^^-3 - 12*_As[15]*T^^-4 - 20*_As[16]*T^^-5)
								 + (-rho^^2/(2*gamma) - (2*gamma^^2)^^-1)*(6*_As[17]*T^^-3 - 12*_As[18]*T^^-4 - 20*_As[19]*T^^-5))
			+ _G[1]/T + _G[2] + _G[3]*T + _G[4]*T^^2 + _G[5]*T^^3 + _G[6]*T^^4;
	
   }
   const double updateEnergy_Prho(double P, double rho, double[2] Tbracket = [300, 1300], double tol = 0.001){
	//will get Temperature based on Riddler's Method - see PJ's notes
	//double e = min(max(getSpecificEnergy(rho,P/rho/Gas.R,Gas),Gas.u0*0.6),2.5e6);//estimate energy based on temperature, don't go too high or too low
	//writeln("Temperature guess", P/rho/Gas.R);
	//writeln("first guess at e: ", e);
	double[2] bracket = [updateEnergy_rhoT(rho, Tbracket[0]),updateEnergy_rhoT(rho, Tbracket[1])];
	double e = RidderSolverXYZ(P,&updatePressure_rhoe, rho, bracket, tol);
	string errorString = "P: " ~ to!string(P) ~ ", rho: " ~ to!string(rho);
	assert(!isNaN(e), errorString);
	//WORKS WELL FOR SUPERCRITICAL
   return e;}
   const double updateRho_PT(double P, double T, double[2] bracket = [2, 1000], double tol = 0.1){
	//when temperature is close to 	critical temperature be very careful
	assert(T > 305, "Temperature too low and close to critical point (305K)");
	double rho = RidderSolverXYZ(P, &updatePressure_Trho, T, bracket, tol);
	string errorString = "P: " ~ to!string(P) ~ ", T: " ~ to!string(T);
	assert(!isNaN(rho), errorString);
	return rho;
   }
   const double updatePressure_rhoe(double rho, double e){
   	   double T = updateTemperature_rhoe(rho, e);
   return updatePressure_rhoT(rho, T);}
   
   const double updatePressure_Trho(double T, double rho){
   return updatePressure_rhoT(rho, T);}//just a re-casting of the same function
      
   alias simpleXYZFunction = const double delegate(double, double);
   const double RidderSolverXYZ(double z, simpleXYZFunction f, double y, double[2] bracket, double tol){//possibly could make a new function (look at therm_perf_gas.d) LINE220
	//solves the function z = f(y,x) for x given a certain z and y 
	//writeln("The bracket is: ", bracket);
	//writeln("rho: ", rho);
	//order of function is rho, T, Gas
	double f1 = f(y, bracket[0]) - z;
	double f2 = f(y, bracket[1]) - z;
	int i = 0;
	double x4;
	while ((bracket[1] - bracket[0])>tol){
		if (f2 == 0){return bracket[1];}//this is to prevent errors when f2 = 0
		if (f1 == 0){return bracket[0];}
		double x3 = (bracket[0] + bracket[1])/2.0;//x3
		
		double f3 = f(y, x3) - z;
		
		double epsilon = (f3 + sgn(f2)*sqrt(f3^^2 - f1*f2))/(f2+1e-12);//had to add 1e-12 to stop the divide by zero error
		//writeln("epsilon: ", epsilon);
		x4 = x3 - f3*epsilon*(bracket[0] - x3)/(f1-epsilon*f3);
		double f4 = f(y,x4) - z;
		//writeln("x1: ", bracket[0], " x2: ", bracket[1], " x3: ", x3, " x4:", x4);
		//writefln("f1: %e f2: %e f3: %e f4: %e",f1, f2*1e40, f3, f4);
		if (f3*f2 < 0.0){
			if (f4*f2 < 0.0) {
				
				bracket[0] = x4; f1 = f4;
			}
			else {
				
				bracket[0] = x3; f1 = f3;
				bracket[1] = x4; f2 = f4;
			}
			}
		else{
			if (f4*f1 < 0.0) {
				
				bracket[1] = x4; f2 = f4;
			}
			else{
				
				bracket[0] = x4; f1 = f4;
				bracket[1] = x3; f2 = f3;}
			}
		i++;
		
		assert(i != 100);
		//writeln("x3 is: ", x3);
		//writeln("x4 is: ", x4);
		//writeln("bracket is: ", bracket);
		
	}
	return x4;
	}
	const double updateSoundSpeed_rhoT(double rho, double T) {
		double gamma = _As[20];
		double dP_drho = _Rgas*T + 2*rho*(_As[1]*T + _As[2] + _As[3]*T^^-1 + _As[4]*T^^-2 + _As[5]*T^^-3) 
						+ 3*rho^^2*(_As[6]*T + _As[7] + _As[8]*T^^-1)
						+ 4*rho^^3*(_As[9]*T + _As[10])
						+ 5*rho^^4*(_As[11]*T + _As[12])
						+ 6*rho^^5*_As[13]
						+ exp(-gamma*rho^^2) * (-2*gamma*rho*(rho^^3*(_As[14]*T^^-2 + _As[15]*T^^-3 + _As[16]*T^^-4)
															 +rho^^5*(_As[17]*T^^-2 + _As[18]*T^^-3 + _As[19]*T^^-4))
												+ 3*rho^^2*(_As[14]*T^^-2 + _As[15]*T^^-3 + _As[16]*T^^-4)
												+ 5*rho^^4*(_As[17]*T^^-2 + _As[18]*T^^-3 + _As[19]*T^^-4));	
		double dP_dT = rho*_Rgas + rho^^2*(_As[1] - _As[3]*T^^-2 - 2*_As[4]*T^^-3 - 3*_As[5]*T^^-4)
					 + rho^^3*(_As[6] - _As[8]*T^^-2)
					 + rho^^4*_As[9]
					 + rho^^5*_As[11]
					 + exp(-gamma*rho^^2)*(rho^^3*(-2*_As[14]*T^^-3 - 3*_As[15]*T^^-4 - 4*_As[16]*T^^-5)
										 + rho^^5*(-2*_As[17]*T^^-3 - 3*_As[18]*T^^-4 - 4*_As[19]*T^^-5));
		double cV = get_de_dT(rho, T);
		return  sqrt(dP_drho + dP_dT ^^2 * T/rho^^2/cV);
	}


} // end class Ideal_gas

unittest {
    import std.stdio;
    auto gm = new IdealGas();
    assert(gm.species_name(0) == "ideal air", "species name list");
    auto gd = new GasState(gm, 100.0e3, 300.0);
    assert(approxEqual(gm.R(gd), 287.086), failedUnitTest());
    assert(gm.n_modes == 1, failedUnitTest());
    assert(gm.n_species == 1, failedUnitTest());
    assert(approxEqual(gd.p, 1.0e5), failedUnitTest());
    assert(approxEqual(gd.T[0], 300.0), failedUnitTest());
    assert(approxEqual(gd.massf[0], 1.0), failedUnitTest());

    gm.update_thermo_from_pT(gd);
    gm.update_sound_speed(gd);
    assert(approxEqual(gd.rho, 1.16109), failedUnitTest());
    assert(approxEqual(gd.e[0], 215314.0), failedUnitTest());
    assert(approxEqual(gd.a, 347.241), failedUnitTest());
    gm.update_trans_coeffs(gd);
    assert(approxEqual(gd.mu, 1.84691e-05), failedUnitTest());
    assert(approxEqual(gd.k[0], 0.0262449), failedUnitTest());

    lua_State* L = init_lua_State("sample-data/ideal-air-gas-model.lua");
    gm = new IdealGas(L);
    lua_close(L);
    assert(approxEqual(gm.R(gd), 287.086), failedUnitTest());
    assert(gm.n_modes == 1, failedUnitTest());
    assert(gm.n_species == 1, failedUnitTest());
    assert(approxEqual(gd.p, 1.0e5), failedUnitTest());
    assert(approxEqual(gd.T[0], 300.0), failedUnitTest());
    assert(approxEqual(gd.massf[0], 1.0), failedUnitTest());

    gm.update_thermo_from_pT(gd);
    gm.update_sound_speed(gd);
    assert(approxEqual(gd.rho, 1.16109), failedUnitTest());
    assert(approxEqual(gd.e[0], 215314.0), failedUnitTest());
    assert(approxEqual(gd.a, 347.241), failedUnitTest());
    gm.update_trans_coeffs(gd);
    assert(approxEqual(gd.mu, 1.84691e-05), failedUnitTest());
    assert(approxEqual(gd.k[0], 0.0262449), failedUnitTest());
}
