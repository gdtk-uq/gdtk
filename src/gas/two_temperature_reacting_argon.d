/**
 * powers_aslam_gas.d
 *
 * Two-temperature reacting argon based off of:
 * "Quasi-One-Dimensional, Nonequilibrium Gas Dynamics of Partially Ionised Two-Temperature Argon"
 * Martin I. Hoffert and Hwachii Lien
 * 
 * The Physics of Fluids 10, 1769 (1967); doi 10.1063/1.1762356
 *
 *
 * Authors: Daniel Smith and Rory Kelly
 * Version: 19-July-2017: initial cut.
 */

module gas.two_temperature_reacting_argon;

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

class TwoTemperatureReactingArgon: GasModel {
public:

    this(lua_State *L) {
	// Some parameters are fixed and some come from the gas model file.
	_n_species = 3;
	_n_modes = 1;
	_species_names.length = 3;
	_species_names[0] = "Ar";
	_species_names[1] = "Ar_plus";
	_species_names[2] = "e_minus";
	lua_getglobal(L, "TwoTemperatureReactingArgon");
	// [TODO] test that we actually have the table as item -1
	// Now, pull out the remaining numeric value parameters.
	_Rgas = getDouble(L, -1, "R");
	_theta_ion = getDouble(L, -1, "theta_ion");
	_theta_A1star = getDouble(L, -1, "theta_A1star");
	_mol_masses.length = 3;
	_mol_masses[0] = 39.948e-3; // Units are kg/mol
	_mol_masses[2] = 5.485799e-7; // Units are kg/mol
	_mol_masses[1] = _mol_masses[0] - _mol_masses[2]; // Units are kg/mol
	_Cp = 520;
	_Cv = 312;
	lua_pop(L, 1); // dispose of the table
	// Compute derived parameters
	create_species_reverse_lookup();
    } // end constructor

    override string toString() const
    {
	char[] repr;
	repr ~= "TwoTemperatureReactingArgon =(";
	repr ~= "species=[\"Ar\", \"Ar_plus\", \"e_minus\"]";
	repr ~= ", Mmass=[" ~ to!string(_mol_masses[0]);
	repr ~= "," ~ to!string(_mol_masses[1]);
	repr ~= "," ~ to!string(_mol_masses[2]) ~ "]";
	repr ~= ")";
	return to!string(repr);
    }

    override void update_thermo_from_pT(GasState Q) const 
    {
    	double alpha = (Q.massf[2]/_mol_masses[2]) / ((Q.massf[2]/_mol_masses[2])+(Q.massf[0]/_mol_masses[0]));
	Q.rho = Q.p/(_Rgas*(Q.Ttr + alpha*Q.T_modes[0]));				//update rho      			Q.p/(Q.Ttr*_Rgas); //ideal gas law
	Q.u = 3/2*_Rgas*Q.Ttr;		 					//update internal energy of heavy particles 	//_Cv*Q.Ttr- Q.massf[1]*_q;
	Q.e_modes[0] = 3/2*_Rgas*alpha*Q.T_modes[0]+alpha*_Rgas*_theta_ion;	//update energy of electronic mode...
    }
    override void update_thermo_from_rhou(GasState Q) const
    {
    	double alpha = (Q.massf[2]/_mol_masses[2]) / ((Q.massf[2]/_mol_masses[2])+(Q.massf[0]/_mol_masses[0]));
	Q.Ttr = 2/3*Q.u/_Rgas;								//Q.u/_Cv;//(Q.u + Q.massf[1]*_q)/_Cv;
	Q.T_modes[0] = 	2/3*(Q.e_modes[0]-alpha*_Rgas*_theta_ion)/(_Rgas*alpha);		//
	Q.p = Q.rho*_Rgas*(Q.Ttr+alpha*Q.T_modes[0]);					// Q.rho*_Rgas*Q.Ttr;
    }
    override void update_thermo_from_rhoT(GasState Q) const
    {
    	double alpha = (Q.massf[2]/_mol_masses[2]) / ((Q.massf[2]/_mol_masses[2])+(Q.massf[0]/_mol_masses[0]));
	Q.p = Q.rho*_Rgas*(Q.Ttr+alpha*Q.T_modes[0]);	//Q.rho*_Rgas*Q.Ttr;
	Q.u = 3/2*_Rgas*Q.Ttr;			//_Cv*Q.Ttr; //- Q.massf[1]*_q;
    	Q.e_modes[0] = 3/2*_Rgas*alpha*Q.T_modes[0] + alpha*_Rgas*_theta_ion;		//
    }
    override void update_thermo_from_rhop(GasState Q) const
    {
    	double alpha = (Q.massf[2]/_mol_masses[2]) / ((Q.massf[2]/_mol_masses[2])+(Q.massf[0]/_mol_masses[0]));
	Q.Ttr = Q.p/Q.rho/_Rgas-alpha*Q.T_modes[0];		// Q.p/(Q.rho*_Rgas);// Assume Q.T_modes[0] is set independently, and correct.
	// Assume Q.T_modes[0] is set independently, and correct.
	Q.u = 3/2*_Rgas*Q.Ttr;			// _Cv*Q.Ttr; //- Q.massf[1]*_q;
	Q.e_modes[0] = 3/2*_Rgas*alpha*Q.T_modes[0] + alpha*_Rgas*_theta_ion;		//
    }
    
    override void update_thermo_from_ps(GasState Q, double s) const
    {
	throw new GasModelException("update_thermo_from_ps not implemented in TwoTemperatureReactingArgon.");
    }
    override void update_thermo_from_hs(GasState Q, double h, double s) const
    {
	throw new GasModelException("update_thermo_from_hs not implemented in TwoTemperatureReactingArgon.");
    }
    override void update_sound_speed(GasState Q) const
    {
    	double gamma = dhdT_const_p(Q)/dudT_const_v(Q);
	Q.a = sqrt(gamma*_Rgas*Q.Ttr);
    }
    override void update_trans_coeffs(GasState Q)
    {
	// TODO (DANIEL) - Should add these in such that viscous effects can be modelled
	Q.mu = 0.0;
	Q.k = 0.0;
	Q.k_modes[0] = 0.0;
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
	return _Rgas*Q.Ttr; //TODO (Daniel) Check this
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
	throw new GasModelException("entropy not implemented in TwoTemperatureReactingArgon.");
    }

private:
    // Thermodynamic constants
    double _Rgas; // J/kg/K
    double _theta_ion;
    double _theta_A1star;
    double _Cp;
    double _Cv;

} // end class PowersAslamGas


// Now, for the reaction update...
//
// It is included here because it is a small amount of code and
// it is specific to this particular gas model.

final class UpdateArgonFrac : ThermochemicalReactor {
    
    this(string fname, GasModel gmodel)
    {
	super(gmodel); // hang on to a reference to the gas model
	// We need to pick a number of pieces out of the gas-model file, again.
	// Although they exist in the GasModel object, they are private.
	auto L = init_lua_State();
	doLuaFile(L, fname);
	lua_getglobal(L, "TwoTemperatureReactingArgon");
	_Rgas = getDouble(L, -1, "R");
	_theta_ion = getDouble(L, -1, "theta_ion");
	_theta_A1star = getDouble(L, -1, "theta_A1star");
	_mol_masses.length = 3;
	_mol_masses[0] = 39.948e-3; // Units are kg/mol
	_mol_masses[2] = 5.485799e-7; // Units are kg/mol
	_mol_masses[1] = _mol_masses[0] - _mol_masses[2]; // Units are kg/mol	
	lua_pop(L, 1); // dispose of the table
	lua_close(L);
    }
    
    override void opCall(GasState Q, double tInterval, ref double dtSuggest,
			 ref double[] params)
    {
 //   	double m_Ar = 6.6335209e-26; //mass of argon (kg)
 //   	double m_e = 9.10938e-31; //mass of electron (kg)

 //   	double alpha;

	//double n_e;
	//double n_Ar;
 //       double kfA;
 //       double kfe;
 //       double krA;
 //       double kre;

 //      	double ne_dot_A;
 //   	double ne_dot_e; 
 //   	double ne_dot;

 //       double Q_ea; 
	//double Q_ei;

	//double v_ea;
 //       double v_ei;

 //       double u_trans_ionisation_e;
	//double u_trans_ionisation_Ar;
 //       double u_trans_collisions;

 //       double NumberSteps = 10000000;
 //       double dt = tInterval/NumberSteps;

 //   	for (int number = 1; number <= NumberSteps; ++number) { // This is a simple euler integration...
 //   		writeln(number);
 //         	 alpha = (Q.massf[2]/_mol_masses[2]) / ((Q.massf[2]/_mol_masses[2])+(Q.massf[0]/_mol_masses[0]));

	//    	//Determine the current number densities
	//        n_e = Q.rho/_mol_masses[2]*Q.massf[2]*Av; // number density of electrons
	//	n_Ar = Q.rho/_mol_masses[0]*Q.massf[0]*Av; // number density of Ar

	//    	//Determine the current rate coefficients
	//        kfA = 1.68e-26*pow(Q.Ttr,1.5)*(_theta_A1star/Q.Ttr+2)*exp(-_theta_A1star/Q.Ttr);
	//        kfe = 3.75e-22*pow(Q.T_modes[0],1.5)*(_theta_A1star/Q.T_modes[0]+2)*exp(-_theta_A1star/Q.T_modes[0]);
	//        krA = 5.8e-49*(_theta_A1star/Q.Ttr+2)*exp((_theta_ion-_theta_A1star)/Q.Ttr);
	//        kre = 1.29e-44*(_theta_A1star/Q.T_modes[0]+2)*exp((_theta_ion-_theta_A1star)/Q.T_modes[0]);

	//    	//determine the current rate of ionisation
	//    	ne_dot_A = kfA*pow(n_Ar,2) - krA*n_Ar*pow(n_e,2);	//production rate of electrons due to argon collisions
	//    	ne_dot_e = kfe*n_Ar*n_e - kre*pow(n_e,3);	//production rate of electrons due to electron collisions
	//    	ne_dot = ne_dot_A + ne_dot_e;

	//    	//determine the new number densities
	//    	n_e = n_e + ne_dot*dt;
	//    	n_Ar = n_Ar - ne_dot*dt;

	//    	//convert back to mass fractions
	//        Q.massf[0] = n_Ar/Av/Q.rho*_mol_masses[0];
	//        Q.massf[1] = n_e/Av/Q.rho*_mol_masses[1]; //number density of Argon+ is the same as electron number density
	//        Q.massf[2] = n_e/Av/Q.rho*_mol_masses[2];

	//    	//find the collision cross sectional area
	    	
	//    	if (Q.T_modes[0] < 1.0e4) {
	//    		Q_ea = (0.39-0.551e-4*Q.T_modes[0]+0.595e-8*pow(Q.T_modes[0],2))*pow(10,-20);
	//	} else if (Q.T_modes[0]<1.0e5) {
	//		Q_ea = (-0.35+0.775e-4*Q.T_modes[0])*pow(10,-20);
	//	} else {
	//		throw new GasModelException("Electron Temperature Too high for collision cross section model");
	//	}
	    	 
	//    	Q_ei = 1.95e-10*pow(Q.T_modes[0],-2)*log(1.53e8*pow(Q.T_modes[0],3)/n_e);
	//        //find the collision frequencies
	//        v_ea = (1-alpha)*Q.rho/m_Ar*pow(8*Kb*Q.T_modes[0]/PI/m_e,0.5)*Q_ea;//electron-Ar collisions
	//        v_ei = alpha*Q.rho/m_Ar*pow(8*Kb*Q.T_modes[0]/PI/m_e,0.5)*Q_ei;//electron-Ar+ collisions

	//        //update the energy of each state
	//        u_trans_ionisation_e = ne_dot_e*Kb*_theta_ion*dt;// energy transferred to electron mode through ionisation
	//        u_trans_ionisation_Ar= ne_dot_A*Kb*_theta_ion*dt;
	//        u_trans_collisions = 3*n_e*m_e/m_Ar*(v_ea+v_ei)*Kb*(Q.Ttr-Q.T_modes[0])*dt;// energy transferred to electron mode through collisions
	//        Q.u = Q.u - u_trans_collisions;
	//        Q.e_modes[0] = Q.e_modes[0] + u_trans_collisions - u_trans_ionisation_e;
	//	_gmodel.update_thermo_from_rhou(Q);

 //   	}

	//// Since the internal energy and density in the (isolated) reactor is fixed,
	//// we need to evaluate the new temperature, pressure, etc.
	_gmodel.update_sound_speed(Q);  
      
    }

private:
    // Reaction rate constant
    double _alpha; // 1/s
    double _Rgas;
    double _theta_ion;
    double _theta_A1star;
    double[] _mol_masses;
} // end class UpdateAB


//// Unit test of the basic gas model...

version(two_temperature_reacting_argon_test) {
    import std.stdio;
    import util.msg_service;
    int main() {
	writeln("THIS IS PERFORMING THE UNIT TEST");
	lua_State* L = init_lua_State();
	doLuaFile(L, "sample-data/two_temperature_reacting_argon-model.lua");
	auto gm = new TwoTemperatureReactingArgon(L);
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
