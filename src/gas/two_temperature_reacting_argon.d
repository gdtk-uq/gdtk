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
	_ion_tol = getDouble(L, -1, "ion_tol");
	_T_modes_ref = getDouble(L, -1, "T_modes_ref");

	_mol_masses.length = 3;
	_mol_masses[0] = 39.948e-3; // Units are kg/mol
	_mol_masses[2] = 5.485799e-7; // Units are kg/mol
	_mol_masses[1] = _mol_masses[0] - _mol_masses[2]; // Units are kg/mol
	_Cp = 520.0;
	_Cv = 312.0;
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
	Q.rho = Q.p/(_Rgas*(Q.Ttr + alpha*Q.T_modes[0]));
	Q.u = 3.0/2.0*_Rgas*Q.Ttr;
	if (alpha<=_ion_tol) {
		Q.T_modes[0] = _T_modes_ref;
		Q.u_modes[0] = 3.0/2.0*_Rgas*alpha*Q.T_modes[0] + alpha*_Rgas*_theta_ion;
	} else {
		Q.u_modes[0] = 3.0/2.0*_Rgas*alpha*Q.T_modes[0] + alpha*_Rgas*_theta_ion;
	}
    }
    override void update_thermo_from_rhou(GasState Q) const
    {
    	double alpha = (Q.massf[2]/_mol_masses[2]) / ((Q.massf[2]/_mol_masses[2])+(Q.massf[0]/_mol_masses[0]));
	Q.Ttr = 2.0/3.0*Q.u/_Rgas;
	if (alpha <= _ion_tol) {
		Q.T_modes[0] = _T_modes_ref;
		Q.u_modes[0] = 3.0/2.0*_Rgas*alpha*Q.T_modes[0] + alpha*_Rgas*_theta_ion;
	} else {
		Q.T_modes[0] = (Q.u_modes[0]/alpha-_Rgas*_theta_ion)*2.0/3.0/_Rgas;
	}
	Q.p = Q.rho*_Rgas*(Q.Ttr+alpha*Q.T_modes[0]);					// Q.rho*_Rgas*Q.Ttr;
    }

    override void update_thermo_from_rhoT(GasState Q) const
    {
    	double alpha = (Q.massf[2]/_mol_masses[2]) / ((Q.massf[2]/_mol_masses[2])+(Q.massf[0]/_mol_masses[0]));
	Q.p = Q.rho*_Rgas*(Q.Ttr+alpha*Q.T_modes[0]);	//Q.rho*_Rgas*Q.Ttr;
	Q.u = 3.0/2.0*_Rgas*Q.Ttr;
	if (alpha <= _ion_tol) {
		Q.T_modes[0] = _T_modes_ref;
		Q.u_modes[0] = 3.0/2.0*_Rgas*alpha*Q.T_modes[0] + alpha*_Rgas*_theta_ion;
	} else {
		Q.u_modes[0] = 3.0/2.0*_Rgas*alpha*Q.T_modes[0] + alpha*_Rgas*_theta_ion;
	}
    }

    override void update_thermo_from_rhop(GasState Q) const
    {
    	double alpha = (Q.massf[2]/_mol_masses[2]) / ((Q.massf[2]/_mol_masses[2])+(Q.massf[0]/_mol_masses[0]));
	Q.Ttr = Q.p/Q.rho/_Rgas - alpha*Q.T_modes[0];
	// Assume Q.T_modes[0] is set independently, and correct.
	Q.u = 3.0/2.0*_Rgas*Q.Ttr;
	if (alpha <= _ion_tol) {
		Q.T_modes[0] = _T_modes_ref;
		Q.u_modes[0] = 3.0/2.0*_Rgas*alpha*Q.T_modes[0] + alpha*_Rgas*_theta_ion;
	} else {
		Q.u_modes[0] = 3.0/2.0*_Rgas*alpha*Q.T_modes[0] + alpha*_Rgas*_theta_ion;
	}
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
    	double _gamma = dhdT_const_p(Q)/dudT_const_v(Q);
	Q.a = sqrt(_gamma*_Rgas*Q.Ttr);
	//[TODO] update the _Cv and _Cp properties to be dependant on alpha...
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
    	double alpha = (Q.massf[2]/_mol_masses[2]) / ((Q.massf[2]/_mol_masses[2])+(Q.massf[0]/_mol_masses[0]));
	return 3.0/2.0*_Rgas*(1+alpha) + alpha-_Rgas*alpha*(1-alpha)/(2-alpha)*pow((3.0/2.0*Q.Ttr+alpha*_theta_ion)/Q.Ttr,2);
    }
    override double dhdT_const_p(in GasState Q) const
    {
    	double alpha = (Q.massf[2]/_mol_masses[2]) / ((Q.massf[2]/_mol_masses[2])+(Q.massf[0]/_mol_masses[0]));
    	return 5.0/2.0*_Rgas*(1+alpha) + _Rgas/2*alpha*(1-pow(alpha,2))*pow((5.0/2.0*Q.Ttr+alpha*_theta_ion)/Q.Ttr,2);
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
    double _ion_tol;
    double _T_modes_ref;

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
	_T_modes_ref = getDouble(L, -1, "T_modes_ref");
	_mol_masses.length = 3;
	_mol_masses[0] = 39.948e-3; // Units are kg/mol
	_mol_masses[2] = 5.485799e-7; // Units are kg/mol
	_mol_masses[1] = _mol_masses[0] - _mol_masses[2]; // Units are kg/mol	
	_ion_tol = getDouble(L, -1, "ion_tol");
	_chem_dt = getDouble(L, -1, "chem_dt");
	lua_pop(L, 1); // dispose of the table
	lua_close(L);
    }
    
    override void opCall(GasState Q, double tInterval, ref double dtSuggest,
			 ref double[] params)
    {
    	if (Q.Ttr > 3000.0) {
    	//writeln("Running a single step of the chemistry update");

    	double m_Ar = 6.6335209e-26; //mass of argon (kg)
    	double m_e = 9.10938e-31; //mass of electron (kg)
        
        double Kb = Boltzmann_constant;
        double Av = Avogadro_number;
    	
    	double alpha;

	double n_e;
	double n_Ar;
        double kfA;
        double kfe;
        double krA;
        double kre;

       	double ne_dot_A;
    	double ne_dot_e; 
    	double ne_dot;

        double Q_ea; 
	double Q_ei;

	double v_ea;
        double v_ei;

        double u_trans_ionisation;
        double u_trans_ionisation_heavy;
        double u_trans_ionisation_electron;
        double u_trans_collisions;

        double _chem_dt = 1.0e-11;
        int NumberSteps = to!int(tInterval/_chem_dt + 1);
        if (NumberSteps == 0) {NumberSteps = 1;}
        _chem_dt = tInterval/NumberSteps;
        //writeln("_chem_dt = ", _chem_dt);
  	//alpha = (Q.massf[2]/_mol_masses[2]) / ((Q.massf[2]/_mol_masses[2])+(Q.massf[0]/_mol_masses[0]));

    	//Determine the current number densities
        n_e = Q.rho/_mol_masses[2]*Q.massf[2]*Av; // number density of electrons
	n_Ar = Q.rho/_mol_masses[0]*Q.massf[0]*Av; // number density of Ar
    	alpha = n_e/(n_e + n_Ar);

    	double orig_n = n_e + n_Ar;
    	double n_sum;

	double internal_energy_initial = 3.0/2.0*_Rgas*(Q.Ttr+alpha*Q.T_modes[0])+alpha*_Rgas*_theta_ion;

        //writeln("mass fraction electrons = ", Q.massf[2]);

	//writeln("n_e = ", n_e);
	//writeln("n_Ar = ", n_Ar);

	//writeln("theta_A1star = ", _theta_A1star);
	//writeln("theta_ion = ", _theta_ion);
	//writeln("theta_A1star - theta_ion", _theta_A1star-_theta_ion);
	//writeln("number of steps = ", NumberSteps);
    	for (int number = 1; number <= NumberSteps; ++number) { // This is a simple euler integration...
    		//writeln("number = ", number);
    		//CHEMISTRY PART----------------------------------------------------------------------------------
	    	//Determine the current rate coefficients
	        kfA = 1.68e-26*Q.Ttr*sqrt(Q.Ttr)*(_theta_A1star/Q.Ttr+2)*exp(-_theta_A1star/Q.Ttr);
	        kfe = 3.75e-22*Q.T_modes[0]*sqrt(Q.T_modes[0])*(_theta_A1star/Q.T_modes[0]+2)*exp(-_theta_A1star/Q.T_modes[0]);
	        krA = 5.8e-49*(_theta_A1star/Q.Ttr+2)*exp((_theta_ion-_theta_A1star)/Q.Ttr);
	        kre = 1.29e-44*(_theta_A1star/Q.T_modes[0]+2)*exp((_theta_ion-_theta_A1star)/Q.T_modes[0]);

	        //writeln("kfA = ", kfA);
	        //writeln("kfe = ", kfe);
	        //writeln("krA = ", krA);
	        //writeln("kre = ", kre);

	    	//determine the current rate of ionisation
	    	ne_dot_A = kfA*pow(n_Ar,2) - krA*n_Ar*pow(n_e,2);	//production rate of electrons due to argon collisions
	    	ne_dot_e = kfe*n_Ar*n_e - kre*pow(n_e,3);	//production rate of electrons due to electron collisions
	    	ne_dot = ne_dot_A + ne_dot_e;

	    	//writeln("ne_dot_A = ", ne_dot_A);
	    	//writeln("ne_dot_e = ", ne_dot_e);
	    	//writeln("ne_dot = ", ne_dot);

	    	//determine the new number densities
	    	n_e = n_e + ne_dot*_chem_dt;
	    	n_Ar = n_Ar - ne_dot*_chem_dt;

	    	alpha = n_e/(n_e + n_Ar);

	    	//writeln("n_e = ", n_e);
	    	//writeln("n_Ar = ", n_Ar);
	    	//writeln("alpha = ", alpha);

		//update energy based on ionisation process
		//Come to this conceptually. 

		if (alpha <= _ion_tol) {
			Q.T_modes[0] = _T_modes_ref;
			Q.u_modes[0] = 3.0/2.0*_Rgas*alpha*Q.T_modes[0] + alpha*_Rgas*_theta_ion;
		} else {
			Q.u_modes[0] = 3.0/2.0*_Rgas*alpha*Q.T_modes[0] + alpha*_Rgas*_theta_ion - ne_dot_e*Kb*_theta_ion*_chem_dt/Q.rho;
			Q.T_modes[0] = (Q.u_modes[0]/alpha-_Rgas*_theta_ion)*2.0/3.0/_Rgas;
		}

		Q.u = internal_energy_initial - Q.u_modes[0];
		Q.Ttr = 2.0/3.0*Q.u/_Rgas;

	    	if (n_e <= 0.0) { // if the number densities of electrons go below zero, then force this to not occur and update thermo.
	    		Q.u = internal_energy_initial;
	    		Q.u_modes[0] = 0.0;
	    		//mass
	    		n_Ar = orig_n;
	    		n_e = 0.0;
	    		alpha = n_e/(n_e + n_Ar);
	    		//energy
	    		//temperature
			Q.Ttr =  2.0/3.0*Q.u/_Rgas;
	    		Q.T_modes[0] = _T_modes_ref;
	    	} else {
	    	n_sum = orig_n/(n_e + n_Ar);

	    	n_e = n_e*n_sum;
	    	n_Ar = n_Ar*n_sum;
	    	//THERMAL RELAXATION PART----------------------------------------------------------------------------------
	    	//find the collision cross sectional area

		if (alpha > _ion_tol) { // no thermal relaxation if below ion tolerance
	    	if (Q.T_modes[0] < 1.0e4) {
	    		Q_ea = (0.39 - 0.551e-4*Q.T_modes[0] + 0.595e-8*pow(Q.T_modes[0],2))*1.0e-20;
		} else if (Q.T_modes[0]<1.0e5) {
			Q_ea = (-0.35 + 0.775e-4*Q.T_modes[0])*1.0e-20;
		} else {
			Q_ea = (-0.35 + 0.775e-4*50000)*1.0e-20;
			//writeln("Q.u = ", Q.u);
			//writeln("Q.u_modes[0] = ", Q.u_modes[0]);
			//writeln("Q.Ttr = ", Q.Ttr);
			//writeln("Q.T_modes[0] = ", Q.T_modes[0]);
			////writeln("u_trans_ionisation = ", u_trans_ionisation);
			//writeln("ne_dot ", ne_dot);
			//writeln("Q.rho = ", Q.rho);
			//writeln("kfA = ", kfA);
			//writeln("kfe = ", kfe);
			//writeln("krA = ", krA);
			//writeln("kre = ", kre);
	  //      	writeln("number = ", number);
			//throw new GasModelException("Electron Temperature Too high for collision cross section model");
		}

	    	Q_ei = 1.95e-10*pow(Q.T_modes[0],-2)*log(1.53e8*pow(Q.T_modes[0],3)/(n_e/1.0e6));
	    	if (Q_ei < 0.0) {Q_ei = 0.0;}

	        //find the collision frequencies
	        v_ea = (1-alpha)*Q.rho/m_Ar*sqrt(8*Kb*Q.T_modes[0]/PI/m_e)*Q_ea; //electron-Ar collisions
	        v_ei = alpha*Q.rho/m_Ar*sqrt(8*Kb*Q.T_modes[0]/PI/m_e)*Q_ei; //electron-Ar+ collisions

	        //update the energy of each state
	        u_trans_collisions    = 3*n_e*m_e/m_Ar*(v_ea+v_ei)*Kb*(Q.Ttr-Q.T_modes[0])*_chem_dt/Q.rho;// energy transferred to electron mode through collisions
	        //writeln("u_trans_collisions = ", u_trans_collisions);
	        Q.u -= u_trans_collisions;
	        Q.u_modes[0] += u_trans_collisions;

	        //update thermo properties based on energy transfer
		Q.Ttr = 2.0/3.0*Q.u/_Rgas;
		Q.T_modes[0] = (Q.u_modes[0]/alpha-_Rgas*_theta_ion)*2.0/3.0/_Rgas;

		} // end if alpha > ion tol
		} // end if statement regarding number density of electrons

    	}

    	//convert back to mass fractions //Density has not changed since finite volume cell with  no flux
        Q.massf[0] = n_Ar/Av/Q.rho*_mol_masses[0];
        Q.massf[1] = n_e/Av/Q.rho*_mol_masses[1]; //number density of Argon+ is the same as electron number density
        Q.massf[2] = n_e/Av/Q.rho*_mol_masses[2];

	if ((Q.massf[0] + Q.massf[1] + Q.massf[2]) < 0.5) {
		writeln("mass fraction sum = ", Q.massf[0] + Q.massf[1] + Q.massf[2]);
		writeln("Q.u = ", Q.u);
		writeln("Q.u_modes[0] = ", Q.u_modes[0]);
		writeln("Q.Ttr = ", Q.Ttr);
		writeln("Q.T_modes[0] = ", Q.T_modes[0]);
		//writeln("u_trans_ionisation = ", u_trans_ionisation);
		writeln("ne_dot =", ne_dot);
		writeln("n_e = ", n_e);
		writeln("n_Ar = ", n_Ar);
		writeln("Q.rho = ", Q.rho);
		writeln("kfA = ", kfA);
		writeln("kfe = ", kfe);
		writeln("krA = ", krA);
		writeln("kre = ", kre);
		writeln("_chem_dt = ", _chem_dt);
		writeln("NumberSteps = ", NumberSteps);
		writeln("tInterval = ", tInterval);
		throw new GasModelException("mass fraction far from 1...");
		} 

	// Since the internal energy and density in the (isolated) reactor is fixed,
	// we need to evaluate the new temperature, pressure, etc.

	if (Q.T_modes[0] < 0) {
		throw new GasModelException("T_modes is negative!!!!!");
	}
	_gmodel.update_thermo_from_rhou(Q);
	_gmodel.update_sound_speed(Q);  
      	} // end if
    }

private:
    // Reaction rate constant
    double _Rgas;
    double _theta_ion;
    double _theta_A1star;
    double[] _mol_masses;
    double _ion_tol;
    double _T_modes_ref;
    double _chem_dt;
} // end class UpdateAB


//// Unit test of the basic gas model...

version(two_temperature_reacting_argon_test) {
    import std.stdio;
    import util.msg_service;
    int main() {
    	writeln("Beginning the unit test...");
    	writeln("Testing the gas state functions...");
	lua_State* L = init_lua_State();
	doLuaFile(L, "sample-data/two-temperature-reacting-argon-model.lua");
	auto gm = new TwoTemperatureReactingArgon(L);
	lua_close(L);
	auto gd = new GasState(3, 1);
	gd.p = 1.0e5;
	gd.Ttr = 300.0;
	gd.T_modes[0] = 300;
	gd.massf[0] = 1.0; gd.massf[1] = 0.0; gd.massf[2] = 0.0;

	assert(approxEqual(gm.R(gd), 208.0, 1.0e-4), failedUnitTest());
	assert(gm.n_modes == 1, failedUnitTest());
	assert(gm.n_species == 3, failedUnitTest());
	assert(approxEqual(gd.p, 1.0e5, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.Ttr, 300.0, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.massf[0], 1.0, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.massf[1], 0.0, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.massf[2], 0.0, 1.0e-6), failedUnitTest());

	gm.update_thermo_from_pT(gd);
	gm.update_sound_speed(gd);
	double my_rho = 1.0e5 / (208.0 * 300.0);
	assert(approxEqual(gd.rho, my_rho, 1.0e-4), failedUnitTest());

	double my_Cv = gm.dudT_const_v(gd);
	double my_u = my_Cv*300.0; 
	assert(approxEqual(gd.u, my_u, 1.0e-3), failedUnitTest());

	double my_Cp = gm.dhdT_const_p(gd);
	double my_a = sqrt(my_Cp/my_Cv*208.0*300.0);
	assert(approxEqual(gd.a, my_a, 1.0e-3), failedUnitTest());

	gm.update_trans_coeffs(gd);
	assert(approxEqual(gd.mu, 0.0, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.k, 0.0, 1.0e-6), failedUnitTest());
	assert(approxEqual(gd.k, 0.0, 1.0e-6), failedUnitTest());

		//writeln("Testing the chemistry...");

		//gd.Ttr = 4000; // Set temperature to 3000K - this should be enmough to initiate some ionisation	
		//gm.update_thermo_from_pT(gd); // update gas state
		//gm.update_sound_speed(gd); // upate sound speed (not necessary)

		//writeln("------------------------------------------");
		//writeln("Before:");
		//writeln("mass fraction argon = ", gd.massf[0]);
		//writeln("mass fraction Argon+ = ", gd.massf[1]);
		//writeln("mass fraction electron = ", gd.massf[2]);
		//writeln("Pressure = ", gd.p);
		//writeln("Temperature = ", gd.Ttr);
		//writeln("T modes = ", gd.T_modes[0]);
		//writeln("density = ", gd.rho);
		//writeln("internal energy = ", gd.u);
		//writeln("U modes = ", gd.u_modes[0]);

		//auto reactor = new UpdateArgonFrac("sample-data/two-temperature-reacting-argon-model.lua", gm);
		//double[] params;
		//double dtSuggest;
		//writeln("----");
		//reactor.opCall(gd,1.0,dtSuggest,params);

		//writeln("----");
		//writeln("After:");
		//writeln("mass fraction argon = ", gd.massf[0]);
		//writeln("mass fraction Argon+ = ", gd.massf[1]);
		//writeln("mass fraction electron = ", gd.massf[2]);
		//writeln("Pressure = ", gd.p);
		//writeln("Temperature = ", gd.Ttr);
		//writeln("T modes = ", gd.T_modes[0]);
		//writeln("density = ", gd.rho);
		//writeln("internal energy = ", gd.u);
		//writeln("U modes = ", gd.u_modes[0]);
		//writeln("Unit test completed");

	writeln("======================================================================");
	writeln("Doing the chemistry relaxation problem");

	auto reactor = new UpdateArgonFrac("sample-data/two-temperature-reacting-argon-model.lua", gm);
	double[] params;
	double dtSuggest;

	//need molar masses to determine alpha
	double[3] _mol_masses;
	_mol_masses[0] = 39.948e-3; // Units are kg/mol
	_mol_masses[2] = 5.485799e-7; // Units are kg/mol
	_mol_masses[1] = _mol_masses[0] - _mol_masses[2]; // Units are kg/mol
	double _theta_ion = 183100.0;
	double _Rgas = 208;
	double alpha;
	double new_vel; 

	//pre-shock conditions
	double rho1 = 0.021366;
	double T1 = 300;
	double u1 = 6.1e3;
	double M1 = 18.9;

	//inital gas properties
	auto GS = new GasState(3, 1);
	GS.p = GS.p = 594945.77;
	GS.Ttr = GS.Ttr = 33750.78;
	GS.T_modes[0] = 10000;
	GS.massf[0] = 1.0; GS.massf[1] = 0.0; GS.massf[2] = 0.0;
	double vel = 1537; // m/s
	double x = 0.0;

	double e_new;

	//update the gas model based on the properties specified above
	gm.update_thermo_from_pT(GS); // update gas state
	gm.update_sound_speed(GS); // upate sound speed (not necessary)

	// some time stepping information
	double maxtime = 4.0e-6; // max time for the simulation
	double dt = 1.312e-11; // time step for chemistry update
	int maxsteps = to!int(maxtime/dt + 1); // number of steps in which to iterate through
	int printfreq = maxsteps/10;
	int writefreq = maxsteps/1000;
	//writeln(maxsteps);

	//initialise the storage arrays
	double[] t_list;
	double[] x_list;
	double[] T_list;
	double[] T_modes_list;
	double[] P_list;
	double[] u_list;
	double[] alpha_list;

	//initial values for storage arrays
	t_list ~= 0.0;
	x_list ~= x;
	T_list ~= GS.Ttr;
	T_modes_list ~= GS.T_modes[0];
	P_list ~= GS.p;
	u_list ~= vel;
	alpha_list ~= 0.0;

	//main for loop for chemistry
	for (int i = 1; i <= maxsteps; ++i) {
		reactor.opCall(GS,dt,dtSuggest,params); // perform the chemistry update
		//new alpha
		alpha = (GS.massf[2]/_mol_masses[2]) / ((GS.massf[2]/_mol_masses[2])+(GS.massf[0]/_mol_masses[0]));
		//update x position
		new_vel = u1/8.*(5 + 3.0/pow(M1,2) - sqrt(9*pow(1-1.0/pow(M1,2),2) + 
			96*alpha/5./pow(M1,2)*_theta_ion/T1));
		//writeln("new vel = ", new_vel);
		x += (vel + new_vel)/2*dt;
		vel = new_vel;
		//update rho, P and u.
		GS.rho = rho1*u1/vel;
		GS.p = rho1*(_Rgas*T1 + pow(u1,2)) - GS.rho*pow(vel,2);
		//new internal energy
		e_new = 5*_Rgas*T1/2. + pow(u1,2)/2 - pow(vel,2)/2 - GS.p/GS.rho;
		//give all of the new energy to the heavy particles
		GS.u = GS.u + (e_new - (GS.u + GS.u_modes[0]));
		//update the temperature of the heavy paritcles based on this
		gm.update_thermo_from_rhou(GS); 



		if (writefreq == 0) {writefreq = 1;}
		if ((i % writefreq) == 0) { // only save 1000 points
			if ((i % printfreq) == 0) {
				writeln(to!double(i)/maxsteps*100, "% completed");
			}

			t_list ~= i*dt;
			x_list ~= x;
			T_list ~= GS.Ttr;
			T_modes_list ~= GS.T_modes[0];
			P_list ~= GS.p;
			u_list ~= vel;
			alpha_list ~= alpha;
		}
	}
	//writeln("T_list = ", T_list);
	//writeln("T_modes_list = ", T_modes_list);
	//writeln("alpha massf list = ", massf_electron_list);

	writeln("writing to Data file... please wait");
	double[][] collateddata = [t_list,x_list,T_list,T_modes_list,P_list,u_list,alpha_list];

	File file = File("/home/uqdsmi31/Documents/ArgonFiniteRateValidation/results","w");
	foreach (i; 0..t_list.length) {
		file.writeln(collateddata[0][i], " ", collateddata[1][i], " ", collateddata[2][i], " ", collateddata[3][i], " ", collateddata[4][i], " ", collateddata[5][i], " ", collateddata[6][i]);
	}

	//writeln("mass fraction sum = ", GS.massf[0] + GS.massf[1] + GS.massf[2]);


	//writeln(collateddata);
	//writeln(x_list);
	//writeln("pow(10.5,-2)",pow(10.5,-2));
	//writeln("1/pow(10.5,2)",1.0/pow(10.5,2));
	//the actual calculation step
	//reactor.opCall(gd,1.0,dtSuggest,params);


	return 0;
    }
}


//	    	if (isNaN(kfA)) {
//		    	writeln("Q.u = ", Q.u);
//		    	writeln("Q.Ttr = ", Q.Ttr);
//	    		throw new GasModelException("kfA not a number.");
//	    	}

//	    	writeln("--------------------------------");
//	    	writeln("ne_dot = ", ne_dot);
//	    	writeln("ne_dot_A = ", ne_dot_A);
//	    	writeln("ne_dot_e = ", ne_dot_e);
//	    	writeln("n_Ar = ", n_Ar);
//	    	writeln("n_e = ", n_e);
//    		writeln("kfe = ", kfe);
//    		writeln("kre = ", kre);
//    		writeln("T_modes[0] = ", Q.T_modes[0]);
//	    	writeln("--------------------------------");


//	    	if (isNaN(v_ea)) {
//	    		writeln("alpha = ", alpha);
//	    		writeln("Q.T_modes[0] = ", Q.T_modes[0]);
//	    		writeln("ne_dot = ", ne_dot);
//	    		//writeln("ne_dot = ", ne_dot);
//	    		//writeln("ne_dot = ", ne_dot);
//	    		writeln("kfe = ", kfe);
//	    		writeln("kre = ", kre);
//	    		writeln("kfA = ", kfA);
//	    		writeln("krA = ", krA);
//	    		writeln("ne_dot_e = ", ne_dot_e);
//	    		writeln("ne_dot_A = ", ne_dot_A);

//	    		writeln("n_e = ", n_e);
//	    		writeln("n_Ar = ", n_Ar);
//	    		throw new GasModelException("v_ea not a number.");
//	    	} 



//	    	if (isNaN(u_trans_collisions)) {
//	    		writeln("Q.T_modes[0] = ",Q.T_modes[0]);
//	    		writeln("Q.Ttr = ",Q.Ttr);
//	    		writeln("v_ea = ", v_ea);
//	    		writeln("v_ei = ", v_ei);
//	    		writeln("alpha = ", alpha);
//	    		writeln("Q_ea = ", Q_ea);
//	    		throw new GasModelException("u_trans_collisions not a number.");
//	    	}

//if (isNaN(Q_ea)) {
//	    		writeln("Q.T_modes[0] = ", Q.T_modes[0]);
//	    		writeln("n_e = ", n_e);
//	    		throw new GasModelException("Q_ea not a number.");
//	    	}



//	    		    	if (isNaN(ne_dot)) {
//	    		writeln("ne_dot_e = ", ne_dot_e);
//	    		writeln("ne_dot_A = ", ne_dot_A);
//	    		writeln("kfe = ", kfe);
//	    		writeln("kre = ", kre);
//	    		writeln("Q.T_modes[0] = ", Q.T_modes[0]);
//	    		throw new GasModelException("ne_dot not a number.");
//	    	}