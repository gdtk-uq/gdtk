/**
 * two_temperature_argon_kinetics.d
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

module kinetics.two_temperature_argon_kinetics;

import std.stdio : writeln;
import std.math : exp, log, pow, sqrt, PI;
import std.conv : to;

import gas;
import util.lua;
import util.lua_service;
import kinetics.thermochemical_reactor;

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


