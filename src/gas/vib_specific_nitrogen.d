/**
 * vib_specific_nitrogen.d
 * Authors: Rowan G., Katrina Sklavos and Peter J.
 *
 * This is a 10-level vibrationally-specific model for nitrogen
 * as descrbied in:
 *
 * Giordano, et al. (1997)
 * Vibrationally Relaxing Flow of N2 past an Infinite Cylinder
 * Journal of Thermophysics and Heat Transfer, vol 11, no 1, pp 27 - 35
 *
 */

module gas.vib_specific_nitrogen;

import std.algorithm.iteration;

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

// [TODO:KS]
// Presently the number of vibrational levels
// is hard-coded. We'd like to test at small number
// then increase this to 10.
// Eventually, we might remove the hard-coded number.
immutable int N_VIB_LEVELS = 3;

class VibSpecificNitrogen: GasModel {
public:
    this()
    {
	_n_species = N_VIB_LEVELS;
	_n_modes = N_VIB_LEVELS;
	_species_names.length = _n_species;
	foreach (isp; 0 .. _n_species) {
            _species_names[isp] = format("N2-vib-%d", isp);
        }
	create_species_reverse_lookup();
    }

    override string toString() const
    {
	char[] repr;
	repr ~= "VibSpecificNitrogen=()";
	return to!string(repr);
    }

    override void update_thermo_from_pT(GasState Q) const 
    {
	Q.rho = Q.p/(_R_N2*Q.Ttr);
	
	foreach (imode; 0 .. _n_modes) {
	     Q.u_modes[imode] = (Avogadro_number/_M_N2) * Q.massf[imode] * vib_energy(imode);
	}
	Q.u = 2.5 * _R_N2 * Q.Ttr;
	
    }
    override void update_thermo_from_rhou(GasState Q) const
    {
	Q.Ttr = 0.4 * Q.u * (1/_R_N2);
	Q.p = Q.rho * _R_N2 * Q.Ttr;
	Q.T_modes[0] = compute_Tvib(Q, 200.0, 1000.0, 1e-4);

    }
    override void update_thermo_from_rhoT(GasState Q) const
    {
	Q.p = Q.rho * _R_N2 * Q.Ttr;
	Q.u = 2.5 * _R_N2 * Q.Ttr;
	
	foreach (imode; 0 .. _n_modes) {
	     Q.u_modes[imode] = (Avogadro_number/_M_N2) * Q.massf[imode] * vib_energy(imode);
	}

    }
    override void update_thermo_from_rhop(GasState Q) const
    {
	Q.Ttr = 0.4 * Q.u * (1/_R_N2);
	Q.u = 2.5 * _R_N2 * Q.Ttr;
	foreach (imode; 0 .. _n_modes) {
	     Q.u_modes[imode] = (Avogadro_number/_M_N2) * Q.massf[imode] * vib_energy(imode);
	}
    }
    
    override void update_thermo_from_ps(GasState Q, double s) const
    {
	throw new Error("ERROR: VibSpecificNitrogen:update_thermo_from_ps NOT IMPLEMENTED.");

    }

    override void update_thermo_from_hs(GasState Q, double h, double s) const
    {
	throw new Error("ERROR: VibSpecificNitrogen:update_thermo_from_hs NOT IMPLEMENTED.");

    }
    override void update_sound_speed(GasState Q) const
    {
	Q.a = (_gamma * _R_N2 * Q.Ttr)^^0.5; 
    }
    override void update_trans_coeffs(GasState Q)
    {
	// The gas is inviscid.
	Q.mu = 0.0;
	Q.k = 0.0;
    }
    override double dudT_const_v(in GasState Q) const
    {
        return _R_N2/(_gamma - 1.0);
    }
    override double dhdT_const_p(in GasState Q) const
    {
	throw new Error("ERROR: VibSpecificNitrogen:dhdT_const_p NOT IMPLEMENTED.");
    }
    override double dpdrho_const_T(in GasState Q) const
    {
	throw new Error("ERROR: VibSpecificNitrogen:dpdrho_const_T NOT IMPLEMENTED.");
    }
    override double gas_constant(in GasState Q) const
    {
	return _R_N2;
    }
    override double internal_energy(in GasState Q) const
    {
	return Q.u + sum(Q.u_modes);
    }
    override double enthalpy(in GasState Q) const
    {
	return Q.u + sum(Q.u_modes) + Q.p/Q.rho;
    }
    override double entropy(in GasState Q) const
    {
	throw new Error("ERROR: VibSpecificNitrogen:entropy NOT IMPLEMENTED.");
    }

private:
    double _R_N2 = 296.805; // gas constant for N2
    double _M_N2 = 0.0280134;
    double _gamma = 7./5.; // ratio of specific heats.
    double kB = Boltzmann_constant;
    double vib_energy(int i) const 
    {
        int I = i+1;
        double w_e = 235857;
        double we_xe = 1432.4;
        double we_ye = -0.226;
        double h = 6.626e-34;
        double c = 2.998e8;
        double e = h*c * (w_e*(I-0.5) - we_xe*(I-0.5)^^2 + we_ye*(I-0.5)^^3);
        return e;
     }
     
     double boltzmann_eq(double Tf1) const
     {
        double summ = 0;
        
        foreach(ei; 0 .. N_VIB_LEVELS) {
             summ += exp(-vib_energy(ei) / (kB*Tf1));
        }        
        double e1 = vib_energy(0);
        double temp_func = (exp(-e1/(kB*Tf1)) / summ);
        return temp_func;   
     }
     
     double compute_Tvib(GasState Q, double x0, double x1, double tol) const
     {
        //this method (secant) is used to compute Tf1
        double init_x0 = x0;
        double init_x1 = x1;
        double fx0 = boltzmann_eq(x0) - Q.massf[0];
        double fx1 = boltzmann_eq(x1) - Q.massf[0];
        double max_it = 100; 
        
        foreach(n; 0 .. max_it) {
            if (abs(fx1) < tol) {return x1;}
        
            double x2 = ((x0*fx1) - (x1*fx0)) / (fx1 - fx0);
            x0 = x1;
            x1 = x2;
            fx0 = fx1;
            fx1 = boltzmann_eq(x2) - Q.massf[0];
        } //end foreach
        return x1;
     }
     
    
     /**
     //BELOW FUNCTION ACTUALLY BELONGS IN THERMOCHEMICALREACTOR CLASS
     // Indices need checking
     double F_coeff(int i, double Temp) const
     {
        double I = i + 1;
        double ft = 10e-6 * Avogadro_number * exp(-3.24093 - (140.69597/Temp^^0.2));
        double delta_t = 0.26679 - (6.99237e5 * Temp) + (4.70073e-9 * Temp^^2);
        double f = (I-1) * ft * exp((I-2)*delta_t);
        return f;     
     }   
     
     double B_coeff(int i, double Temp) const
     {
        double B = F_coeff(i,Temp) * exp(-(vib_energy(i) - vib_energy(i-1)) / (kB * Temp));
        // F_coeff already uses i+1
        return B;
     } 
     **/

} // end class VibSpecificNitrogen

version(vib_specific_nitrogen_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
	auto gm = new VibSpecificNitrogen();
	auto Q = new GasState(N_VIB_LEVELS, N_VIB_LEVELS);
	Q.p = 1.0e5;
	Q.Ttr = 300.0;
	foreach (imode; 0 .. gm.n_modes()) {
	    Q.T_modes[imode] = 300.0;
	}
	Q.massf[] = 1.0/10;

        double R_N2 = 296.805; // gas constant for N2
        double M_N2 = 0.0280134;
        double gamma = 7./5.; // ratio of specific heats.

	gm.update_thermo_from_pT(Q);
	double my_rho = 1.0e5 / (R_N2 * 300.0);
	assert(approxEqual(Q.rho, my_rho, 1.0e-6), failedUnitTest());
	
	double my_u = 2.5 * R_N2 * 300.0;
	assert(approxEqual(Q.u, my_u, 1.0e-6), failedUnitTest());
	
	//double my_u_modes = (Avogadro_number/M_N2) * 1.0 * vib_energy(1); //using mode 1 for test (massf = 1.0 here)
	//assert(approxEqual(Q.u_modes[0], my_u_modes, 1.0e-6), failedUnitTest());
	
	gm.update_trans_coeffs(Q);
	assert(approxEqual(Q.mu, 0.0, 1.0e-6), failedUnitTest());
	assert(approxEqual(Q.k, 0.0, 1.0e-6), failedUnitTest());
	
	gm.update_sound_speed(Q);
	double my_a = (gamma * R_N2 * 300.0)^^0.5;
	assert(approxEqual(Q.a, my_a, 1.0e-6), failedUnitTest());
	//assert(approxEqual(1.0,2.0,1.0e-6), failedUnitTest()); //fail on purpose for check
	
	return 0;
    }
}


