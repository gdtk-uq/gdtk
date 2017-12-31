/**
 * vib_specific_nitrogen_kinetics.d
 * Authors: Rowan G., Katrina Sklavos and Peter J.
 *
 * This is a 10-level vibrationally-specific model for nitrogen
 * as descrbied in:
 *
 * Giordano, et al. (1997)
 * Vibrationally Relaxing Flow of N2 past an Infinite Cylinder
 * Journal of Thermophysics and Heat Transfer, vol 11, no 1, pp 27 - 35
 *
 *
 * This particular module details with the kinetics of vibrational
 * population exchanges.
 *
 */

module kinetics.vib_specific_nitrogen_kinetics;

import std.math : exp;

import gas;
import gas.vib_specific_nitrogen;
import kinetics.thermochemical_reactor;

final class VibSpecificNitrogenRelaxtion : ThermochemicalReactor {
    
    this(string fname, GasModel gmodel)
    {
	super(gmodel); // hang on to a reference to the gas model
    }
    
    override void opCall(GasState Q, double tInterval,
			 ref double dtChemSuggest, ref double dtThermSuggest, 
			 ref double[] params)
    {
        int nsteps = 10;
        double dt = tInterval/nsteps;
        foreach(step; 0 .. nsteps) {
        
        
        // Replenishing and depleting equations 
        foreach(imode; 1 .. N_VIB_LEVELS-1) {
            double rho_dot = (Q.rho*Q.rho/_M_N2) * (B_coeff(imode,Q.T)*Q.massf[imode-1] +
                                (-F_coeff(imode, Q.T) - B_coeff(imode+1, Q.T))*Q.massf[imode] 
                                 + F_coeff(imode+1, Q.T)*Q.massf[imode+1]);
            
           
            double rho_i = (rho_dot * dt) + Q.massf[imode]*Q.rho;
            Q.massf[imode] = rho_i / Q.rho;
        }
        //Hard coding the replenishing equation (ground state)
        double rho_dot_0 = (Q.rho*Q.rho/_M_N2) * ((-B_coeff(1, Q.T)*Q.massf[0]) 
                                 + (F_coeff(1, Q.T)*Q.massf[1])); 
        double rho_0 = (rho_dot_0 * dt) + Q.massf[0]*Q.rho;
        Q.massf[0] = rho_0 / Q.rho;
        
        
        //Hard coding the depleting equation (last quantum level)
        double rho_dot_l = (Q.rho*Q.rho/_M_N2) * (B_coeff(N_VIB_LEVELS-1, Q.T)*Q.massf[N_VIB_LEVELS-1-1] +
                                (-F_coeff(N_VIB_LEVELS-1, Q.T))*Q.massf[N_VIB_LEVELS-1]); 
                                             
        double rho_l = (rho_dot_l * dt) + Q.massf[N_VIB_LEVELS-1]*Q.rho;
        Q.massf[N_VIB_LEVELS-1] = rho_l / Q.rho;
        //writeln("Gas state: ");
        //writeln(Q);
        scale_mass_fractions(Q.massf, 1.0e-6, 1.0e-3);
        }
        
        //_gmodel.update_thermo_from_rhop(Q);
        _gmodel.update_thermo_from_pT(Q);
        //assert(fabs(1.0 - sum(Q.massf)) < 1e-6, "Oops, total mass fraction is not 1.0");
	_gmodel.update_sound_speed(Q);
        
    } // end opCall

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
     
     // Indices need checking
     double F_coeff(int i, double Temp) const
     {
        double I = i + 1;
        double ft = 1e-6 * Avogadro_number * exp(-3.24093 - (140.69597/Temp^^0.2));
        double delta_t = 0.26679 - (6.99237e-5 * Temp) + (4.70073e-9 * Temp^^2);
        double f = (I-1) * ft * exp((I-2)*delta_t);
        return f;     
     }   
     
     double B_coeff(int i, double Temp) const
     {
        double B = F_coeff(i,Temp) * exp(-(vib_energy(i) - vib_energy(i-1)) / (kB * Temp));
        // F_coeff already uses i+1
        return B;
     }

} // end class 
