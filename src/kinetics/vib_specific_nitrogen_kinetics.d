/**
 * vib_specific_nitrogen_kinetics.d
 * Authors: Rowan G., Katrina Sklavos, Christopher Chang and Peter J.
 *
 * Edited by Christopher C. to include V-V transfer
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

import std.math;
import nm.complex;
import nm.number;

import gas;
import gas.vib_specific_nitrogen;
import kinetics.thermochemical_reactor;

final class VibSpecificNitrogenRelaxation : ThermochemicalReactor {
    
    this(string fname, GasModel gmodel)
    {
        super(gmodel); // hang on to a reference to the gas model
    }
    
    override void opCall(GasState Q, double tInterval,
                         ref double dtChemSuggest, ref double dtThermSuggest, 
                         ref number[] params)
    {
        int nsteps = 10;
        double dt = tInterval/nsteps;
        number vvRepRate;            // V-V Replenishing rate   
	number vvDepRate;            // V-V Depletion Rate
	number vtRhoDot0, vvRhoDot0; // First level density time derivatives
	number vtRhoDotL, vvRhoDotL; //Last level density time derivatives
	number vtRhoDot, vvRhoDot;   //Intermediate level density time derivatives
	number rho_i, rho_0, rho_l;  //Total density time derivatives

        foreach(step; 0 .. nsteps) {
                
            // Replenishing and depleting equations 
            foreach(imode; 1 .. N_VIB_LEVELS-1) {
                vtRhoDot = (Q.rho*Q.rho/_M_N2)*(BCoeff(imode,Q.T)*Q.massf[imode-1] - 
                                                (FCoeff(imode,Q.T) + BCoeff(imode+1,Q.T)) * Q.massf[imode] + 
                                                FCoeff(imode+1,Q.T)*Q.massf[imode+1]);
                //V-V Exchange Reaction velocities 
                vvRepRate = 0.0; 
                vvDepRate = 0.0;

                foreach (jmode; 2 .. N_VIB_LEVELS) {
                    vvRepRate += (Q.rho*Q.rho/_M_N2)*(vvFCoeff(imode+1,jmode,Q.T)
                                                      *Q.massf[imode+1]*Q.massf[jmode-1]-
                                                      vvBCoeff(imode+1,jmode,Q.T)*Q.massf[imode]*Q.massf[jmode]);

                    vvDepRate += (Q.rho*Q.rho/_M_N2)*(vvFCoeff(imode,jmode,Q.T)
                                                      *Q.massf[imode]*Q.massf[jmode-1]-
                                                      vvBCoeff(imode,jmode,Q.T)*Q.massf[imode-1]*Q.massf[jmode]);
                }

                vvRhoDot = vvRepRate - vvDepRate;
                rho_i = ((vtRhoDot + vvRhoDot) * dt) + Q.massf[imode]*Q.rho;
                Q.massf[imode] = rho_i / Q.rho;
            }


            // Ground state and upper-most state treated specially.
            vtRhoDot0 = (Q.rho*Q.rho/_M_N2)*((-BCoeff(1,Q.T)*Q.massf[0]) 
                                             +(FCoeff(1,Q.T)*Q.massf[1])); 
 
            vvRepRate = 0.0; 

            foreach (jmode; 2 .. N_VIB_LEVELS) {
                vvRepRate += (Q.rho*Q.rho/_M_N2)*(vvFCoeff(1,jmode,Q.T)
                                                  *Q.massf[1]*Q.massf[jmode-1]-
                                                  vvBCoeff(1,jmode,Q.T)*Q.massf[0]*Q.massf[jmode]);
            }

            vvRhoDot0 = vvRepRate;
            rho_0 = ((vtRhoDot0 + vvRhoDot0) * dt) + Q.massf[0]*Q.rho;
            Q.massf[0] = rho_0 / Q.rho;
        
            vtRhoDotL = (Q.rho*Q.rho/_M_N2) * (BCoeff(N_VIB_LEVELS-1,Q.T)
                                               *Q.massf[N_VIB_LEVELS-2]+(-FCoeff(N_VIB_LEVELS-1,Q.T))
                                               *Q.massf[N_VIB_LEVELS-1]); 

            vvDepRate = 0.0; 

            foreach (jmode; 2 .. N_VIB_LEVELS) {
                vvDepRate += (Q.rho*Q.rho/_M_N2)*(vvFCoeff(N_VIB_LEVELS-1,jmode,Q.T)
                                                  *Q.massf[N_VIB_LEVELS-1]*Q.massf[jmode-1]-
                                                  vvBCoeff(N_VIB_LEVELS-1,jmode,Q.T)*Q.massf[N_VIB_LEVELS-2]*Q.massf[jmode]);
            }

            vvRhoDotL = -vvDepRate;
            
            rho_l = ((vtRhoDotL + vvRhoDotL) * dt) + Q.massf[N_VIB_LEVELS-1]*Q.rho;
        
            Q.massf[N_VIB_LEVELS-1] = rho_l / Q.rho;
            scale_mass_fractions(Q.massf, 1.0e-6, 1.0e-3);
        }
        // RJG, CC: Think about this thermo updat.
        // I think it should be a rhou update.
        _gmodel.update_thermo_from_pT(Q);
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
    number FCoeff(int i, number T) const
    {
        double I = i + 1;
        number ft = 1e-6 * Avogadro_number * exp(-3.24093 - (140.69597/T^^0.2));
        number delta_t = 0.26679 - (6.99237e-5 * T) + (4.70073e-9 * T^^2);
        number f = (I-1) * ft * exp((I-2)*delta_t);
        return f;     
    }   
    
    number vvFCoeff(int i, int j, number T) const 
    {
        number vvdelta_t = exp((-6.8/sqrt(T))*abs(i-j));
	number f = 2.5e-20*Avogadro_number*(i-1)*(j-1)*(T/300.0)^^(3.0/2.0)
            *vvdelta_t*(1.5-0.5*vvdelta_t);
	return f;
    }
    
    number BCoeff(int i, number T) const
    {
        number B = FCoeff(i,T) * exp(-(vib_energy(i) - vib_energy(i-1)) / (kB * T));
        // F_coeff already uses i+1
        return B;
    }

    number vvBCoeff(int i, int j, number T) const 
    {
        number B = vvFCoeff(i,j,T)*exp(-1.0*(vib_energy(i)-vib_energy(i-1)
                                             +vib_energy(j-1)-vib_energy(j))/(kB*T));
	return B;
    }

} // end class 
