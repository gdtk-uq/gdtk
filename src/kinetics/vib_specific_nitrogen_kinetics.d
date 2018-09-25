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
 * This particular module deals with the kinetics of vibrational
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
    // Keep a reference to the specific gas model
    // so that we can dip into it for specialized data.
    VibSpecificNitrogen gm;
    
    this(string fname, GasModel gmodel)
    {
        super(gmodel); // hang on to a reference to the gas model
        gm = cast(VibSpecificNitrogen) gmodel;
        assert(gm, "Oops, wrong gas model; should have been VibSpecificNitrogen.");
    }

    @nogc
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
        //
        number rho2M = Q.rho*Q.rho/gm._M_N2;
        foreach(step; 0 .. nsteps) {
            // Replenishing and depleting equations 
            foreach(imode; 1 .. numVibLevels-1) {
                vtRhoDot = rho2M *
                    (BCoeff(imode,Q.T)*Q.massf[imode-1] -
                     (FCoeff(imode,Q.T) + BCoeff(imode+1,Q.T)) * Q.massf[imode] + 
                     FCoeff(imode+1,Q.T)*Q.massf[imode+1]);
                //V-V Exchange Reaction velocities 
                vvRepRate = 0.0; 
                vvDepRate = 0.0;
                foreach (jmode; 2 .. numVibLevels) {
                    vvRepRate += rho2M *
                        (vvFCoeff(imode+1,jmode,Q.T)*Q.massf[imode+1]*Q.massf[jmode-1] -
                         vvBCoeff(imode+1,jmode,Q.T)*Q.massf[imode]*Q.massf[jmode]);

                    vvDepRate += rho2M *
                        (vvFCoeff(imode,jmode,Q.T)*Q.massf[imode]*Q.massf[jmode-1] -
                         vvBCoeff(imode,jmode,Q.T)*Q.massf[imode-1]*Q.massf[jmode]);
                }
                vvRhoDot = vvRepRate - vvDepRate;
                rho_i = ((vtRhoDot + vvRhoDot) * dt) + Q.massf[imode]*Q.rho;
                Q.massf[imode] = rho_i / Q.rho;
            }
            // Ground state and upper-most state treated specially.
            vtRhoDot0 = rho2M *
                ((-BCoeff(1,Q.T)*Q.massf[0]) + (FCoeff(1,Q.T)*Q.massf[1])); 
            vvRepRate = 0.0; 
            foreach (jmode; 2 .. numVibLevels) {
                vvRepRate += rho2M *
                    (vvFCoeff(1,jmode,Q.T)*Q.massf[1]*Q.massf[jmode-1] -
                     vvBCoeff(1,jmode,Q.T)*Q.massf[0]*Q.massf[jmode]);
            }
            vvRhoDot0 = vvRepRate;
            rho_0 = ((vtRhoDot0 + vvRhoDot0) * dt) + Q.massf[0]*Q.rho;
            Q.massf[0] = rho_0 / Q.rho;
            vtRhoDotL = rho2M *
                (BCoeff(numVibLevels-1,Q.T)*Q.massf[numVibLevels-2] +
                 (-FCoeff(numVibLevels-1,Q.T))*Q.massf[numVibLevels-1]); 
            vvDepRate = 0.0; 
            foreach (jmode; 2 .. numVibLevels) {
                vvDepRate += rho2M *
                    (vvFCoeff(numVibLevels-1,jmode,Q.T)*Q.massf[numVibLevels-1]*Q.massf[jmode-1] -
                     vvBCoeff(numVibLevels-1,jmode,Q.T)*Q.massf[numVibLevels-2]*Q.massf[jmode]);
            }
            vvRhoDotL = -vvDepRate;
            rho_l = ((vtRhoDotL + vvRhoDotL) * dt) + Q.massf[numVibLevels-1]*Q.rho;
            Q.massf[numVibLevels-1] = rho_l / Q.rho;
            scale_mass_fractions(Q.massf, 1.0e-6, 1.0e-3);
        }
        // [TODO] RJG, CC: Think about the following thermo update.
        // I think it should be a rhou update.
        _gmodel.update_thermo_from_pT(Q);
        _gmodel.update_sound_speed(Q);
    } // end opCall

     
    // [TODO] Indices need checking
    @nogc
    number FCoeff(int i, number T) const
    {
        double I = i + 1;
        number ft = 1e-6 * Avogadro_number * exp(-3.24093 - (140.69597/T^^0.2));
        number delta_t = 0.26679 - (6.99237e-5 * T) + (4.70073e-9 * T^^2);
        number f = (I-1) * ft * exp((I-2)*delta_t);
        return f;     
    }   

    @nogc
    number vvFCoeff(int i, int j, number T) const 
    {
        number vvdelta_t = exp((-6.8/sqrt(T))*abs(i-j));
	number f = 2.5e-20*Avogadro_number*(i-1)*(j-1)*(T/300.0)^^(3.0/2.0)
            *vvdelta_t*(1.5-0.5*vvdelta_t);
	return f;
    }

    @nogc
    number BCoeff(int i, number T) const
    {
        number B = FCoeff(i,T) * exp(-(gm._vib_energy[i] - gm._vib_energy[i-1]) / (gm.kB * T));
        // F_coeff already uses i+1
        return B;
    }

    @nogc
    number vvBCoeff(int i, int j, number T) const 
    {
        number B = vvFCoeff(i,j,T)*exp(-1.0*(gm._vib_energy[i] - gm._vib_energy[i-1] +
                                             gm._vib_energy[j-1] - gm._vib_energy[j])/(gm.kB*T));
	return B;
    }
} // end class 
