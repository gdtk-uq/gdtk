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
    number[] dRhoDt; // Time derivative of the species densities, kg/s/m^^3.

    this(string fname, GasModel gmodel)
    {
        super(gmodel); // hang on to a reference to the gas model
        gm = cast(VibSpecificNitrogen) gmodel;
        if (!gm) { throw new Error("Oops, wrong gas model; should have been VibSpecificNitrogen."); }
        dRhoDt.length = gm.numVibLevels;
    }

    @nogc
    override void opCall(GasState Q, double tInterval,
                         ref double dtChemSuggest, ref double dtThermSuggest,
                         ref number[maxParams] params)
    {
        // First, set a time step size that does not cause too large a change in species.
        number rhoErr = computeDrhoDt(Q, dRhoDt);
        double bigDrhoDt = 0.0;
        foreach (rd; dRhoDt) {
            if (fabs(rd) > bigDrhoDt) { bigDrhoDt = fabs(rd); }
        }
        double dt = (dtChemSuggest > 0.0) ? dtChemSuggest : tInterval;
        if (bigDrhoDt > 0.0) {
            dt = 0.001 * Q.rho/bigDrhoDt;
        }
        int nsteps = cast(int)(ceil(tInterval/dt));
        dt = tInterval/nsteps;
        // Now, do the time integration.
        // [TODO] We should upgrade the integration step to have better order of error.
        foreach(step; 0 .. nsteps) {
            foreach (i; 0 .. gm.numVibLevels) { Q.massf[i] += dRhoDt[i]/Q.rho * dt; }
            scale_mass_fractions(Q.massf, 1.0e-6, 1.0e-3);
            _gmodel.update_thermo_from_rhou(Q);
            _gmodel.update_sound_speed(Q);
            if (step < nsteps-1) {
                // For the next step, if there is one.
                rhoErr = computeDrhoDt(Q, dRhoDt);
            }
        }
        dtChemSuggest = dt;
    } // end opCall

    @nogc override void eval_source_terms(GasModel gmodel, GasState Q, ref number[] source)
    {
        number rhoErr = computeDrhoDt(Q, source);
    }

    @nogc number computeDrhoDt(GasState Q, ref number[] drhodt)
    // Computes the array of time derivatives for the species densities.
    // Returns the sum of those derivatives as a measure of error.
    {
        number vvRepRate;            // V-V Replenishing rate
	number vvDepRate;            // V-V Depletion Rate
	number vtRhoDot0, vvRhoDot0; // First level density time derivatives
	number vtRhoDotL, vvRhoDotL; // Last level density time derivatives
	number vtRhoDot, vvRhoDot;   // Intermediate level density time derivatives
        //
        number rho2M = Q.rho*Q.rho/gm._M_N2;
        int L = gm.numVibLevels;
        //
        // Replenishing and depleting equations for all but the first and last species.
        foreach(i; 1 .. L-1) {
            vtRhoDot = rho2M * (BCoeff(i,Q.T)*Q.massf[i-1]
                                - (FCoeff(i,Q.T) + BCoeff(i+1,Q.T)) * Q.massf[i]
                                + FCoeff(i+1,Q.T)*Q.massf[i+1]);
            //V-V Exchange Reaction velocities
            vvRepRate = 0.0;
            vvDepRate = 0.0;
            foreach (j; 1 .. L) {
                vvRepRate += rho2M * (vvFCoeff(i+1,j,Q.T)*Q.massf[i+1]*Q.massf[j-1]
                                      - vvBCoeff(i+1,j,Q.T)*Q.massf[i]*Q.massf[j]);
                //
                vvDepRate += rho2M * (vvFCoeff(i,j,Q.T)*Q.massf[i]*Q.massf[j-1]
                                      - vvBCoeff(i,j,Q.T)*Q.massf[i-1]*Q.massf[j]);
            }
            vvRhoDot = vvRepRate - vvDepRate;
            drhodt[i] = vtRhoDot + vvRhoDot;
        }
        // Ground state and upper-most state treated specially.
        vtRhoDot0 = rho2M * ((-BCoeff(1,Q.T)*Q.massf[0]) + (FCoeff(1,Q.T)*Q.massf[1]));
        vvRepRate = 0.0;
        foreach (j; 1 .. L) {
            vvRepRate += rho2M * (vvFCoeff(1,j,Q.T)*Q.massf[1]*Q.massf[j-1]
                                  - vvBCoeff(1,j,Q.T)*Q.massf[0]*Q.massf[j]);
        }
        vvRhoDot0 = vvRepRate;
        drhodt[0] = vtRhoDot0 + vvRhoDot0;
        //
        vtRhoDotL = rho2M * (BCoeff(L-1,Q.T)*Q.massf[L-2] + (-FCoeff(L-1,Q.T))*Q.massf[L-1]);
        vvDepRate = 0.0;
        foreach (j; 1 .. L) {
            vvDepRate += rho2M * (vvFCoeff(L-1,j,Q.T)*Q.massf[L-1]*Q.massf[j-1]
                                  - vvBCoeff(L-1,j,Q.T)*Q.massf[L-2]*Q.massf[j]);
        }
        vvRhoDotL = -vvDepRate;
        drhodt[L-1] = vtRhoDotL + vvRhoDotL;
        //
        // As a measure of error, return the sum of the density derivatives.
        // For a good calculation, it should be zero.
        number err = 0.0; foreach (dr; drhodt) { err += dr; }
        return err;
    } // end computeDrhoDt()

    @nogc
    number FCoeff(int i, number T) const
    {
        number ft = 1e-6 * Avogadro_number * exp(-3.24093 - (140.69597/T^^0.2));
        number delta_t = 0.26679 - (6.99237e-5 * T) + (4.70073e-9 * T^^2);
        number f = i * ft * exp((i-1)*delta_t);
        return f;
    }

    @nogc
    number vvFCoeff(int i, int j, number T) const
    {
        number vvdelta_t = exp((-6.8/sqrt(T))*abs(i-j));
	number f = 2.5e-20*Avogadro_number*i*j*(T/300.0)^^(3.0/2.0)
            *vvdelta_t*(1.5-0.5*vvdelta_t);
	return f;
    }

    // [TODO] It might be good to inline these functions because these functions
    // will be repeating the forward-reaction coefficient evaluation.

    @nogc
    number BCoeff(int i, number T) const
    {
        number B = FCoeff(i,T) * exp(-(gm._vib_energy[i] - gm._vib_energy[i-1]) / (gm.kB * T));
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
