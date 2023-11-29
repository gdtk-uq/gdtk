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
import std.algorithm;
import ntypes.complex;
import nm.number;
import nm.bbla;

import gas;
import gas.vib_specific_nitrogen;
import kinetics.thermochemical_reactor;

final class VibSpecificNitrogenRelaxation : ThermochemicalReactor {

    this(string fname, GasModel gmodel)
    {
        super(gmodel); // hang on to a reference to the gas model
        gm = cast(VibSpecificNitrogen) gmodel;
        if (!gm) { throw new Error("Oops, wrong gas model; should have been VibSpecificNitrogen."); }
        int L = gm.numVibLevels;
        dRhoDt0.length = L; dRhoDt1.length = L;
        mf0.length = L; mf1.length = L;
        crhs = new Matrix!(double)(L, L+1);
        vtF.length = L; vtB.length = L;
        vvF.length = L; foreach (i; 0 .. L) { vvF[i].length = L; }
        vvB.length = L; foreach (i; 0 .. L) { vvB[i].length = L; }
    }

    @nogc
    override void opCall(ref GasState Q, double tInterval, ref double dtSuggest,
                         ref number[maxParams] params)
    {
        int L = gm.numVibLevels;
        foreach (i; 0 .. L) { mf0[i] = Q.massf[i]; }
        // Set a time step size.
        // Note that, presently, we ignore dtSuggest.
        number rhoErr = computeDrhoDt(Q.rho, Q.T, mf0, dRhoDt0);
        // Limit the stepsize to allow only small changes in mass fractions per step
        // by looking at just the ground-state population with the largest reaction rate.
        // This might avoid the minor species from driving the time step to tiny values.
        double dt = tInterval;
        number bigDrhoDt = 0.0;
        foreach (rd; dRhoDt0) {
            if (fabs(rd) > bigDrhoDt) { bigDrhoDt = fabs(rd); }
        }
        if (bigDrhoDt > 0.0) {
            dt = min(0.001*Q.rho.re*Q.massf[0].re/bigDrhoDt.re, tInterval);
        }
        int nsteps = cast(int)(ceil(tInterval/dt));
        dt = tInterval/nsteps;
        // Now, do the time integration.
        // [TODO] We should upgrade the integration step to have better order of error.
        // This might allow us to take bigger steps and save some computational effort.
        foreach(step; 0 .. nsteps) {
            // Do an implicit (backward-Euler) update of the mass-fraction vector.
            //
            // Assemble the linear system for computing the mass-fraction increments,
            // solve for the mass-fraction increments and add them.
            double h = 1.0e-6;
            foreach (i; 0 .. L) {
                // Perturbed state vector, changing one element for this pass..
                foreach (j; 0 .. L) { mf1[j] = mf0[j]; if (i == j) { mf1[j] += h; } }
                // Sensitivity coefficients computed via finite-differences.
                rhoErr = computeDrhoDt(Q.rho, Q.T, mf1, dRhoDt1);
                foreach (j; 0 .. L) {
                    number dFdy = (dRhoDt1[j]-dRhoDt0[j])/Q.rho/h;
                    crhs[i,j] = ((i == j) ? 1.0/dt : 0.0) - dFdy.re;
                }
                // Right-hand side vector is also packed into the augmented matrix.
                crhs[i,L] = dRhoDt0[i].re/Q.rho.re;
            }
            gaussJordanElimination!double(crhs);
            foreach (i; 0 .. L) {
                mf1[i] = mf0[i] + crhs[i,L]; if (mf1[i] < 0.0) { mf1[i] = 0.0; }
            }
            scale_mass_fractions(mf1, 1.0e-6, 1.0e-3);
            //
            foreach (i; 0 .. L) { Q.massf[i] = mf1[i]; }
            _gmodel.update_thermo_from_rhou(Q);
            _gmodel.update_sound_speed(Q);
            //
            if (step < nsteps-1) {
                // For the next step, if there is one.
                foreach (i; 0 .. L) { mf0[i] = mf1[i]; }
                rhoErr = computeDrhoDt(Q.rho, Q.T, mf0, dRhoDt0);
            }
        }
        dtSuggest = dt;
    } // end opCall

    @nogc override void eval_source_terms(GasModel gmodel, ref GasState Q, ref number[] source)
    {
        number rhoErr = computeDrhoDt(Q.rho, Q.T, Q.massf, source);
    }

    @nogc number computeDrhoDt(number rho, number T, number[] massf, ref number[] drhodt)
    // Computes the array of time derivatives for the species densities.
    // Returns the sum of those derivatives as a measure of error.
    {
        number rho2M = rho*rho/gm._M_N2;
        int L = gm.numVibLevels;
        //
        // V-T exchange reactions.
        //
        // Start by computing the arrays of forward and backward coefficients.
        foreach (i; 1 .. L) {
            // Forward rate coefficient for reaction level i to level i-1, equation 29
            number ft = 1e-6 * Avogadro_number * exp(-3.24093 - (140.69597/T^^0.2));
            number delta_t = 0.26679 - (6.99237e-5 * T) + (4.70073e-9 * T^^2);
            vtF[i] = i * ft * exp((i-1)*delta_t);
        }
        foreach (i; 1 .. L) {
            // Backward rate coefficient for reaction level i-1 to level i, equation 28.
            vtB[i] = vtF[i] * exp(-(gm._vib_energy[i] - gm._vib_energy[i-1]) / (gm.kB * T));
        }
        // V-T reaction velocities, equations 26 and 27.
        foreach(i; 1 .. L-1) {
            drhodt[i] = rho2M * (vtB[i]*massf[i-1]
                                - (vtF[i] + vtB[i+1])*massf[i]
                                + vtF[i+1]*massf[i+1]);
        }
        // Ground state and upper-most state treated specially.
        drhodt[0] = rho2M * (-vtB[1]*massf[0] + vtF[1]*massf[1]);
        drhodt[L-1] = rho2M * (vtB[L-1]*massf[L-2] - vtF[L-1]*massf[L-1]);
        //
        // V-V exchange reactions.
        //
        foreach (i; 1 .. L) {
            foreach (j; 1 .. L) {
                // F^(j-1,j)_(i,i-1) Forward rate coefficient, equation 41.
                number vvdelta_t = exp((-6.8/sqrt(T))*abs(i-j));
                vvF[i][j] = 2.5e-20*Avogadro_number*(i-1)*(j-1)
                    * (T/300.0)^^(3.0/2.0) * vvdelta_t * (1.5-0.5*vvdelta_t);
            }
        }
        foreach (i; 1 .. L) {
            foreach (j; 1 .. L) {
                // B^(j,j-1)_(i-1,i) Backward rate coefficient, equation 40.
                vvB[i][j] = vvF[i][j]
                    * exp(-1.0*(gm._vib_energy[i] - gm._vib_energy[i-1] +
                                gm._vib_energy[j-1] - gm._vib_energy[j])/(gm.kB*T));
            }
        }
        //
        // V-V Exchange Reaction velocities
        //
        number vvRepRate; // V-V Replenishing rate
	number vvDepRate; // V-V Depletion Rate
        //
        foreach(i; 1 .. L-1) {
            vvRepRate = 0.0;
            vvDepRate = 0.0;
            foreach (j; 1 .. L) {
                vvRepRate += vvF[i+1][j]*massf[i+1]*massf[j-1] - vvB[i+1][j]*massf[i]*massf[j];
                vvDepRate += vvF[i][j]*massf[i]*massf[j-1] - vvB[i][j]*massf[i-1]*massf[j];
            }
            drhodt[i] += rho2M * (vvRepRate - vvDepRate);
        }
        vvRepRate = 0.0;
        foreach (j; 1 .. L) {
            vvRepRate += vvF[1][j]*massf[1]*massf[j-1] - vvB[1][j]*massf[0]*massf[j];
        }
        drhodt[0] += rho2M * vvRepRate;
        //
        vvDepRate = 0.0;
        foreach (j; 1 .. L) {
            vvDepRate += vvF[L-1][j]*massf[L-1]*massf[j-1] - vvB[L-1][j]*massf[L-2]*massf[j];
        }
        drhodt[L-1] -= rho2M * vvDepRate;
        //
        // As a measure of error, return the sum of the density derivatives.
        // For a good calculation, it should be zero.
        number err = 0.0; foreach (dr; drhodt) { err += dr; }
        debug {
            import std.stdio;
            if (err > 1.0e-9) {
                writefln("err=%e, rho=%e T=%e massf=%s drhodt=%s", err, rho, T, massf, drhodt);
            }
        }
        return err;
    } // end computeDrhoDt()

private:
    VibSpecificNitrogen gm;  // Keep a reference to the specific gas model.
    number[] mf0, mf1;  // Mass-fraction arrays
    number[] dRhoDt0, dRhoDt1;  // Time derivatives of the species densities, kg/s/m^^3.
    Matrix!double crhs;  // Augmented matrix for the linear solve in the implicit-update.
    number[] vtF, vtB;  // Coefficients for the vibrational-translation exchanges.
    number[][] vvF, vvB;  // Coefficients for the vib-vib exchanges.
} // end class
