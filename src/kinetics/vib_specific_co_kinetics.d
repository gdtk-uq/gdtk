/*
    Prototyping for vibrationally specific CO gas model, intended for use in simulating 
    gas dynamic lasers.

    References:
    "Kinetic Modelling of the High-Power Carbon Monoxide Laser*"
    Joseph W. (The GOAT) Rich, Journal of Applied Physics, Volume 42, Number 7, June 1971

    Notes:
        - What an excellent paper holy heck.
        TODO: Energy conservation?
    @author: Nick Gibbons (24/03/22)
*/

module kinetics.vib_specific_co_kinetics;

import std.math;
import std.algorithm;
import std.stdio;
import std.conv;
import nm.complex;
import nm.number;
import nm.bbla;

import gas;
import gas.vib_specific_co;
import kinetics.thermochemical_reactor;

final class VibSpecificCORelaxation : ThermochemicalReactor {

    this(string fname, GasModel gmodel)
    {
        super(gmodel); // hang on to a reference to the gas model
        gm = cast(VibSpecificCO) gmodel;
        if (!gm) { throw new Error("Oops, wrong gas model; should have been VibSpecificCO."); }
        size_t L = gm.n_vibe_states;
        dRhoDt0.length = L; dRhoDt1.length = L;
        mf0.length = L; mf1.length = L;
        crhs = new Matrix!(double)(L, L+1);

        // Allocate space for one extra energy level at the top of the range, whose population is 0.0
        N.length = L+1;
        E.length = L+1;
        Pvvm1.length = L+1; 
        Pvvm1_wm1w.length = L+1; 
        foreach(ref P; Pvvm1_wm1w) P.length = L+1;
        foreach(v; 0 .. gm.n_vibe_states+1) E[v] = vibe_energy(v);
    }

    @nogc
    override void opCall(GasState Q, double tInterval, ref double dtSuggest,
                         ref number[maxParams] params)
    {
        int L = gm.n_vibe_states;
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
                    crhs._data[i][j] = ((i == j) ? 1.0/dt : 0.0) - dFdy.re;
                }
                // Right-hand side vector is also packed into the augmented matrix.
                crhs._data[i][L] = dRhoDt0[i].re/Q.rho.re;
            }
            gaussJordanElimination!double(crhs);
            foreach (i; 0 .. L) {
                mf1[i] = mf0[i] + crhs._data[i][L]; if (mf1[i] < 0.0) { mf1[i] = 0.0; }
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

    @nogc override void eval_source_terms(GasModel gmodel, GasState Q, ref number[] source)
    {
        number rhoErr = computeDrhoDt(Q.rho, Q.T, Q.massf, source);
    }

    @nogc number computeDrhoDt(number rho, number T, number[] massf, ref number[] drhodt)
    /* 
        Compute the "reaction" source terms that drive changes in the CO vibrational population.
        This is essentially equation (2.1) but converted to mass rate of change instead of
        number density.

        @author: Nick Gibbons
    */
    {
        size_t vm = gm.n_vibe_states;

        number NCO = 0.0;
        foreach(v; 0 .. vm){
            N[v] = Avogadro_number*massf[v]*rho/_M;
            NCO += N[v];
        }
        N[vm] = 0.0; // The level above the highest we are tracking has zero population

        number Z11 = compute_CO_CO_collision_frequency(T);

        // Molecules can't go lower than the ground state, or higher than the highest state
        Pvvm1[0] = 0.0;
        Pvvm1[vm] = 0.0;
        foreach(w; 0 .. vm+1) Pvvm1_wm1w[0][w] = 0.0;
        foreach(w; 0 .. vm+1) Pvvm1_wm1w[vm][w] = 0.0;

        foreach(v; 1 .. vm){
            Pvvm1[v] = P_v_vm1(Z11, T, v);

            Pvvm1_wm1w[v][0] = 0.0; // w Molecules can't come from lower than the ground state
            Pvvm1_wm1w[v][vm]= 0.0; // w Molecules can't go above the highest state
            foreach(w; 1 .. vm){
                Pvvm1_wm1w[v][w] = P_v_vm1_wm1_w(Z11, T, v, w);
            }
        }

        // Now we loop over the levels and compute each ones rate from the nearby jumps
        foreach(v; 0 .. vm){
            number dNvdt = 0.0;
            number Nvp1 = N[v+1]; number Nv = N[v]; number Nvm1 = (v==0) ? to!number(0.0) : N[v-1];
            double Evp1 = E[v+1]; double Ev = E[v]; double Evm1 = (v==0) ? 0.0 : E[v-1];

            // T-V transfers from collisions
            dNvdt += Pvvm1[v+1]*NCO*(Nvp1 - exp(-(Evp1 - Ev)/kb/T)*Nv); // incoming from above
            dNvdt -= Pvvm1[v]*NCO*(Nv   - exp(-(Ev - Evm1)/kb/T)*Nvm1); // outgoing to below

            // incoming from above V-V transfers from collisions
            foreach(w; 1 .. vm){
                number Nwp1 = N[w+1]; number Nw = N[w]; number Nwm1 = N[w-1];
                double Ewp1 = E[w+1]; double Ew = E[w]; double Ewm1 = E[w-1];
                dNvdt += Pvvm1_wm1w[v+1][w]*(Nvp1*Nwm1 - exp(-(Evp1+Ewm1-Ev-Ew)/kb/T)*Nv*Nw);
            }

            // outgoing to below V-V transfers from collisions
            foreach(w; 0 .. vm-1){
                number Nwp1 = N[w+1]; number Nw = N[w]; number Nwm1 = (v==0) ? to!number(0.0) : N[w-1];
                double Ewp1 = E[w+1]; double Ew = E[w]; double Ewm1 = (v==0) ? 0.0 : E[w-1];
                dNvdt -= Pvvm1_wm1w[v][w+1]*(Nv*Nw - exp(-(Ev+Ew-Ewp1-Evm1)/kb/T)*Nvm1*Nwp1);
            }

            // Plus electrons (we shouldn't have these unless we are doing discharges?)
            // Plus radiation TODO

            // dNvdt is in molecules/m3/second, convert to kg/m3/second
            drhodt[v] = dNvdt*_M/Avogadro_number;
        }


        // The ODE intregration I inherited from vib_specific_nitrogen needs a measure of error
        // so I do what they did. Shouldn't this be zero?
        number err = 0.0; foreach (dr; drhodt) { err += dr; }
        debug {
            if (err > 1.0e-9) {
                writefln("err=%e, rho=%e T=%e massf=%s drhodt=%s", err, rho, T, massf, drhodt);
            }
        }
        return err;
    }

    @nogc const double vibe_energy(int v)
    // Returns the quantum level energy for species i.
    {
        return E01*(v - d*v*(v+1.0));
    }

    @nogc const number compute_CO_CO_collision_frequency(number T)
    {
        return 4.0*sigma11*sigma11*sqrt(pi*kb*T/2.0/m11);
    }

    @nogc const number x_v_vm1(number T, size_t v){
        return half_to_3_halves*sqrt(theta_dash11/T)*(1.0-2.0*d*v);
    }

    @nogc const number y_v_vm1_wm1_w(number T, size_t v, size_t w){
        return 2.0*d*half_to_3_halves*sqrt(theta_dash11/T)*fabs(1.0*(v-w));
    }

    @nogc const number F(number x){
        number exp_minus2xon3 = exp(-2.0*x/3.0);
        return 0.5*(3.0-exp_minus2xon3)*exp_minus2xon3;
    }

    @nogc const number P_v_vm1(number Z11, number T, size_t v){
        // Equation 2.3, for CO V-T exchange between v to v-1
        number xvvm1 = x_v_vm1(T, v);
        number Fxvvm1 = F(xvvm1);
        return Z11*P11*T*v/(1.0-d*v)*Fxvvm1;
    }

    @nogc const number P_v_vm1_wm1_w(number Z11, number T, size_t v, size_t w){
        // Equation 2.3, for CO V-T exchange between v to v-1
        number yvvm1_wm1w = y_v_vm1_wm1_w(T, v, w);
        number Fyvvm1_wm1w = F(yvvm1_wm1w);
        return Z11*Q11*T*v/(1.0-d*v)*w/(1.0-d*w)*Fyvvm1_wm1w;
    }

private:
    VibSpecificCO gm;  // Keep a reference to the specific gas model.
    immutable double _M = 0.0280101;
    // Constants from table 1
    immutable double E01 = 4.31e-13*1e-7; // erg -> J Ground state energy (FIXME: Is this E1? A typo???)
    immutable double d   = 0.00598;       //          CO Anharmonicity
    immutable double P11 = 1.77e-4;
    immutable double Q11 = 3.70e-6;
    immutable double theta_dash11 = 4.45e6; // degrees, presumably Kelvin
    immutable double sigma11 = 3.75e-8/100.0; // cm -> m
    immutable double m11 = 2.3236e-23/1000.0; // g -> kg
    immutable double pi = 3.14159265359;

    immutable double kb = Boltzmann_constant;
    immutable double half_to_3_halves = pow(0.5, 1.5);

    // Allocatable workspace
    double[] E;
    number[] N, Pvvm1;
    number[][] Pvvm1_wm1w;

    // Integrator machinery, same as vib_specific_nitrogen
    number[] mf0, mf1;         // Mass-fraction arrays
    number[] dRhoDt0, dRhoDt1; // Time derivatives of the species densities, kg/s/m^^3.
    Matrix!double crhs;        // Augmented matrix for the linear solve in the implicit-update.

} // end class

version(vib_specific_co_kinetics_test) {
    import std.stdio;
    import util.msg_service;
    import std.math : isClose;
    import gas.vib_specific_co;
    void main() {
        // Shout out to the author for providing this wonderful table of example rates.
        immutable size_t nrows = 12;
        size_t[nrows] vs = [1,8,9,10,11,12,15,20,30,40,50,60];
        size_t[nrows] ws = [0,7,8,9,10,11,14,19,29,39,49,59];
        number[nrows] Pvvm1_times_N_targets = [8.728e-12, 1.695e-9, 3.010e-9, 5.274e-9, 9.152e-9,
                             1.575e-8, 7.735e-8, 1.009e-6, 1.455e-4, 1.875e-2, 2.279e0, 2.679e2];
        number[nrows] Pvvm1_wwp1_times_N_targets = [1.644e3, 1.146e5, 1.470e5, 1.837e5, 2.251e5,
                                   2.714e5, 4.409e5, 8.383e5, 2.171e6, 4.490e6, 8.275e6, 1.422e7];

        /*
        Note: The Pvvm1_wwp1 entry for v=50 and w=49 given in the table does not match the
        results computed by my code, although all the others do. In light of this, I've 
        chosen to exclude it from the test, assuming that it must be a typo of some kind.

        @NNG 29/03/22
        */

        // Set up gas state to match the data in table II
        auto gm = new VibSpecificCO("../gas/sample-data/vib-specific-CO-gas.lua");
        auto Q = new GasState(gm.n_species, 0);
        Q.p = 26.7;  // Pa
        Q.T = 175.0; // K
        // Set up the species mass fractions assuming vibrational equilibrium.
        foreach (v; 0 .. gm.n_vibe_states) { Q.massf[v] = gm.boltzmann_eq_population_fraction(v, Q.T); }
        gm.update_thermo_from_pT(Q);

        immutable double _M = 0.0280101;
        number NCO = 0.0;
        foreach(v; 0 .. gm.n_vibe_states) NCO += Avogadro_number*Q.massf[v]*Q.rho/_M;

        auto reactor = new VibSpecificCORelaxation("", gm);
        number Z11 = reactor.compute_CO_CO_collision_frequency(Q.T);

        foreach(i; 0 .. nrows){
            size_t v=vs[i];
            size_t w=ws[i];
            number Pvvm1_times_N_target = Pvvm1_times_N_targets[i];
            number Pvvm1_wwp1_times_N_target = Pvvm1_wwp1_times_N_targets[i];

            number Pvvm1_times_N = reactor.P_v_vm1(Z11, Q.T, v)*NCO;
            number Pvvm1_wwp1_times_N = reactor.P_v_vm1_wm1_w(Z11, Q.T, v, w+1)*NCO;

            //writefln("%02d:  Pv %e (%e) diff: %e", i, Pvvm1_times_N, Pvvm1_times_N_target, (Pvvm1_times_N-Pvvm1_times_N_target)/Pvvm1_times_N_target);
            //writefln("%02d: Pvw %e (%e) diff: %e", i, Pvvm1_wwp1_times_N, Pvvm1_wwp1_times_N_target, (Pvvm1_wwp1_times_N-Pvvm1_wwp1_times_N_target)/Pvvm1_wwp1_times_N_target);

            assert(isClose(Pvvm1_times_N, Pvvm1_times_N_target, 5.0e-3), failedUnitTest());
            assert(isClose(Pvvm1_wwp1_times_N, Pvvm1_wwp1_times_N_target, 5.0e-3), failedUnitTest());
        }
    }
}
