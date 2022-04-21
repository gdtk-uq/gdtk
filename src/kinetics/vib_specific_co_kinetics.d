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

        dcu.length = L;
        dcd.length = L;
        c.length = L;
        E.length = L;
        Pvvm1.length = L;
        Pvvm1_wm1w.length = L;
        foreach(ref P; Pvvm1_wm1w) P.length = L;
        for(int v=0; v<gm.n_vibe_states; v++) E[v] = vibe_energy(v);
    }

    @nogc
    override void opCall(GasState Q, double tInterval, ref double dtSuggest,
                         ref number[maxParams] params)
    {
        throw new Error("Substepping does not work for vibe_specific_co!");
    } // end opCall

    @nogc override void eval_source_terms(GasModel gmodel, GasState Q, ref number[] source)
    {
        number rhoErr = computeDrhoDt(Q.rho, Q.T, Q.massf, source);
    }

    @nogc number computeDrhoDt(number rho, number T, number[] massf, ref number[] drhodt)
    /*
        Compute the "reaction" source terms that drive changes in the CO vibrational population.
        This is essentially equation (2.1) but converted to mole concentration instead of
        number density. I did this to reduce the large amount of floating point error in the
        number density based formulation.

        Notes: There are some old fashioned for loops here because we don't wave to use
        ulong's (which are basically what size_t is) in any math. We always want to use
        int's when there's math involved, or else your can get underflow issues, like
        I did in the early days of this file.

        @author: Nick Gibbons
    */
    {
        size_t vm = gm.n_vibe_states;
        foreach(ref d; dcu) d= 0.0;
        foreach(ref d; dcd) d= 0.0;

        number cCO = 0.0;
        foreach(v; 0 .. vm){
            c[v] = massf[v]*rho/_M;
            cCO += c[v];
        }

        number Z11 = compute_CO_CO_collision_frequency(T);

        // Molecules can't go lower than the ground state
        Pvvm1[0] = 0.0;
        foreach(w; 0 .. vm) Pvvm1_wm1w[0][w] = 0.0;

        for(int v=0; v<vm; v++){
            Pvvm1[v] = P_v_vm1(Z11, T, v);

            Pvvm1_wm1w[v][0] = 0.0; // w Molecules can't come from lower than the ground state
            for(int w=1; w<vm; w++){
                Pvvm1_wm1w[v][w] = P_v_vm1_wm1_w(Z11, T, v, w);
            }
        }

        // Now we loop over the levels and compute each ones rate from the nearby jumps
        foreach(v; 1 .. vm){
            number dcvdt = 0.0;
            number cv = c[v]; number cvm1 = c[v-1];
            double Ev = E[v]; double Evm1 = E[v-1];

            // T-V transfers from collisions
            number tv = Pvvm1[v]*cCO*(cv   - exp(-(Ev - Evm1)/kb/T)*cvm1);
            dcd[v] -= tv;
            dcu[v-1] += tv;

            // V-V transfers are a bit a of tricky accounting problem. Essentially we need to loop
            // over each unique combination of collisions. You can think about the upper triangle
            // of the matrix, excluding the diagonal row, which has no net effect on dcdt
            // This is NOT properly explained in the paper, for some reason.
            foreach(w; v+1 .. vm){
                number cw = c[w]; number cwm1 = c[w-1];
                double Ew = E[w]; double Ewm1 = E[w-1];
                number vv = Pvvm1_wm1w[v][w]*(cv*cwm1 - exp(-(Ev+Ewm1-Ew-Evm1)/kb/T)*cvm1*cw);
                dcd[v] -= vv;
                dcu[v-1] += vv;
                dcd[w] += vv;
                dcu[w-1] -= vv;
            }
        }
        foreach(v; 0 .. vm){
            drhodt[v] = (dcu[v]+dcd[v])*_M;
        }

        // The ODE intregration I inherited from vib_specific_nitrogen needs a measure of error
        // so I do what they did. Shouldn't this be zero? No, there's actually a fair bit of
        // rounding error.
        number err = 0.0; foreach (dr; drhodt) { err += dr; }
        debug {
            if (fabs(err) > 1.0e-6) {
                writeln("\n     dcu      dcd");
                number dsum = 0.0;
                foreach(v; 0 .. gm.n_vibe_states){
                    dsum+= dcu[v] - dcd[v];
                    writefln(" %e %e (%e) + %e", dcu[v], dcd[v], dcu[v] + dcd[v], drhodt[v]);
                }
                writefln("err=%e, rho=%e T=%e massf=%s drhodt=%s", err, rho, T, massf, drhodt);
                throw new Error("Bad mass conservation");
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
        /*
           Collision rate has been modified to be per mole from the formula
           in the paper, so that we can work in mole concentrations and
           have smaller numbers with less round off error
        */
        return 4.0*Avogadro_number*sigma11*sigma11*sqrt(pi*kb*T/2.0/m11);
    }

    @nogc const number x_v_vm1(number T, int v){
        return half_to_3_halves*sqrt(theta_dash11/T)*(1.0-2.0*d*v);
    }

    @nogc const number y_v_vm1_wm1_w(number T, int v, int w){
        return 2.0*d*half_to_3_halves*sqrt(theta_dash11/T)*fabs(1.0*(v-w));
    }

    @nogc const number F(number x){
        number exp_minus2xon3 = exp(-2.0*x/3.0);
        return 0.5*(3.0-exp_minus2xon3)*exp_minus2xon3;
    }

    @nogc const number P_v_vm1(number Z11, number T, int v){
        // Equation 2.3, for CO V-T exchange between v to v-1
        number xvvm1 = x_v_vm1(T, v);
        number Fxvvm1 = F(xvvm1);
        number P = Z11*P11*T*v/(1.0-d*v)*Fxvvm1;
        return P;
    }

    @nogc const number P_v_vm1_wm1_w(number Z11, number T, int v, int w){
        // Equation 2.5, for CO V-V exchange between v to v-1 and w-1 to w
        number yvvm1_wm1w = y_v_vm1_wm1_w(T, v, w);
        number Fyvvm1_wm1w = F(yvvm1_wm1w);
        number P = Z11*Q11*T*v/(1.0-d*v)*w/(1.0-d*w)*Fyvvm1_wm1w;
        return P;
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
    number[] dcu, dcd;
    number[] c, Pvvm1;
    number[][] Pvvm1_wm1w;

} // end class

version(vib_specific_co_kinetics_test) {
    import std.stdio;
    import util.msg_service;
    import std.math : isClose;
    import gas.vib_specific_co;
    void main() {
        // Shout out to the author for providing this wonderful table of example rates.
        immutable size_t nrows = 12;
        int[nrows] vs = [1,8,9,10,11,12,15,20,30,40,50,60];
        int[nrows] ws = [0,7,8,9,10,11,14,19,29,39,49,59];
        number[nrows] Pvvm1_times_c_targets = [8.728e-12, 1.695e-9, 3.010e-9, 5.274e-9, 9.152e-9,
                             1.575e-8, 7.735e-8, 1.009e-6, 1.455e-4, 1.875e-2, 2.279e0, 2.679e2];
        number[nrows] Pvvm1_wwp1_times_c_targets = [1.644e3, 1.146e5, 1.470e5, 1.837e5, 2.251e5,
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
        number cCO = 0.0;
        foreach(v; 0 .. gm.n_vibe_states) cCO += Q.massf[v]*Q.rho/_M;

        auto reactor = new VibSpecificCORelaxation("", gm);
        number Z11 = reactor.compute_CO_CO_collision_frequency(Q.T);

        foreach(i; 0 .. nrows){
            int v=vs[i];
            int w=ws[i];
            number Pvvm1_times_c_target = Pvvm1_times_c_targets[i];
            number Pvvm1_wwp1_times_c_target = Pvvm1_wwp1_times_c_targets[i];

            number Pvvm1_times_c = reactor.P_v_vm1(Z11, Q.T, v)*cCO;
            number Pvvm1_wwp1_times_c = reactor.P_v_vm1_wm1_w(Z11, Q.T, v, w+1)*cCO;

            //writefln("%02d:  Pv %e (%e) diff: %e", i, Pvvm1_times_c, Pvvm1_times_c_target, (Pvvm1_times_c-Pvvm1_times_c_target)/Pvvm1_times_c_target);
            //writefln("%02d: Pvw %e (%e) diff: %e", i, Pvvm1_wwp1_times_c, Pvvm1_wwp1_times_c_target, (Pvvm1_wwp1_times_c-Pvvm1_wwp1_times_c_target)/Pvvm1_wwp1_times_c_target);

            assert(isClose(Pvvm1_times_c, Pvvm1_times_c_target, 5.0e-3), failedUnitTest());
            assert(isClose(Pvvm1_wwp1_times_c, Pvvm1_wwp1_times_c_target, 5.0e-3), failedUnitTest());
        }
    }
}
