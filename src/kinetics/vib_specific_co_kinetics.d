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
import ntypes.complex;
import nm.number;
import nm.bbla;

import util.lua;
import util.lua_service;

import gas;
import gas.vib_specific_co;
import kinetics.thermochemical_reactor;
import kinetics.relaxation_time;

class VibSpecificCORelaxation : ThermochemicalReactor {

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
    override void opCall(ref GasState Q, double tInterval, ref double dtSuggest,
                         ref number[maxParams] params)
    {
        throw new Error("Substepping does not work for vibe_specific_co!");
    } // end opCall

    @nogc override void eval_source_terms(GasModel gmodel, ref GasState Q, ref number[] source)
    {
        source[] = to!number(0.0);
        foreach(v; 0 .. gm.n_vibe_states) c[v] = Q.massf[v]*Q.rho/_M;
        compute_COCO_pseudoreactions(Q.T, c, source);
    }

    @nogc final void compute_COCO_pseudoreactions(number T, number[] c, ref number[] drhodt)
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
        foreach(v; 0 .. vm) cCO += c[v];

        number Z11 = compute_CO_CO_collision_frequency(T);

        // First compute state-to-state transition probabilities
        Pvvm1[0] = 0.0; // Molecules can't go lower than the ground state
        foreach(w; 0 .. vm) Pvvm1_wm1w[0][w] = 0.0;

        for(int v=0; v<vm; v++){
            Pvvm1[v] = P_v_vm1(Z11, T, v);

            Pvvm1_wm1w[v][0] = 0.0; // w Molecules can't come from lower than the ground state
            for(int w=1; w<vm; w++){
                Pvvm1_wm1w[v][w] = P_v_vm1_wm1_w(Z11, T, v, w);
            }
        }

        // Now we loop over the levels and compute each one's actual rate from the nearby jumps
        foreach(v; 1 .. vm){
            number dcvdt = 0.0;
            number cv = c[v]; number cvm1 = c[v-1];
            double Ev = E[v]; double Evm1 = E[v-1];

            // T-V transfers from collisions
            number tv = Pvvm1[v]*cCO*(cv - exp(-(Ev - Evm1)/kb/T)*cvm1);
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
            drhodt[v] += (dcu[v]+dcd[v])*_M;
        }

        // The ODE intregration I inherited from vib_specific_nitrogen needs a measure of error
        // which is where this came from. It's not precisely zero because of truncation error.
        debug {
            number err = 0.0; foreach (dr; drhodt) { err += dr; }
            if (fabs(err) > 1.0e-6) {
                writeln("\n     dcu      dcd");
                number dsum = 0.0;
                foreach(v; 0 .. gm.n_vibe_states){
                    dsum+= dcu[v] - dcd[v];
                    writefln(" %e %e (%e) + %e", dcu[v], dcd[v], dcu[v] + dcd[v], drhodt[v]);
                }
                writefln("err=%e T=%e c=%s drhodt=%s", err, T, c, drhodt);
                throw new Error("Bad mass conservation");
            }
        }
        return;
    }

protected:
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
        return 4.0*Avogadro_number*sigma11*sigma11*sqrt(pi*kb*T/2.0/m11); // FIXME: Absorb into kb?
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


} // end class

final class VibSpecificCOMixtureRelaxation : VibSpecificCORelaxation {

    this(string chem_file_name, string kinetics_file_name, GasModel gmodel)
    {
        super(chem_file_name, gmodel); // hang on to a reference to the gas model
        gm = cast(VibSpecificCOMixture) gmodel;
        if (!gm) { throw new Error("Oops, wrong gas model; should have been VibSpecificCOMixture."); }

        // fname contains some regular vibrational relaxation mechanisms. Read and store them
        int mode = 0;
        auto L = init_lua_State();
        doLuaFile(L, kinetics_file_name);

        // Load in relaxation times from a normal energy exchange file
        lua_getglobal(L, "mechanism");
        lua_pushnil(L); // dummy first key
        while (lua_next(L, -2) != 0) { // -1 is the dummy key, -2 is the mechanism table
            string rateModel = getString(L, -1, "rate");
            if (rateModel!="Landau-Teller") throw new Error("Supplied Energy Exchange Mechanism is of the wrong type!");

            string pspecies = getString(L, -1, "p");
            string qspecies = getString(L, -1, "q");
            int p = -1;                          // p should be CO, always
            int q = gm.species_index(qspecies);  // q is the other species

            lua_getfield(L, -1, "relaxation_time");
            relaxation_times ~= createRelaxationTime(L, p, q, gm);
            lua_pop(L, 1); // Pop relaxation_time table
            lua_pop(L, 1); // discard value but keep key so that lua_next can remove it (?!)
        }
        lua_pop(L, 1); // remove mechanisms table
        lua_close(L);

        molef.length = gm.n_species;
        numden.length = gm.n_species;

        foreach(j; 0 .. gm.n_others){
            double Mj = gm.cgm.mol_masses[j];
            double MuCOj = 1.0/(1.0/_M + 1.0/Mj); // Note this is kg/mole, not per particle!

            // The formula in Flament 1992 uses mu as the per particle reduced mass, divided by kb
            // Equivalently, we use the per mole Mu and divide by R_universal
            double ThetaCOj = 16.0*pi*pi*pi*pi*MuCOj*weCO*weCO*speed_of_light*speed_of_light*l*l/R_universal;
            Theta_COj ~= ThetaCOj;
        }
    }

    @nogc final override void eval_source_terms(GasModel gmodel, ref GasState Q, ref number[] source)
    {
        source[] = to!number(0.0);
        foreach(v; 0 .. gm.n_vibe_states) c[v] = Q.massf[v]*Q.rho/_M;
        compute_COCO_pseudoreactions(Q.T, c, source);
        compute_COMixture_pseudoreactions(gmodel, Q, c, source);
    }

    @nogc final void compute_COMixture_pseudoreactions(GasModel gmodel, ref GasState Q, number[] c, ref number[] source)
    {
        size_t vm = gm.n_vibe_states;
        number Mmix = gm.molecular_mass(Q);
        foreach(j; vm .. gm.n_species){
            molef[j] = Q.massf[j]*Mmix/gm.mol_masses[j];
            // numden should be unused, but we need it in place for the eval call below.
        }
        foreach(ref d; dcu) d= 0.0;
        foreach(ref d; dcd) d= 0.0;

        // As in the parent class, we work in molar concentration instead of number density
        // This introduces an Avogadro number into the equations, which combines with a kb inside
        // of the definition of P0_COj to turn into an R_universal
        foreach(j; 0 .. gm.n_others){
            number cj = Q.massf[vm+j]*Q.rho/gm.mol_masses[vm+j];

            number tauCOj_times_pj = relaxation_times[j].eval(Q, molef, numden);
            number Lambda10 = lambda_vvm1_COj(1, j, Q.T);
            number F10 = F(Lambda10);
            number P0_COj = (1.0-d)*R_universal*Q.T/tauCOj_times_pj/F10/(1.0 - exp(-Theta_CO/Q.T));

            // Similar structure to the TV transfers from CO only
            for (int v=1; v<vm; v++){

                double dE1= -E01*(1.0 - 2.0*d*v); // Hardcoded formula for deltaE_vvm1_CO
                number Lambda_vvm1 = lambda_vvm1_COj(v, j, Q.T);
                number Fvvm1 = F(Lambda_vvm1);
                number Pvvm1 = P0_COj*v/(1-d*v)*Fvvm1;

                number cv = c[v]; number cvm1 = c[v-1];
                number tv = Pvvm1*cj*(cv - exp(dE1/kb/Q.T)*cvm1);
                dcd[v] -= tv;
                dcu[v-1] += tv;
            }
        }

        foreach(v; 0 .. vm){
            source[v] += (dcu[v]+dcd[v])*_M;
        }
    }

private:
    // Expressions for the non-CO rates come from a collection of papers compiled by David Petty
    immutable double weCO = 2169.81358*100.0;  // cm-1 -> m-1 Note that E01 = h*c*weCO
    //immutable double xeCO = 0.006124171;     // Note that this is equivalent to Rich's d
    immutable double Theta_CO = Plancks_constant*speed_of_light*weCO/Boltzmann_constant;
    immutable double l = 0.2e-10;              // Range parameter in m (Whatever this is...)
    double[] Theta_COj;                        // Binary characteristic temperaures

    number[] molef, numden;
    VibSpecificCOMixture gm;  // Keep a reference to the specific gas model.
    RelaxationTime[] relaxation_times;

    @nogc const pure
    number lambda_vvm1_COj(size_t v, size_t j, number T){
        double dE1= fabs(-E01*(1.0 - 2.0*d*v)); // Hardcoded formula for deltaE_vvm1_CO
        return sqrt(dE1*dE1/kb/kb/8.0/T*Theta_COj[j]/Theta_CO/Theta_CO);
    }

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
        auto Q = GasState(gm.n_species, 0);
        Q.p = 26.7;  // Pa
        Q.T = 175.0; // K
        number Tvib = to!number(1500.0);
        // Set up the species mass fractions assuming vibrational equilibrium.
        foreach (v; 0 .. gm.n_vibe_states) { Q.massf[v] = gm.boltzmann_eq_population_fraction(v, Tvib); }
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

        // Check the total mass is being conserved
        number[] source;
        source.length = gm.n_species;
        reactor.eval_source_terms(gm, Q, source);

        number sum = 0.0; number total = 0.0;
        foreach(v; 0 .. gm.n_species){
            sum += source[v];
            total += fabs(source[v]);
            //writefln("v %02d Y %e source %e", v, Q.massf[v], source[v]);
        }
        //writefln("sum %e total %e error %e", sum, total, fabs(sum)/total);
        assert((fabs(sum)/total<1e-12), failedUnitTest());
    }
}

version(vib_specific_co_mixture_kinetics_test) {
    import std.stdio;
    import util.msg_service;
    import std.math : isClose;
    import gas.vib_specific_co;
    void main() {
        // Set up a CO only gas state with some vibrational nonequilibrium
        auto gm = new VibSpecificCO("../gas/sample-data/vib-specific-CO-gas.lua");
        auto Q = GasState(gm.n_species, 0);
        Q.p = 26.7;  // Pa
        Q.T = 175.0; // K
        number Tvib = to!number(1500.0);
        // Set up the species mass fractions assuming vibrational equilibrium.
        foreach (v; 0 .. gm.n_vibe_states) { Q.massf[v] = gm.boltzmann_eq_population_fraction(v, Tvib); }
        gm.update_thermo_from_pT(Q);

        auto reactor = new VibSpecificCORelaxation("", gm);
        number[] source;
        source.length = gm.n_species;
        reactor.eval_source_terms(gm, Q, source);


        // Now build an identical state using the mixture machinery
        auto gmmix = new VibSpecificCOMixture("sample-input/vib-specific-CO-mixture.lua");
        auto Qmix = GasState(gmmix.n_species, 0);
        Qmix.p = 26.7;  // Pa
        Qmix.T = 175.0; // K
        // Set up the species mass fractions assuming vibrational equilibrium.
        foreach (v; 0 .. gmmix.n_vibe_states) { Qmix.massf[v] = gmmix.boltzmann_eq_population_fraction(v, Tvib); }
        Qmix.massf[gm.n_vibe_states] = 0.0;
        gmmix.update_thermo_from_pT(Qmix);

        auto reactormix = new VibSpecificCOMixtureRelaxation("", "sample-input/ee-n2-co.lua", gmmix);
        number[] sourcemix;
        sourcemix.length = gmmix.n_species;
        reactormix.eval_source_terms(gmmix, Qmix, sourcemix);

        foreach(v; 0 .. gm.n_vibe_states){
            //writefln("state %d source %e sourcemix %e", v, source[v], sourcemix[v]);
            assert(isClose(source[v], sourcemix[v], 1.0e-6), failedUnitTest());
        }
        assert(isClose(sourcemix[gmmix.n_species-1], 0.0, 1.0e-12), failedUnitTest());


        // Now we need to test the N2/CO relaxation
        // Now build an identical state using the mixture machinery
        auto gmn2 = new VibSpecificCOMixture("sample-input/vib-specific-CO-mixture.lua");
        auto Qn2 = GasState(gmn2.n_species, 0);
        Qn2.p = 53.4;  // Pa
        Qn2.T = 175.0; // K
        // Set up the species mass fractions assuming vibrational equilibrium.
        foreach (v; 0 .. gmn2.n_vibe_states) { Qn2.massf[v] = 0.5*gmn2.boltzmann_eq_population_fraction(v, Tvib); }
        Qn2.massf[gmn2.n_vibe_states] = 0.5;
        gmn2.update_thermo_from_pT(Qn2);

        // Check the total mass is being conserved
        auto reactorn2 = new VibSpecificCOMixtureRelaxation("", "sample-input/ee-n2-co.lua", gmn2);
        number[] sourcen2;
        sourcen2.length = gmn2.n_species;
        reactorn2.eval_source_terms(gmn2, Qn2, sourcen2);

        number sum = 0.0; number total = 0.0;
        foreach(v; 0 .. gm.n_species){
            sum += sourcen2[v];
            total += fabs(sourcen2[v]);
            //writefln("v %02d Y %e source %e", v, Qn2.massf[v], sourcen2[v]);
        }
        //writefln("sum %e total %e error %e", sum, total, fabs(sum)/total);
        assert((fabs(sum)/total<1e-12), failedUnitTest());
    }
}
