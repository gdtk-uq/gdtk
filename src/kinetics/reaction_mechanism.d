/**
 * reaction_mechanism.d
 *
 * Author: Rowan G.
 */

module kinetics.reaction_mechanism;

import std.stdio;
import std.algorithm;
import std.math;
import std.conv;
import ntypes.complex;
import nm.number;

import gas;
import util.lua;
import util.lua_service;
import util.msg_service;
import kinetics.rate_constant;
import kinetics.reaction;

class ReactionMechanism
{
public:
    @property @nogc size_t n_reactions() const { return _reactions.length; }
    this(Reaction[] reactions, size_t n_species, double T_lower_limit, double T_upper_limit)
    {
        foreach ( ref r; reactions) {
            _reactions ~= r.dup();
        }
        // Allocate work arrays
        if (!_work_arrays_initialised) {
            _q.length = n_species;
            _L.length = n_species;
            _gibbs_energies.length = n_species;
            _conc_for_source_terms.length = n_species;
            _rates_for_source_terms.length = n_species;
            _work_arrays_initialised = true;
        }
        _T_lower_limit = T_lower_limit;
        _T_upper_limit = T_upper_limit;
    }
    ReactionMechanism dup()
    {
        return new ReactionMechanism(_reactions, _q.length, _T_lower_limit, _T_upper_limit);
    }

    @nogc
    final void eval_rate_constants(GasModel gm, ref GasState Q)
    {

        number T_save = Q.T;
        if ( Q.T < _T_lower_limit ) { Q.T = _T_lower_limit; }
        if ( Q.T > _T_upper_limit ) { Q.T = _T_upper_limit; }
        eval_gibbs_free_energies(gm, Q);

        foreach ( ref r; _reactions ) {
            r.eval_rate_constants(Q, _gibbs_energies);
        }
        // Always reset Q.T on exit
        Q.T = T_save;
    }

    @nogc
    final void eval_tickrates(in number[] conc, number[] forward, number[] backward)
    {
        foreach ( ref r; _reactions ) r.eval_rates(conc);
        foreach ( i, ref r; _reactions ) {
            forward[i] = r.forward_tickrate;
            backward[i] = r.backward_tickrate;
        }
    }
    @nogc
    final void eval_rates(in number[] conc, number[] rates)
    {
        eval_split_rates(conc, _q, _L);
        foreach ( isp; 0..conc.length ) {
            rates[isp] = _q[isp] - _L[isp];
        }
    }
    @nogc
    final void eval_split_rates(in number[] conc, number[] q, number[] L)
    {
        foreach ( ref r; _reactions ) r.eval_rates(conc);
        foreach ( isp; 0..conc.length ) {
            q[isp] = 0.0;
            L[isp] = 0.0;
        }
        foreach ( ir, ref r; _reactions ) {
            foreach ( isp; r.participants ) {
                q[isp] += r.production(isp);
                L[isp] += r.loss(isp);
            }
        }
    }

    @nogc
    final void eval_source_terms(GasModel gmodel, ref GasState Q, ref number[] source)
    {
        gmodel.massf2conc(Q, _conc_for_source_terms);
        eval_rate_constants(gmodel, Q);
        eval_rates(_conc_for_source_terms, _rates_for_source_terms);
        gmodel.rates2source(_rates_for_source_terms, source);
    }

    @nogc
    final number k_f(int ir)
    {
        return _reactions[ir].k_f;
    }
    @nogc
    final number k_b(int ir)
    {
        return _reactions[ir].k_b;
    }
    @nogc
    final number rate(int ir, int isp)
    {
        return _reactions[ir].production(isp) - _reactions[ir].loss(isp);
    }
    @nogc
    final const number production_rate(int ir, int isp)
    {
        return _reactions[ir].production(isp);
    }
    @nogc
    final const number loss_rate(int ir, int isp)
    {
        return _reactions[ir].loss(isp);
    }
    @nogc
    final bool reactionHasParticipant(int ir, int isp)
    {
        foreach (p; _reactions[ir].participants) {
            if (p == isp) {
                return true;
            }
        }
        // If we fall out of the loop, then isp was not found among
        // the participants.
        return false;
    }
    /++
     + Selects an approriate chemistry time step based on rate of change information.
     +
     + This selection algorithm first came to my attention in the paper by Young and Boris.
     + In this implementation, we use the modification suggested by Mott.
     +
     + References:
     + Young and Boris (1977)
     + A Numerical Technique for Solving Stiff Ordinary Differential Equations
     + Associated with the Chemical Kinetics Reactive-Flow Problem
     + The Journal of Physical Chemistry, 81:25 pp. 2424--2427
     +
     + Mott, D.R. (1999)
     + New Quasi-Steady-State and Partial-Equilibrium Methods for
     + Integrating Chemically Reacting Systems.
     + PhD Thesis, The University of Michigan
     +/
    @nogc
    final double estimateStepSize(in number[] conc)
    {
        immutable double MUCH_LARGER_FACTOR = 10000.0;
        immutable double ZERO_EPS = 1.0e-50;
        immutable double ZERO_TOL = 1.0e-30;
        immutable double EPS1 = 1.0e-3;
        immutable double CHEM_STEP_UPPER_INITIAL_LIMIT = 1.0e-10;
        immutable double CHEM_STEP_LOWER_INITIAL_LIMIT = 1.0e-15;

        double min_dt = 1.0e6; // something massive to get us started
        double dt;
        eval_split_rates(conc, _q, _L);
        foreach ( isp, c; conc ) {
            //writefln("c= %e  fabs= %e", c, fabs(_q[isp] - _L[isp]));
            if ( c > 0.0 && (fabs(_q[isp] - _L[isp]) > ZERO_TOL) ) {
                //writefln("q=%e  L=%e", _q[isp], _L[isp]);
                if ( _q[isp] > MUCH_LARGER_FACTOR*_L[isp] ) {
                    // See Mott's thesis, Eq 3.40
                    // In it, he gives this means of estimating the timestep
                    // when q >> L. But he doesn't give a value for how
                    // much larger that factor should be.
                    // In the end, it's only an estimate so it shouldn't matter much.
                    // (hopefully).
                    // The value should be 1/p_i, so that is....
                    dt = c.re/(_L[isp].re + ZERO_EPS);
                    //writefln("Mott's dt= %e ", dt);
                }
                else {
                    dt = fabs( c.re / (_q[isp].re - _L[isp].re) );
                    //writefln("Young and Boris dt= %e", dt);
                }
                min_dt = min(min_dt, dt);
            }
        }
        double dt_chem = EPS1*min_dt;
        //      writefln("dt_chem= %e", dt_chem);
        dt_chem = min(dt_chem, CHEM_STEP_UPPER_INITIAL_LIMIT);
        //      writefln("dt_chem= %e", dt_chem);
        dt_chem = max(dt_chem, CHEM_STEP_LOWER_INITIAL_LIMIT);
        //      writefln("dt_chem= %e", dt_chem);
        return dt_chem;
    }
    @nogc
    final void eval_gibbs_free_energies(GasModel gm, ref GasState Q)
    {
        number[5] T_modes_save;
        number p_save = Q.p;
        version(multi_T_gas) { foreach(imode; 0 .. gm.n_modes) T_modes_save[imode] = Q.T_modes[imode]; }

        Q.p = P_atm;
        version(multi_T_gas) { foreach(imode; 0 .. gm.n_modes) Q.T_modes[imode] = Q.T; }

        gm.gibbs_free_energies(Q, _gibbs_energies);

        Q.p = p_save;
        version(multi_T_gas) { foreach(imode; 0 .. gm.n_modes) Q.T_modes[imode] = T_modes_save[imode]; }
    }

private:
    Reaction[] _reactions;
    // Working array space
    bool _work_arrays_initialised = false;
    number[] _conc_for_source_terms;
    number[] _rates_for_source_terms;
    number[] _gibbs_energies;
    number[] _q;
    number[] _L;
    double _T_lower_limit;
    double _T_upper_limit;
}

/++
 + Creates a ReactionMechanism object from information contained in a LuaTable.
 +
 + The expected format of the table is a of Reaction objects, eg.
 + tab = { [1]={... reaction 1 info ...},
 +         [2]={... reaction 2 info ...},
 +              .......................
 +         [n]={... reaction n info ...}
 +       }
 +/
ReactionMechanism createReactionMechanism(lua_State* L, GasModel gmodel, double T_lower_limit, double T_upper_limit)
{
    auto n_species = gmodel.n_species;
    auto n_reactions = to!int(lua_objlen(L, -1));
    Reaction[] reactions;
    foreach ( i; 1..n_reactions+1 ) {
        lua_rawgeti(L, -1, i);
        reactions ~= createReaction(L, gmodel);
        lua_pop(L, 1);
    }
    return new ReactionMechanism(reactions, n_species, T_lower_limit, T_upper_limit);
}

ReactionMechanism createReactionMechanism(string fname, GasModel gmodel, double T_lower_limit, double T_upper_limit)
{
    auto L = init_lua_State();
    doLuaFile(L, fname);
    lua_getglobal(L, "reaction");
    return createReactionMechanism(L, gmodel, T_lower_limit, T_upper_limit);
}

version(reaction_mechanism_test) {
    import std.math;
    import gas.therm_perf_gas;

    int main() {
        auto gmodel = new ThermallyPerfectGas("sample-input/H2-I2-HI.lua");
        // Test the rate of concentration change at the initial
        // condition for the H2 + I2 reaction system.
        double[] conc = [4.54, 4.54, 0.0];
        auto rc = new ArrheniusRateConstant(1.94e14, 0.0, 20620.0, -1);
        auto gd = GasState(3, 1);
        gd.T = 700.0;
        auto reaction = new ElementaryReaction(rc, rc, gmodel, [0, 1], [1, 1],
                                               [2], [2], 3, -1, -1);
        auto reacMech = new ReactionMechanism([reaction], 3, 100.0, 10000.0);
        double[] rates;
        rates.length = 3;
        reacMech.eval_rate_constants(gmodel, gd);
        reacMech.eval_rates(conc, rates);
        assert(isClose([-643.9303, -643.9303, 1287.8606], rates, 1.0e-6), failedUnitTest());

        return 0;
    }
}
