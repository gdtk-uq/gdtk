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

import gas;
import util.lua;
import util.lua_service;
import util.msg_service;
import kinetics.rate_constant;
import kinetics.reaction;

class ReactionMechanism
{
public:
    @property ulong n_reactions() const { return _reactions.length; }
    this(in Reaction[] reactions, size_t n_species)
    {
	foreach ( ref r; reactions) {
	    _reactions ~= r.dup();
	}
	// Allocate work arrays
	if (!_work_arrays_initialised) {
	    _q.length = n_species;
	    _L.length = n_species;
	    _work_arrays_initialised = true;
	}
    }
    ReactionMechanism dup() const
    {
	return new ReactionMechanism(_reactions, _q.length);
    }
    
    final void eval_rate_constants(in GasState Q)
    {
	foreach ( ref r; _reactions ) r.eval_rate_constants(Q);
    }

    final void eval_rates(in double[] conc, double[] rates)
    {
	eval_split_rates(conc, _q, _L);
	foreach ( isp; 0..conc.length ) rates[isp] = _q[isp] - _L[isp];
    }
    final void eval_split_rates(in double[] conc, double[] q, double[] L)
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
    final double k_f(int ir)
    {
	return _reactions[ir].k_f;
    }
    final double k_b(int ir)
    {
	return _reactions[ir].k_b;
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
    final double estimateStepSize(in double[] conc)
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
		    dt = c/(_L[isp] + ZERO_EPS);
		    //writefln("Mott's dt= %e ", dt);
		}
		else {
		    dt = fabs( c / (_q[isp] - _L[isp]) );
		    //writefln("Young and Boris dt= %e", dt);
		}
		min_dt = min(min_dt, dt);
	    }
	}
	double dt_chem = EPS1*min_dt;
	//	writefln("dt_chem= %e", dt_chem);
	dt_chem = min(dt_chem, CHEM_STEP_UPPER_INITIAL_LIMIT);
	//	writefln("dt_chem= %e", dt_chem);
	dt_chem = max(dt_chem, CHEM_STEP_LOWER_INITIAL_LIMIT);
	//	writefln("dt_chem= %e", dt_chem);
	return dt_chem;
    }

private:
    Reaction[] _reactions;
    // Working array space
    static bool _work_arrays_initialised = false;
    static double[] _q;
    static double[] _L;
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
ReactionMechanism createReactionMechanism(lua_State* L, in GasModel gmodel)
{
    auto n_species = gmodel.n_species;
    auto n_reactions = to!int(lua_objlen(L, -1));
    Reaction[] reactions;
    foreach ( i; 1..n_reactions+1 ) {
	lua_rawgeti(L, -1, i);
	reactions ~= createReaction(L, n_species);
	lua_pop(L, 1);
    }
    return new ReactionMechanism(reactions, n_species);
}

ReactionMechanism createReactionMechanism(string fname, in GasModel gmodel)
{
    auto L = init_lua_State(fname);
    lua_getglobal(L, "reaction");
    return createReactionMechanism(L, gmodel);
}

unittest
{
    import std.math;
    // Test the rate of concentration change at the initial
    // condition for the H2 + I2 reaction system.
    double[] conc = [4.54, 4.54, 0.0];
    auto rc = new ArrheniusRateConstant(1.94e14, 0.0, 20620.0);
    auto gd = new GasState(3, 1);
    gd.T[0] = 700.0;
    auto reaction = new ElementaryReaction(rc, rc, [0, 1], [1, 1],
					   [2], [2], 3);
    auto reacMech = new ReactionMechanism([reaction], 3);
    double[] rates;
    rates.length = 3;
    reacMech.eval_rate_constants(gd);
    reacMech.eval_rates(conc, rates);
    assert(approxEqual([-643.9303, -643.9303, 1287.8606], rates), failedUnitTest());
}
