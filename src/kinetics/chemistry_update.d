/**
 * This module provides routines to update a gas phase
 * chemistry system.
 *
 * Author: Rowan G. and Kyle D.
 */

module kinetics.chemistry_update;

import std.conv;
import std.stdio;
import std.string;
import std.math;
import std.algorithm;
import ntypes.complex;
import nm.number;
import util.lua;
import util.lua_service;
import gas;
import kinetics.thermochemical_reactor;
import kinetics.reaction_mechanism;

immutable double DT_INCREASE_PERCENT = 10.0; // allowable percentage increase on succesful step
immutable double DT_DECREASE_PERCENT = 50.0; // allowable percentage decrease on succesful step
                                             // Yes, you read that right. Sometimes a step is succesful
                                             // but the timestep selection algorithm will suggest
                                             // a reduction. We limit that reduction to no more than 50%.
immutable double DT_REDUCTION_FACTOR = 10.0; // factor by which to reduce timestep
                                             // after a failed attempt
immutable double H_MIN = 1.0e-15; // Minimum allowable step size
immutable double ALLOWABLE_MASSF_ERROR = 1.0e-3; // Maximum allowable error in mass fraction over an update

enum ResultOfStep { success, failure };

final class ChemistryUpdate : ThermochemicalReactor {
    ReactionMechanism rmech;
    ChemODEStep cstep;
    bool tightTempCoupling;
    int maxSubcycles;
    int maxAttempts;

    this(string fname, GasModel gmodel)
    {
        super(gmodel);
        // Allocate memory
        _Qinit = GasState(gmodel.n_species, gmodel.n_modes);
        _conc0.length = gmodel.n_species;
        _concOut.length = gmodel.n_species;

        // Configure other parameters via Lua state.
        auto L = init_lua_State();
        doLuaFile(L, fname);

        // Check on species order before proceeding.
        lua_getglobal(L, "species");
        if (lua_isnil(L, -1)) {
            string errMsg = format("There is no species listing in your chemistry input file: '%s'\n", fname);
            errMsg ~= "Please be sure to use an up-to-date version of prep-chem to generate your chemistry file.\n";
            throw new Error(errMsg);
        }
        foreach (isp; 0 .. gmodel.n_species) {
            lua_rawgeti(L, -1, isp);
            auto sp = to!string(luaL_checkstring(L, -1));
            if (sp != gmodel.species_name(isp)) {
                string errMsg = "Species order is incompatible between gas model and chemistry input.\n";
                errMsg ~= format("Chemistry input file is: '%s'.\n", fname);
                errMsg ~= format("In gas model: index %d ==> species '%s'\n", isp, gmodel.species_name(isp));
                errMsg ~= format("In chemistry: index %d ==> species '%s'\n", isp, sp);
                throw new Error(errMsg);
            }
            lua_pop(L, 1);
        }
        lua_pop(L, 1);

        lua_getglobal(L, "config");
        lua_getfield(L, -1, "tempLimits");
        lua_getfield(L, -1, "lower");
        double T_lower_limit = lua_tonumber(L, -1);
        lua_pop(L, 1);
        lua_getfield(L, -1, "upper");
        double T_upper_limit = lua_tonumber(L, -1);
        lua_pop(L, 1);
        lua_pop(L, 1);
        lua_pop(L, 1);

        lua_getglobal(L, "reaction");
        rmech = createReactionMechanism(L, gmodel, T_lower_limit, T_upper_limit);
        lua_pop(L, 1);

        lua_getglobal(L, "config");
        lua_getfield(L, -1, "odeStep");
        lua_getfield(L, -1, "method");
        string ode_method = to!string(lua_tostring(L, -1));
        lua_pop(L, 1);
        switch (ode_method) {
        case "rkf":
             lua_getfield(L, -1, "errTol");
             double err_tol = lua_tonumber(L, -1);
             lua_pop(L, 1);
             cstep = new RKFStep(gmodel, rmech, err_tol);
             break;
        case "alpha-qss":
             lua_getfield(L, -1, "eps1");
             double eps1 = lua_tonumber(L, -1);
             lua_pop(L, 1);
             lua_getfield(L, -1, "eps2");
             double eps2 = lua_tonumber(L, -1);
             lua_pop(L, 1);
             lua_getfield(L, -1, "delta");
             double delta = lua_tonumber(L, -1);
             lua_pop(L, 1);
             lua_getfield(L, -1, "maxIters");
             int max_iters = to!int(lua_tointeger(L, -1));
             lua_pop(L, 1);
             cstep = new AlphaQssStep(gmodel, rmech, eps1, eps2, delta, max_iters);
             break;
        default:
             string errMsg = format("ERROR: The odeStep '%s' is unknown.\n", ode_method);
             throw new Error(errMsg);
        }
        lua_pop(L, 1); // pops 'odeStep'

        lua_getfield(L, -1, "tightTempCoupling");
        tightTempCoupling = lua_toboolean(L, -1);
        lua_pop(L, 1);

        lua_getfield(L, -1, "maxSubcycles");
        maxSubcycles = to!int(luaL_checkinteger(L, -1));
        lua_pop(L, 1);

        lua_getfield(L, -1, "maxAttempts");
        maxAttempts = to!int(luaL_checkinteger(L, -1));
        lua_pop(L, 1);

        lua_pop(L, 1); // pops 'config'
        lua_close(L);
    }

    @nogc
    override void opCall(ref GasState Q, double tInterval, ref double dtSuggest,
                         ref number[maxParams] params)
    {
        _Qinit.copy_values_from(Q);
        _gmodel.massf2conc(Q, _conc0);

        // 0. Evaluate the rate constants.
        //    It helps to have these computed before doing other setup work.
        rmech.eval_rate_constants(_gmodel, Q);
        // 1. Sort out the time step for possible subcycling.
        double t = 0.0;
        double h;
        double dtSave;
        if ( dtSuggest > tInterval )
            h = tInterval;
        else if ( dtSuggest <= 0.0 )
            h = rmech.estimateStepSize(_conc0);
        else
            h = dtSuggest;
        //
        // 2. Now do the interesting stuff, increment species change
        int cycle = 0;
        int attempt = 0;
        for ( ; cycle < maxSubcycles; ++cycle ) {
            /* Copy the last good timestep suggestion before moving on.
             * If we're near the end of the interval, we want
             * the penultimate step. In other words, we don't
             * want to store the fractional step of (tInterval - t)
             * that is taken as the last step.
             */
            dtSave = dtSuggest;
            h = min(h, tInterval - t);
            attempt= 0;
            for ( ; attempt < maxAttempts; ++attempt ) {
                ResultOfStep result = cstep(_conc0, h, _concOut, dtSuggest);
                bool passesMassFractionTest = true;
                if ( result  == ResultOfStep.success ) {
                    /* We succesfully took a step of size h according to the ODE method.
                     * However, we need to test that that mass fractions have remained ok.
                     */
                    _gmodel.conc2massf(_concOut, Q);
                    auto massfTotal = sum(Q.massf);
                    if ( fabs(massfTotal - 1.0) > ALLOWABLE_MASSF_ERROR ) {
                        passesMassFractionTest = false;
                    }
                }

                if ( result == ResultOfStep.success && passesMassFractionTest ) {
                    t += h;
                    // Copy succesful values in concOut to conc0, ready for next step
                    _conc0[] = _concOut[];
                    /* We can now make some decision about how to
                     * increase the timestep. We will take some
                     * guidance from the ODE step, but also check
                     * that it's sensible. For example, we won't
                     * increase by more than 10% (or INCREASE_PERCENT).
                     * We'll also balk if the ODE step wants to reduce
                     * the stepsize on a successful step. This is because
                     * if the step was successful there shouldn't be any
                     * need (stability wise or accuracy related) to warrant
                     * a reduction. Thus if the step is successful and
                     * the dtSuggest comes back lower, let's just set
                     * h as the original value for the successful step.
                     */
                    double hMax = h*(1.0 + DT_INCREASE_PERCENT/100.0);
                    if (dtSuggest > h) {
                        /* It's possible that dtSuggest is less than h.
                         * When that occurs, we've made the decision to ignore
                         * the new timestep suggestion supplied by the ODE step.
                         * Our reasoning is that we'd like push for an aggressive
                         * timestep size. If that fails, then a later part of
                         * the timestepping algorithm will catch that and reduce
                         * the timestep.
                         *
                         * It's also possible that the suggested timestep
                         * is much larger than what we just used. For that
                         * case, we cap it at hMax. That's the following line.
                         */
                        h = min(dtSuggest, hMax);
                    }
                    break;
                }
                else { // in the case of failure...
                    /* We now need to make some decision about
                     * what timestep to attempt next. We follow
                     * David Mott's suggestion in his thesis (on p. 51)
                     * and reduce the timestep by a factor of 2 or 3.
                     * (The actual value is set as DT_REDUCTION_FACTOR).
                     * In fact, for the types of problems we deal with
                     * I have found that a reduction by a factor of 10
                     * is more effective.
                     */
                    h /= DT_REDUCTION_FACTOR;
                    if ( h < H_MIN ) {
                        string errMsg = "Hit the minimum allowable timestep in chemistry update.";
                        debug { errMsg ~= format("\ndt= %.4e", H_MIN); }
                        Q.copy_values_from(_Qinit);
                        throw new ThermochemicalReactorUpdateException(errMsg);
                    }
                }
            } // end attempts at single subcycle.
            if ( attempt == maxAttempts ) {
                string errMsg = "Hit maximum number of step attempts within a subcycle for chemistry update.";
                // We did poorly. Let's put the original GasState back in place,
                // and let the outside world know via an Exception.
                Q.copy_values_from(_Qinit);
                throw new ThermochemicalReactorUpdateException(errMsg);
            }
            /* Otherwise, we've done well.
             * If tight temperature coupling is requested, we can reevaluate
             * the temperature at this point. With regards to tight coupling,
             * we follow Oran and Boris's advice on pp. 140-141.
             * To paraphrase: solving a separate differential equation for
             * temperature is computationally expensive, however, it usually
             * suffices to reevaluate the temperature assuming that total internal
             * energy of the system has not changed but has been redistributed
             * among chemical states. Since the chemistry has not altered much,
             * the iteration for temperature should converge in one or two
             * iterations.
             *
             * My own additional argument for avoiding a temperature differential
             * equation is that it does not play nicely with the special
             * ODE methods for chemistry that exploit structure in the species
             * differential equations.
             */
            if ( tightTempCoupling ) {
                _gmodel.conc2massf(_conc0, Q);
                _gmodel.update_thermo_from_rhou(Q);
                rmech.eval_rate_constants(_gmodel, Q);
            }

            if ( t >= tInterval ) { // We've done enough work.
                // If we've only take one cycle, then we would like to use
                // dtSuggest rather than using the dtSuggest from the
                // penultimate step.
                if (cycle == 0) {
                    dtSave = dtSuggest;
                }
                break;
            }
        }
        if ( cycle == maxSubcycles ) {
            // If we make it out here, then we didn't complete within
            // the maximum number of subscyles... we are taking too long.
            // Let's return the gas state to how it was before we failed
            // and throw an Exception.
            string errMsg = "Hit maximum number of subcycles while attempting chemistry update.";
            Q.copy_values_from(_Qinit);
            throw new ThermochemicalReactorUpdateException(errMsg);
        }
        // At this point, it appears that everything has gone well.
        // We'll tweak the mass fractions in case they are a little off.
        _gmodel.conc2massf(_concOut, Q);
        auto massfTotal = sum(Q.massf);
        foreach (ref mf; Q.massf) mf /= massfTotal;
        /*
         * RJG 2018-11-16
         * Let's remove this piece of complication.
         * Let's have the chemistry module *only* update mass fractions.
         * The caller can decide on the thermodynamic constraint.
         * This allows re-use as a chemistry updater in multi-temperature models.
         *
        if (_gmodel.n_modes >= 1) {
            // Changing mass fractions does not change internal energies
            // but it might have been redistributed.
            try {
                _gmodel.update_thermo_from_rhoT(Q);
            }
            catch (GasModelException err) {
                string msg = "The call to update_thermo_from_rhoT in the chemistry update failed.";
                debug { msg ~= format("\ncaught %s", err.msg); }
            throw new ThermochemicalReactorUpdateException(msg);
            }
            auto uOther = sum(Q.u_modes);
            Q.u = uTotal - uOther;
            try {
                _gmodel.update_thermo_from_rhou(Q);
            }
            catch (GasModelException err) {
                string msg = "The call to update_thermo_from_rhou in the chemistry update failed.";
                debug { msg ~= format("\ncaught %s", err.msg); }
            throw new ThermochemicalReactorUpdateException(msg);
            }
        }
        */
        dtSuggest = dtSave;
    }

    @nogc override void eval_source_terms(GasModel gmodel, ref GasState Q, ref number[] source) {
        rmech.eval_source_terms(gmodel, Q, source);
    }

private:
    // Some memory workspace
    GasState _Qinit;
    number[] _conc0;
    number[] _concOut;
} // end class ChemistryUpdate

/++
 + ChemODEStep is an abstract class that provides the interface
 + for a chemistry ODE update step. This class provides just one
 + public service so we make that happen when the object is
 + called like a function with opCall().
 +
 + You might rightly ask "Why use a class when a function would do?"
 + It turns out that there is a lot of extra state that needs to be
 + stored and a lot of working array space. This is possible as
 + module level static data. However, it just plain looks neater
 + in a class.
 +/

class ChemODEStep
{
public:
    this(ReactionMechanism rmech)
    {
        _rmech = rmech;
    }
    @nogc
    abstract ResultOfStep opCall(number[] y0, double h, number[] yOut, ref double hSuggest);
private:
    ReactionMechanism _rmech;
}

/++
 + The Runge-Kutta-Fehlberg method, specialised to work
 + with a chemical kinetic ODE.
 +
 + References:
 + Fehlberg, E. (1969)
 + Low-order classical Runge-Kutta formulas with stepsize control
 + and their application to some heat transfer problems.
 + NASA Technical Report, R-315
 +
 + Cash, J. R. and Karp, A. H. (1990)
 + A Variable Order Runge-Kutta Method for Initial Value
 + Problems with Rapidly Varying Right-Hand Sides
 + ACM Transactions on Mathematical Software, 16:3, pp. 201--222
 +
 + Press, W. H., Teukolsky, S. A., Vetterling, W. T. and Flannery, B. P. (2007)
 + Numerical Recipes: The Art of Scientific Computing, Third Edition
 + Cambridge University Press, New York, USA
 +/

class RKFStep : ChemODEStep {
public:
    this(in GasModel gmodel, ReactionMechanism rmech, double tolerance)
    {
        super(rmech);
        _tolerance = tolerance;
        // Allocate working arrays
        _ndim = gmodel.n_species;
        _yTmp.length = _ndim;
        _yErr.length = _ndim;
        _k1.length = _ndim;
        _k2.length = _ndim;
        _k3.length = _ndim;
        _k4.length = _ndim;
        _k5.length = _ndim;
        _k6.length = _ndim;
    }

    @nogc
    override ResultOfStep opCall(number[] y0, double h, number[] yOut, ref double hSuggest)
    {
        // 0. Set up the constants associated with the update formula
        // We place them in here for visibility reasons... no one else
        // needs to be using the coefficients a21, a31, etc.
        immutable double a21=1./5., a31=3./40., a32=9./40.,
            a41=3./10., a42=-9./10., a43 = 6./5.,
            a51=-11./54., a52=5./2., a53=-70./27., a54=35./27.,
            a61=1631./55296., a62=175./512., a63=575./13824., a64=44275./110592., a65=253./4096.;
        immutable double b51=37./378., b53=250./621., b54=125./594., b56=512./1771.,
            b41=2825./27648., b43=18575./48384., b44=13525./55296., b45=277./14336., b46=1.0/4.0;

        // 1. Apply the formula to evaluate intermediate points
        _rmech.eval_rates(y0, _k1);
        foreach (i; 0 .. _yTmp.length) { _yTmp[i] = y0[i] + h*(a21*_k1[i]); }

        _rmech.eval_rates(_yTmp, _k2);
        foreach (i; 0 .. _yTmp.length) { _yTmp[i] = y0[i] + h*(a31*_k1[i] + a32*_k2[i]); }

        _rmech.eval_rates(_yTmp, _k3);
        foreach (i; 0 .. _yTmp.length) { _yTmp[i] = y0[i] + h*(a41*_k1[i] + a42*_k2[i] + a43*_k3[i]); }

        _rmech.eval_rates(_yTmp, _k4);
        foreach (i; 0 .. _yTmp.length) { _yTmp[i] = y0[i] + h*(a51*_k1[i] + a52*_k2[i] + a53*_k3[i] + a54*_k4[i]); }

        _rmech.eval_rates(_yTmp, _k5);
        foreach (i; 0 .. _yTmp.length) { _yTmp[i] = y0[i] + h*(a61*_k1[i] + a62*_k2[i] + a63*_k3[i] + a64*_k4[i] + a65*_k5[i]); }

        _rmech.eval_rates(_yTmp, _k6);

        // 2. Compute new value and error esimate
        foreach (i; 0 .. yOut.length) { yOut[i] = y0[i] + h*(b51*_k1[i] + b53*_k3[i] + b54*_k4[i] + b56*_k6[i]); }
        foreach (i; 0 .. _yErr.length) { _yErr[i] = yOut[i] - (y0[i] + h*(b41*_k1[i] + b43*_k3[i] + b44*_k4[i] + b45*_k5[i] + b46*_k6[i])); }

        // 3. Lastly, use error estimate as a means to suggest
        //    a new timestep.
        //    Compute error using tol as atol and rtol as suggested in
        //    Press et al. (2007)
        double err = 0.0;
        double sk = 0.0;
        double atol = _tolerance;
        double rtol = _tolerance;

        foreach ( i; 0.._ndim ) {
            sk = atol + rtol*max(fabs(y0[i].re), fabs(yOut[i].re));
            err += (_yErr[i].re/sk)*(_yErr[i].re/sk);
        }
        err = sqrt(err/_ndim);
        // Now use error as an estimate for new step size
        double scale = 0.0;
        const double maxscale = 10.0;
        const double minscale = 0.2;
        const double safe = 0.9;
        const double alpha = 0.2;
        if ( err <= 1.0 ) { // A succesful step...
            if ( err == 0.0 )
                scale = maxscale;
            else {
                // We are NOT using the PI control version as
                // given by Press et al. as I do not want to store
                // the old error from previous integration steps.
                scale = safe * pow(err, -alpha);
                if ( scale < minscale )
                    scale = minscale;
                if ( scale > maxscale )
                    scale = maxscale;
            }
            hSuggest = scale*h;
            return ResultOfStep.success;
        }
        // else, failed step
        scale = max(safe*pow(err, -alpha), minscale);
        hSuggest = scale*h;
        return ResultOfStep.failure;
    }
private:
    int _ndim;
    double _tolerance;
    number[] _yTmp, _yErr;
    number[] _k1, _k2, _k3, _k4, _k5, _k6;
}

/++
 + The Alpha-QSS method, specialised to work
 + with a chemical kinetic ODE.
 +
 + author: K.Damm 2015
 + adapted from R.Gollan's eilmer3 implementation
 +
 + References:
 + Mott, D. (1999)
 + Integrates a system of chemical kinetics ODEs which contains both
 + stiff and non-stiff equations. It is proven to be A-stable and has
 + second order accuracy. The method is single-step self starting,
 + thus has little start-up costs.
 +
 + S.R. Qureshi, R. Prosser (2007)
 + Delta term added into Mott algorithm to prevent excessively small
 + concentrations forcing the solver to take miniscule steps.
 + Default delta = 10.0e-10.
 +/

class AlphaQssStep : ChemODEStep {
public:
    this(in GasModel gmodel, ReactionMechanism rmech, double eps1, double eps2, double delta, int max_iters)
    {
        super(rmech);
        _eps1 = eps1;
        _eps2 = eps2;
        _delta = delta;
        _max_iter = max_iters;
        // Allocate working arrays
        _ndim = gmodel.n_species;
        _yp.length = _ndim;
        _yp0.length = _ndim;
        _yc.length = _ndim;
        _q0.length = _ndim;
        _L0.length = _ndim;
        _p0.length = _ndim;
        _qtilda.length = _ndim;
        _pbar.length = _ndim;
        _qp.length = _ndim;
        _Lp.length = _ndim;
        _pp.length = _ndim;
        _alpha.length = _ndim;
        _alphabar.length = _ndim;
    }

    @nogc
    override ResultOfStep opCall(number[] y0, double h, number[] yOut, ref double hSuggest)
    {
        // 1. Predictor Step
        _rmech.eval_split_rates(y0, _q0, _L0);
        p_on_y(_L0, y0, _p0);
        alpha_compute(_p0, _alpha, h);
        update_conc(_yp, y0, _q0, _p0, _alpha, h);
        // put aside the first predictor values for
        // later use in the convergence test
        foreach ( i; 0 .. _ndim )
          _yp0[i] = _yp[i];

        // 2. Corrector Step
        foreach ( corr; 0 .. _max_iter ) {
            _rmech.eval_split_rates(_yp, _qp, _Lp);
            p_on_y(_Lp, _yp, _pp);
            alpha_compute(_pp, _alphabar, h);
            foreach ( i; 0.._ndim ) {
                _qtilda[i] = _alphabar[i]*_qp[i] + (1-_alphabar[i])*_q0[i];
                _pbar[i] = 0.5*(_p0[i]+_pp[i]);
            }

            // actual corrector step
            update_conc(_yc, y0, _qtilda, _pbar, _alphabar, h);

            bool converged = test_converged(_yc, _yp0, h);

            if ( converged ) {
                // the solver has converged, let's set the current
                // corrector to the yOut vector and suggest a new
                // time step.
                foreach ( i; 0.._ndim ) {
                    yOut[i] = _yc[i];
                }
                hSuggest = step_suggest(h, _yc, _yp0);
                return ResultOfStep.success;
            }
            // if we haven't converged yet make the corrector the new predictor and try again
            foreach ( i; 0 .. _ndim ) {
                _yp[i] = _yc[i];
            }
        }
        // if we have failed the convergence test too many times let's start again with a smaller h
        // This suggested h will be passed through generic step size  check on outer loop
        hSuggest = step_suggest(h, _yc, _yp0);
        return ResultOfStep.failure;
    }

private:
    int _ndim;
    double _eps1;
    double _eps2;
    double _delta;
    int _max_iter;
    immutable double _ZERO_EPS = 1.0e-50;
    number[] _yTmp, _yp, _yp0,  _yc, _q0, _L0, _p0, _pp, _qtilda, _pbar, _qp, _Lp, _alpha, _alphabar;

    // Private functions.
    @nogc
    void p_on_y(number[] L, number[] y, number[] p_y) {
        foreach( i; 0.._ndim)
            p_y[i] = L[i] / (y[i] + _ZERO_EPS);
        return;
    }

    @nogc
    void alpha_compute(number[] p, number[] alpha, double h) {
        foreach ( i; 0.._ndim ) {
            number r = 1.0/(p[i]*h+_ZERO_EPS);  // ZERO_EPS prevents division by 0
            alpha[i] = (180.0*r*r*r+60.0*r*r+11.0*r+1.0)/(360.0*r*r*r+60.0*r*r+12.0*r+1.0);
        }
    }

    @nogc
    void update_conc(number[] yTmp, number[] y0, number[] q, number[] p, number[] alpha, double h) {
        foreach ( i; 0.._ndim )
            yTmp[i] = y0[i] + (h*(q[i]-p[i]*y0[i]))/(1.0+alpha[i]*h*p[i]);
    }

    @nogc
    bool test_converged(in number[] yc, in number[] yp, double h) {
        bool passesTest = true;
        double test = 0.0;
        foreach ( i; 0.._ndim ) {
            if ( yc[i] < _ZERO_EPS )
                continue;
            test = fabs(yc[i].re - yp[i].re);
            // +delta from Qureshi and Prosser (2007)
            if ( test >= (_eps1 * (yc[i].re + _delta)) ) {
                passesTest = false;
                break; // no need to continue testing remaining species
            }
        }
        return passesTest;
    }

    @nogc
    double step_suggest(double h, number[] yc, number[] yp)
    {
        double test = 0.0;
        double sigma = 0.0;
        double h_new = 0.0;

        foreach (i; 0.._ndim) {
            if ( yc[i] < _ZERO_EPS )
                continue;
            test = fabs( yc[i].re - yp[i].re) / (_eps2*(yc[i].re+ _delta));
            if ( test > sigma )
                sigma = test;
        }

        if ( sigma <= 0.0 ) {
            h_new = 10.0*h;
        }
        else {
            double x = sigma;
            x = x - 0.5*(x*x - sigma)/x;
            x = x - 0.5*(x*x - sigma)/x;
            x = x - 0.5*(x*x - sigma)/x;
            h_new = h * ( (1.0 / x) + 0.005 );
        }
    return h_new;
    }
}

version(chemistry_update_test) {
    import std.stdio;
    import util.msg_service;
    import kinetics.rate_constant;
    import kinetics.reaction;
    import gas.therm_perf_gas;

    double numericalEstimate(double dt, double tMax, double[] conc0, ChemODEStep step)
    {
        double t = 0.0;
        double[] conc1;
        conc1.length = conc0.length;
        double dtDummy;
        int count = 0;
        while ( (tMax - t) > 1.0e-9 ) {
            dt = min(dt, tMax - t);
            step(conc0, dt, conc1, dtDummy);
            t += dt;
            conc0 = conc1.dup;
            count++;
        }
        return conc1[2];
    }

    int main() {
        import gas.therm_perf_gas;
        auto gmodel = new ThermallyPerfectGas("sample-input/H2-I2-HI.lua");
        auto rmech = createReactionMechanism("sample-input/H2-I2-inp.lua", gmodel, 100.0, 10000.0);

        GasState gd = GasState(gmodel);
        gd.T = 700.0;
        double c0 = 4.54;
        double[] conc0 = [c0, c0, 0.0];
        gd.p = 2.0*c0*R_universal*gd.T;
        double[] molef = [0.5, 0.5, 0.0];
        gmodel.molef2massf(molef, gd);
        gmodel.update_thermo_from_pT(gd);
        rmech.eval_rate_constants(gmodel, gd);
        double tInterval = 60000.0;
        double analyticalVal = 7.1420197868416215;

        /* 1. Test RKFStep
         * We'll check that the order of accuracy for the method
         * is as expected. We'll choose a timestep and run the problem
         * once and compare the numerical result to an analytical
         * calculation. Then we'll run the problem with a reduced
         * timestep h*pow(2.0, -1/5). For this 5th-order accurate
         * method, we should expect to see the error drop by a factor
         * of 2.0 if we reduce the timestep by 2.0**(-1/5).
         */
        auto rkfStep = new RKFStep(gmodel, rmech, 1.0e-3);
        double dt = 287.175;
        double numVal0 = numericalEstimate(dt, tInterval, conc0, rkfStep);
        double err0 = analyticalVal - numVal0;
        // Reduce timestep and repeat test
        dt *= 2.0^^(-1./5.);
        conc0 = [c0, c0, 0.0];
        double numVal1 = numericalEstimate(dt, tInterval, conc0, rkfStep);
        double err1 = analyticalVal - numVal1;
        assert(isClose(1.96809, err0/err1, 1.0e-4), failedUnitTest());

        /* 2. Test alpha-QSS step
         * We'll check that the order of accuracy for the method
         * is as expected. We'll choose a timestep and run the problem
         * once and compare the numerical result to an analytical
         * calculation. Then we'll run the problem with a reduced
         * timestep h*0.5. For this 2nd-order accurate
         * method, we should expect to see the error drop by a factor
         * of 4.0 if we halve the timestep.
         */
        auto alphaStep = new AlphaQssStep(gmodel, rmech, 0.001, 0.0005, 1.0e-10, 10);
        dt = 2.0;
        numVal0 = numericalEstimate(dt, tInterval, conc0, alphaStep);
        err0 = analyticalVal - numVal0;
        // Reduce timestep and repeat test
        dt *= 0.5;
        conc0 = [c0, c0, 0.0];
        numVal1 = numericalEstimate(dt, tInterval, conc0, alphaStep);
        err1 = analyticalVal - numVal1;
        assert(isClose(7.1420197868416215, numVal1, 1.0e-4), failedUnitTest());
        assert(isClose(4.001, err0/err1, 1.0e-4), failedUnitTest());

        /* 3. Test the complete update algorithm as used
         *    by the flow solver. This might be stretching a
         *    bit what a *unit* is. In this test, we'll exercise
         *    the chemistry update routine by solving the
         *    complete hydrogen-iodine system.
         */

        double dtSuggest = 200.0;
        auto chemUpdate = new ChemistryUpdate("sample-input/H2-I2-inp.lua", gmodel);
        double[maxParams] params;
        chemUpdate(gd, tInterval, dtSuggest, params);
        double[] conc;
        conc.length = 3;
        gmodel.massf2conc(gd, conc);
        assert(isClose(7.14201983840, conc[2], 1.0e-9), failedUnitTest());

        return 0;
    }
}
