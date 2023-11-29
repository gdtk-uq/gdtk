/**
 * yee_kotov_kinetics.d
 *
 * Binary reacting gas as appears in:
 * Yee, Kotov, Wang and Shu (2013)
 * Spurious behavior of shock-capturing methods by the fractional step
 * approach: Problems containing stiff source terms and discontinuities.
 * Journal of Computational Physics, 241:pp. 266--291
 *
 * Kotov, Yee, Panesi, Prabhu and Wray (2014)
 * Computational challenges for simulations related to the NASA
 * electric arc shock tube (EAST) experiments.
 * Journal of Computational Physics, 269:pp. 215--233
 *
 * This kinetics file accompanies the gas model in gas/ideal_gas_ab.d
 *
 * Authors: Rowan G., Jamie B. and Peter J.
 * Version: 2020-05-10, initial cut derived from powers_aslam_kinetics.d
 */

module kinetics.yee_kotov_kinetics;

import std.stdio;
import std.math;
import std.format;
import ntypes.complex;
import nm.number;

import gas;
import gas.ideal_gas_ab;
import util.lua;
import util.lua_service;
import kinetics.thermochemical_reactor;

immutable double DT_INCREASE_PERCENT = 10.0; // allowable percentage increase on succesful step
immutable double DT_DECREASE_PERCENT = 50.0; // allowable percentage decrease on succesful step
                                             // Yes, you read that right. Sometimes a step is succesful
                                             // but the timestep selection algorithm will suggest
                                             // a reduction. We limit that reduction to no more than 50%.
immutable double DT_REDUCTION_FACTOR = 10.0; // factor by which to reduce timestep
                                             // after a failed attempt
immutable double H_MIN = 1.0e-15; // Minimum allowable step size
enum ResultOfStep { success, failure };

final class UpdateAB_YeeKotov : ThermochemicalReactor {

    this(string fname, GasModel gmodel)
    {
        super(gmodel); // hang on to a reference to the gas model
        _igABmodel = cast(IdealGasAB) gmodel;
        if (_igABmodel is null) {
            string errMsg = "Error in construction of UpdateAB_YeeKotov.\n";
            errMsg ~= "The supplied gas model must be an IdealGasAB model.\n";
            throw new ThermochemicalReactorUpdateException(errMsg);
        }
        _Qinit = GasState(gmodel);
        _Q0 = GasState(gmodel);
        // We need to pick a number of pieces out of the gas-model file, again.
        // Although they exist in the GasModel object, they are private.
        auto L = init_lua_State();
        doLuaFile(L, fname);
        lua_getglobal(L, "IdealGasAB");
        // Now, pull out the numeric value parameters.
        _K_0 = getDouble(L, -1, "K_0");
        _T_ign = getDouble(L, -1, "T_ign");
        lua_pop(L, 1); // dispose of the table
        lua_close(L);
    }

    @nogc
    override void opCall(ref GasState Q, double tInterval, ref double dtSuggest,
                         ref number[maxParams] params)
    {
        try {
            massfUpdate(Q, tInterval, dtSuggest);
        }
        catch (GasModelException err) {
            string msg = "The mass fraction update in the Yee-Kotov kinetics module failed.\n";
            debug { msg ~= format("\ncaught %s", err.msg); }
            throw new ThermochemicalReactorUpdateException(msg);
        }
        _gmodel.update_thermo_from_rhou(Q);
        _gmodel.update_sound_speed(Q);
    }

    @nogc override void eval_source_terms(GasModel gmodel, ref GasState Q, ref number[] source) {
        string errMsg = "eval_source_terms not implemented for yee_kotov_kinetics.";
        throw new ThermochemicalReactorUpdateException(errMsg);
    }

private:
    double _K_0;   // reaction rate constant
    double _T_ign; // ignition temperature

    IdealGasAB _igABmodel;
    GasState _Qinit;
    GasState _Q0;
    double _massfAbsTolerance = 1.0e-9;
    double _massfRelTolerance = 1.0e-9;
    int _maxSubcycles = 10000;
    int _maxAttempts = 3;

    @nogc
    void massfUpdate(ref GasState Q, double tInterval, ref double dtSuggest)
    {
        // We borrow the algorithm from ChemistryUpdate.opCall()
        // Take a copy of what's passed in, in case something goes wrong
        _Qinit.copy_values_from(Q);

        // 1. Sort out the timestep for possible subcyling
        double t = 0.0;
        double h;
        double dtSave;
        if ( dtSuggest > tInterval )
            h = tInterval;
        else if ( dtSuggest <= 0.0 )
            h = 0.01*tInterval;
        else
            h = dtSuggest;

        // 2. Now do the interesting stuff, increment the mass fractions
        int cycle = 0;
        int attempt = 0;
        for ( ; cycle < _maxSubcycles; ++cycle) {
            /* Copy the last good timestep suggestion before moving on.
             * If we're near the end of the interval, we want
             * the penultimate step. In other words, we don't
             * want to store the fractional step of (tInterval - t)
             * that is taken as the last step.
             */
            dtSave = dtSuggest;
            h = fmin(h, tInterval - t);
            attempt = 0;
            for ( ; attempt < _maxAttempts; ++attempt) {
                ResultOfStep result = updateStep(Q, h, dtSuggest);
                if (result == ResultOfStep.success) {
                    t += h;
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
                        h = fmin(dtSuggest, hMax);
                    }
                    break;
                }
                else { // in the case of failure...
                    /* We now need to make some decision about
                     * what timestep to attempt next.
                     */
                    h /= DT_REDUCTION_FACTOR;
                    if (h < H_MIN) {
                        string errMsg = "Hit the minimum allowable timestep in Yee-Kotov mass fraction update.";
                        debug { errMsg ~= format("\ndt= %.4e", H_MIN); }
                        Q.copy_values_from(_Qinit);
                        throw new ThermochemicalReactorUpdateException(errMsg);
                    }
                }
            } // end attempts at single subcycle.
            if (attempt == _maxAttempts) {
                string errMsg = "Hit maximum number of step attempts within a subcycle for Yee-Kotov mass fraction update.";
                // We did poorly. Let's put the original GasState back in place,
                // and let the outside world know via an Exception.
                Q.copy_values_from(_Qinit);
                throw new ThermochemicalReactorUpdateException(errMsg);
            }
            /* Otherwise, we've done well. */
            if (t >= tInterval) {
                // We've done enough work.
                // If we've only take one cycle, then we would like to use
                // dtSuggest rather than using the dtSuggest from the
                // penultimate step.
                if (cycle == 0) {
                    dtSave = dtSuggest;
                }
                break;
            }
        }
        if (cycle == _maxSubcycles) {
            // If we make it out here, then we didn't complete within
            // the maximum number of subscyles... we are taking too long.
            // Let's return the gas state to how it was before we failed
            // and throw an Exception.
            string errMsg = "Hit maximum number of subcycles while attempting Yee-Kotov mass fraction update.";
            Q.copy_values_from(_Qinit);
            throw new ThermochemicalReactorUpdateException(errMsg);
        }
        // At this point, it appears that everything has gone well.
        dtSuggest = dtSave;
    }

    @nogc
    number evalRate(ref GasState Q)
    {
        return -_K_0 * exp(-_T_ign/Q.T) * Q.massf[0];
    }

    @nogc
    ResultOfStep updateStep(ref GasState Q, double h, ref double hSuggest)
    {
        _Q0.copy_values_from(Q);
        immutable double a21=1./5., a31=3./40., a32=9./40.,
            a41=3./10., a42=-9./10., a43 = 6./5.,
            a51=-11./54., a52=5./2., a53=-70./27., a54=35./27.,
            a61=1631./55296., a62=175./512., a63=575./13824., a64=44275./110592., a65=253./4096.;
        immutable double b51=37./378., b53=250./621., b54=125./594., b56=512./1771.,
            b41=2825./27648., b43=18575./48384., b44=13525./55296., b45=277./14336., b46=1.0/4.0;

        number k1 = evalRate(Q);
        Q.massf[0] = _Q0.massf[0] + h*(a21*k1);
        Q.massf[1] = 1.0 - Q.massf[0];
        _igABmodel.update_thermo_from_rhou(Q);

        number k2 = evalRate(Q);
        Q.massf[0] = _Q0.massf[0] + h*(a31*k1 + a32*k2);
        Q.massf[1] = 1.0 - Q.massf[0];
        _igABmodel.update_thermo_from_rhou(Q);

        number k3 = evalRate(Q);
        Q.massf[0] = _Q0.massf[0] + h*(a41*k1 + a42*k2 + a43*k3);
        Q.massf[1] = 1.0 - Q.massf[0];
        _igABmodel.update_thermo_from_rhou(Q);

        number k4 = evalRate(Q);
        Q.massf[0] = _Q0.massf[0] + h*(a51*k1 + a52*k2 + a53*k3 + a54*k4);
        Q.massf[1] = 1.0 - Q.massf[0];
        _igABmodel.update_thermo_from_rhou(Q);

        number k5 = evalRate(Q);
        Q.massf[0] = _Q0.massf[0] + h*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5);
        Q.massf[1] = 1.0 - Q.massf[0];
        _igABmodel.update_thermo_from_rhou(Q);

        number k6 = evalRate(Q);
        Q.massf[0] = _Q0.massf[0] + h*(b51*k1 + b53*k3 + b54*k4 + b56*k6);
        Q.massf[1] = 1.0 - Q.massf[0];
        _igABmodel.update_thermo_from_rhou(Q);

        // Compute error estimate.
        number errEst = Q.massf[0] - (_Q0.massf[0] + h*(b41*k1 + b43*k3 + b44*k4 + b45*k5 + b46*k6));

        // And use error estimate as a means to suggest a new timestep.
        double atol = _massfAbsTolerance;
        double rtol = _massfRelTolerance;
        double sk = atol + rtol*fmax(fabs(_Q0.massf[0].re), fabs(Q.massf[0].re));
        double err = fabs(errEst.re/sk);
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
        scale = fmax(safe*pow(err, -alpha), minscale);
        hSuggest = scale*h;
        return ResultOfStep.failure;
    }

} // end class UpdateAB_YeeKotov

