/**
 * Author: Rowan G.
 * Date: 2021-03-28
 */

module kinetics.multi_temperature_thermochemical_reactor;

import std.format;
import std.math;
import std.algorithm;
import std.conv;
import std.stdio;
import ntypes.complex;
import nm.number;

import util.lua;
import util.lua_service;
import gas;
import gas.composite_gas;
import kinetics.thermochemical_reactor;
import kinetics.reaction_mechanism;
import kinetics.energy_exchange_system;

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

class MultiTemperatureThermochemicalReactor : ThermochemicalReactor {
public:

    this(string fname1, string fname2, GasModel gmodel)
    {
        super(gmodel);
        mGmodel = gmodel;
        mNSpecies = mGmodel.n_species;
        mNModes = mGmodel.n_modes;
        mGsInit = GasState(gmodel);

        auto L = init_lua_State();
        doLuaFile(L, fname1);
        lua_getglobal(L, "reaction");
        mRmech = createReactionMechanism(L, gmodel, 300.0, 30000.0);

        // Initialise energy exchange mechanism
        if (typeid(gmodel) is typeid(CompositeGas)) {
            auto cg = cast(CompositeGas) gmodel;
            switch (cg.physicalModel) {
                case "two-temperature-gas":
                    mEES = new TwoTemperatureEnergyExchange(fname2, gmodel);
                    break;
                case "three-temperature-gas":
                    mEES = new ThreeTemperatureEnergyExchange(fname2, gmodel);
                    break;
                case "multi-temperature-gas":
                    mEES = new MultiTemperatureEnergyExchange(fname2, gmodel);
                    break;
                default:
                    throw new Exception("Unkown physical model");
            }
        }
        else {
            mEES = new TwoTemperatureEnergyExchange(fname2, gmodel);
        }

        // Set up the rest of the parameters.
        lua_getglobal(L, "config");
        lua_getfield(L, -1, "odeStep");
        lua_getfield(L, -1, "method");
        string ode_method = to!string(lua_tostring(L, -1));
        lua_pop(L, 1);
        switch (ode_method) {
        case "rkf":
             lua_getfield(L, -1, "errTol");
             if (!lua_isnil(L, -1)) mTol = lua_tonumber(L, -1);
             lua_pop(L, 1);
             break;
        default:
             string errMsg = format("ERROR: The odeStep '%s' cannot be used with TwoTemperatureThermochemicalReactor.\n", ode_method);
             throw new Error(errMsg);
        }
        lua_pop(L, 1); // pops 'odeStep'

        lua_getfield(L, -1, "maxSubcycles");
        if (!lua_isnil(L, -1)) mMaxSubcycles = to!int(luaL_checkinteger(L, -1));
        lua_pop(L, 1);

        lua_getfield(L, -1, "maxAttempts");
        if (!lua_isnil(L, -1)) mMaxAttempts = to!int(luaL_checkinteger(L, -1));
        lua_pop(L, 1);

        lua_pop(L, 1); // pops 'config'
        lua_close(L);

        initialiseWorkSpace();
    }

    void initialiseWorkSpace()
    {
        therm_rates.length = mNModes;
        mConc0.length = mNSpecies;
        mConc.length = mNSpecies;
        m_dCdt.length = mNSpecies;
        m_duvedt.length = mNModes;
        m_y0.length =  mNSpecies + mNModes;
        m_yOut.length = mNSpecies + mNModes;
        m_yTmp.length = mNSpecies + mNModes;
        m_yErr.length = mNSpecies + mNModes;
        m_k1.length = mNSpecies + mNModes;
        m_k2.length = mNSpecies + mNModes;
        m_k3.length = mNSpecies + mNModes;
        m_k4.length = mNSpecies + mNModes;
        m_k5.length = mNSpecies + mNModes;
        m_k6.length = mNSpecies + mNModes;
    }

    @nogc
    override void opCall(ref GasState gs, double tInterval, ref double dtSuggest,
                         ref number[maxParams] params)
    {
        mGsInit.copy_values_from(gs);
        mGmodel.massf2conc(gs, mConc0);
        m_uTotal = mGmodel.internal_energy(gs);

        // Evaluate temperature dependent rate parameters
        mRmech.eval_rate_constants(mGmodel, gs);
        mEES.evalRelaxationTimes(gs);


        // Sort out stepsize for possible subcycling.
        double t = 0.0;
        double h;
        double dtSave;
        if (dtSuggest > tInterval) {
            h = tInterval;
        }
        else if (dtSuggest <= 0.0) {
            h = mRmech.estimateStepSize(mConc0);
        }
        else {
            h = dtSuggest;
        }

        // Now the interesting stuff: increment change to species composition and energy
        int cycle = 0;
        int attempt = 0;
        // Populate y0
        // Species concentrations in first part of array
        foreach (isp; 0 .. mNSpecies) m_y0[isp] = mConc0[isp];
        //Final entry is energy change.
        m_y0[mNSpecies .. $] = gs.u_modes[];
        for ( ; cycle < mMaxSubcycles; ++cycle) {
            mRmech.eval_rate_constants(mGmodel, gs);
            mEES.evalRelaxationTimes(gs);
            /* Copy the last good timestep suggestion before moving on.
             * If we're near the end of the interval, we want
             * the penultimate step. In other words, we don't
             * want to store the fractional step of (tInterval - t)
             * that is taken as the last step.
             */
            dtSave = dtSuggest;
            h = min(h, tInterval - t);
            attempt = 0;
            for ( ; attempt < mMaxAttempts; ++attempt) {
                ResultOfStep result = step(gs, m_y0, h, m_yOut, dtSuggest);
                // Unpack m_yOut
                foreach (isp; 0 .. mNSpecies) mConc0[isp] = m_yOut[isp];

                // Now check step is alright.
                bool passedMassFractionTest = true;
                if (result == ResultOfStep.success) {
                    /* We successfully took a step of size h according to the ODE method.
                     * However, we need to test that the mass fractions have remained ok.
                     */
                    mGmodel.conc2massf(mConc0, gs);
                    auto massfTotal = sum(gs.massf);
                    if (fabs(massfTotal - 1.0) > ALLOWABLE_MASSF_ERROR) {
                        passedMassFractionTest = false;
                    }
                }

                if (result == ResultOfStep.success && passedMassFractionTest) {
                    t += h;
                    // Copy successful values into m_y0 for next step.
                    m_y0[] = m_yOut[];
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
                else { // in the case of failure
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
                        string errMsg = "Hit the minimum allowable timestep in 2-T thermochemical update.";
                        debug { errMsg ~= format("\ndt= %.4e", H_MIN); }
                        gs.copy_values_from(mGsInit);
                        throw new ThermochemicalReactorUpdateException(errMsg);
                    }
                }
            } // end attempts at single step of subcycle.
            if (attempt == mMaxAttempts) {
                string errMsg = "Hit maxmium number of step attempts within a subcycle in multi-T thermochemical update.";
                // We did poorly. Return GasState to how we found it, and throw an Exception.
                gs.copy_values_from(mGsInit);
                throw new ThermochemicalReactorUpdateException(errMsg);
            }
            // We employ a tight temperature coupling approach here.
            // Since 2-T models are typically involved with high-temperature gases,
            // the temperature changes can be large over total tInterval.
            // Given this is the case, it will help robustness to re-evaluate
            // the temperature-dependent rate constants and relaxation times.
            gs.u_modes[] = m_y0[mNSpecies .. $];
            gs.u = m_uTotal - sum(gs.u_modes);
            try {
                mGmodel.update_thermo_from_rhou(gs);
            } catch (GasModelException err) {
                gs.copy_values_from(mGsInit);
                string msg = "Failed to update thermodynamic state in multi-T thermochemical update. ";
                debug {
                    msg ~= err.toString();
                    msg ~= " The initial gas state was\n";
                    msg ~= mGsInit.toString();
                }
                throw new ThermochemicalReactorUpdateException(msg);
            }

            if (t >= tInterval) { // We've done enough cycling.
                // If we've only taken one cycle, then we would like to use
                // the this new suggested dtSuggest rather than the one from
                // the penultimate step.
                if (cycle == 0) {
                    dtSave = dtSuggest;
                }
                break;
            }
        } // end cycles loop
        if (cycle == mMaxSubcycles) {
            // If we make it out here, then we didn't complete within the maximum number of subcycles.
            // That is, we are taking too long.
            // Let's return the GasState to how it was before we failed, and throw an Exception.
            string errMsg = "Hit maximum number of subcycles while attempting 2-T thermochemical update.";
            gs.copy_values_from(mGsInit);
            throw new ThermochemicalReactorUpdateException(errMsg);
        }
        // At this point, it appears that everything has gone well.
        // We'll tweak the mass fractions in case they are a little off.
        auto massfTotal = sum(gs.massf);
        foreach (ref mf; gs.massf) mf /= massfTotal;

        dtSuggest = dtSave;
    }

    @nogc override void eval_source_terms(GasModel gmodel, ref GasState Q, ref number[] source) {
        // species source terms
        auto chem_source = source[0..mNSpecies]; // these are actually references not copies
        mRmech.eval_source_terms(gmodel, Q, chem_source);
        // energy modes source terms
        auto therm_source = source[mNSpecies..source.length];
        mEES.eval_source_terms(mRmech, Q, therm_rates, therm_source);
    }

private:
    number[] therm_rates;
    int mNSpecies;
    int mNModes;
    int mMaxSubcycles=10000;
    int mMaxAttempts=4;
    double mTol=1e-3;
    number[] mConc0;
    number[] mConc;
    number[] m_dCdt;
    number[] m_duvedt;
    number[] m_y0;
    number[] m_yOut;
    number m_uTotal;
    GasModel mGmodel;
    GasState mGsInit;
    ReactionMechanism mRmech;
    EnergyExchangeSystem mEES;
    // Working arrays for RKF step method.
    number[] m_yTmp;
    number[] m_yErr;
    number[] m_k1;
    number[] m_k2;
    number[] m_k3;
    number[] m_k4;
    number[] m_k5;
    number[] m_k6;

    @nogc
    void evalRates(ref GasState gs, number[] y, ref number[] rates) {
        mConc[] = y[0 .. mNSpecies];
        mRmech.eval_rates(mConc, m_dCdt);
        rates[0 .. mNSpecies] = m_dCdt[];

        gs.u_modes[] = y[mNSpecies .. $];
        gs.u = m_uTotal - sum(gs.u_modes);
        mEES.evalRates(gs, mRmech, m_duvedt);
        rates[mNSpecies .. $] = m_duvedt[];
    }

    @nogc
    ResultOfStep step(ref GasState gs, number[] y0, double h, ref number[] yOut, ref double hSuggest)
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
        evalRates(gs, y0, m_k1);
        foreach (i; 0 .. y0.length) m_yTmp[i] = y0[i] + h*a21*m_k1[i];

        evalRates(gs, m_yTmp, m_k2);
        foreach (i; 0 .. y0.length) m_yTmp[i] = y0[i] + h*(a31*m_k1[i] + a32*m_k2[i]);

        evalRates(gs, m_yTmp, m_k3);
        foreach (i; 0 .. y0.length) m_yTmp[i] = y0[i] + h*(a41*m_k1[i] + a42*m_k2[i] + a43*m_k3[i]);

        evalRates(gs, m_yTmp, m_k4);
        foreach (i; 0 .. y0.length) m_yTmp[i] = y0[i] + h*(a51*m_k1[i] + a52*m_k2[i] + a53*m_k3[i] + a54*m_k4[i]);

        evalRates(gs, m_yTmp, m_k5);
        foreach (i; 0 .. y0.length) m_yTmp[i] = y0[i] + h*(a61*m_k1[i] + a62*m_k2[i] + a63*m_k3[i] + a64*m_k4[i] + a65*m_k5[i]);

        evalRates(gs, m_yTmp, m_k6);

        // 2. Compute new value and error esimate
        foreach (i; 0 .. y0.length) yOut[i] = y0[i] + h*(b51*m_k1[i] + b53*m_k3[i] + b54*m_k4[i] + b56*m_k6[i]);
        foreach (i; 0 .. y0.length) m_yErr[i] = yOut[i] - (y0[i] + h*(b41*m_k1[i] + b43*m_k3[i] + b44*m_k4[i] + b45*m_k5[i] + b46*m_k6[i]));

        // 3. Lastly, use error estimate as a means to suggest
        //    a new timestep.
        //    Compute error using tol as atol and rtol as suggested in
        //    Press et al. (2007)
        double err = 0.0;
        double sk = 0.0;
        double atol = mTol;
        double rtol = mTol;

        foreach (i; 0 .. m_yErr.length) {
            sk = atol + rtol*max(fabs(y0[i].re), fabs(yOut[i].re));
            err += (m_yErr[i].re/sk)*(m_yErr[i].re/sk);
        }
        err = sqrt(err/m_yErr.length);
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

}
