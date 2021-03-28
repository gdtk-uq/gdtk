/**
 * Author: Rowan G.
 * Date: 2021-03-28
 */

module kinetics.two_temperature_thermochemical_reactor;

import std.math;
import std.conv;
import nm.complex;
import nm.number;

import util.lua;
import util.lua_service;
import gas;
import kinetics.thermochemical_reactor;

class TwoTemperatureThermochemicalReactor : ThermoChemicalReactor {
public:
    @nogc
    override void opCall(GasState gs, double tInterval,
                         ref double dtChemSuggest, ref double dtThermSuggest,
                         ref number[maxParams] params)
    {
        mGsInit.copy_values_from(gs);
        mMgmodel.massf2conc(gs, mConc0);
        m_uTotal = mGmodel.internal_energy(gs);
        m_uVib0 = gs.u_modes[0];

        // Sort out stepsize for possible subcycling.
        number t = 0.0;
        double h;
        double dtSave;
        double dtSuggest = dtChemSuggest;
        if (dtChemSuggest > tInterval) {
            h = tInterval;
        }
        else if (dtChemSuggest <= 0.0) {
            h = mRmech.estimateStepSize(mConc0);
        }
        else {
            h = dtChemSuggest;
        }
        
        // Now the interesting stuff: increment change to species composition and energy
        int cycle = 0;
        int attempt = 0;
        // Populate y0
        // Species concentrations in first part of array
        foreach (isp; 0 .. mNSpecies) m_y0[isp] = mConc0[isp];
        // Final entry is energy change.
        m_y0[mNSpecies] = gs.u_modes[0];
        foreach ( ; cycle < mMaxSubcycles; ++cycle) {
            /* Copy the last good timestep suggestion before moving on.
             * If we're near the end of the interval, we want
             * the penultimate step. In other words, we don't
             * want to store the fractional step of (tInterval - t)
             * that is taken as the last step.
             */
            dtSave = dtChemSuggest;
            h = min(h, tInterval - t);
            attempt = 0;
            foreach ( ; attempt < mMaxAttempts; ++attempt) {
                ResultOfStep result = step(m_y0, h, m_yOut, dtChemSuggest);
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
                     * the dtChemSuggest comes back lower, let's just set
                     * h as the original value for the successful step.
                     */
                    double hMax = h*(1.0 + DT_INCREASE_PERCENT/100.0);
                    if (dtChemSuggest > h) {
                        /* It's possible that dtChemSuggest is less than h.
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
                        h = min(dtChemSuggest, hMax);
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
                string errMsg = "Hit maxmium number of step attempts within a subcycle in 2-T thermochemical update.";
                // We did poorly. Return GasState to how we found it, and throw an Exception.
                gs.copy_values_from(mGsInit);
                throw new ThermochemicalReactorUpdateException(errMsg);
            }
            // We employ a tight temperature coupling approach here.
            // Since 2-T models are typically involved with high-temperature gases,
            // the temperature changes can be large over total tInterval.
            // Given this is the case, it will help robustness to re-evaluate
            // the temperature-dependent rate constants and relaxation times.
            gs.u_modes[0] = m_y0[mNSpecies];
            gs.u = m_uTotal - gs.u_modes[0];
            mGmodel.update_thermo_from_rhou(gs);
            mRmech.eval_rate_constants(gs);

            if (t >= tInterval) { // We've done enough cycling.
                // If we've only taken one cycle, then we would like to use
                // the this new suggested dtChemSuggest rather than the one from
                // the penultimate step.
                if (cycle == 0) {
                    dtSave = dtChemSuggest;
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

        dtChemSuggest = dtSave;
    }

private:
    int mMaxSubcycles;
    int mMaxAttempts;
    number[] mConc0;
    number[] m_y0;
    number[] m_yOut;
    number[][] mRk45WS;
    number m_uTotal;
    GasModel mGmodel;
    ReactionMechanism mRmech;
    EnergyExchangeSystem mEES;
    
}
