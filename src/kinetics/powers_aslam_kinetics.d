/**
 * powers_aslam_kinetics.d
 *
 * Two-component reacting gas as described in.
 * JM Powers and TD Aslam (2006)
 * Exact solution for multidimensional compressible reactive flow
 * for verifying numerical algorithms.
 * AIAA Journal Vol. 44 No. 2 pages 337-344
 *
 * This kinetics file accompanies the gas model in gas/ideal_gas_ab.d
 *
 * Authors: Peter J. and Rowan G.
 * Version: 2017-01-07: Initial cut.
 *          2021-06-26: Add eval_source_terms() code.
 */

module kinetics.powers_aslam_kinetics;

import std.math;
import std.conv;
import ntypes.complex;
import nm.number;

import gas;
import util.lua;
import util.lua_service;
import kinetics.thermochemical_reactor;

final class UpdateAB : ThermochemicalReactor {

    this(string fname, GasModel gmodel)
    {
        super(gmodel); // hang on to a reference to the gas model
        // We need to pick a number of pieces out of the gas-model file, again.
        // Although they exist in the GasModel object, they are private.
        auto L = init_lua_State();
        doLuaFile(L, fname);
        lua_getglobal(L, "IdealGasAB");
        // Now, pull out the numeric value parameters.
        _alpha = getDouble(L, -1, "alpha");
        _Ti = getDouble(L, -1, "Ti");
        lua_pop(L, 1); // dispose of the table
        lua_close(L);
    }

    override void opCall(ref GasState Q, double tInterval, ref double dtSuggest,
                         ref number[maxParams] params)
    {
        if (Q.T > _Ti) {
            // We are above the ignition point, proceed with reaction.
            number massfA = Q.massf[0];
            number massfB = Q.massf[1];
            // This gas has a very simple reaction scheme that can be integrated explicitly.
            massfA = massfA*exp(-_alpha*tInterval);
            massfB = 1.0 - massfA;
            Q.massf[0] = massfA;
            Q.massf[1] = massfB;
        } else {
            // do nothing, since we are below the ignition temperature
        }
        // Since the internal energy and density in the (isolated) reactor is fixed,
        // we need to evaluate the new temperature, pressure, etc.
        _gmodel.update_thermo_from_rhou(Q);
        _gmodel.update_sound_speed(Q);
    } // end opCall()

    @nogc override void eval_source_terms(GasModel gmodel, ref GasState Q, ref number[] source) {
        if (Q.T > _Ti) {
            // We are above the ignition point, proceed with reaction.
            number rhoA = Q.massf[0] * Q.rho;
            number drhoB_dt = _alpha*rhoA;
            source[0] = -drhoB_dt;
            source[1] = drhoB_dt;
        } else {
            // Zero reaction rates, since we are below the ignition temperature.
            source[0] = to!number(0.0);
            source[1] = to!number(0.0);
        }
    } // end eval_source_terms()

private:
    // Reaction rate constant
    double _alpha; // 1/s
    // Ignition temperature
    double _Ti; // degrees K
} // end class UpdateAB
