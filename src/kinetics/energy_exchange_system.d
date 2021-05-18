/**
 * Interface and implementation of energy exchange system.
 *
 * Author: Rowan G.
 * Date: 2021-03-28
 */

module kinetics.energy_exchange_system;

import std.format;
import std.string;
import std.stdio;
import std.math;
import std.conv;

import nm.complex;
import nm.number;

import util.lua;
import util.lua_service;
import gas;

import kinetics.energy_exchange_mechanism;
import kinetics.relaxation_time;
import kinetics.reaction_mechanism;

interface EnergyExchangeSystem {
    @nogc void evalRelaxationTimes(in GasState gs);
    @nogc void evalRates(in GasState gs, in ReactionMechanism, ref number[] rates);
}

class TwoTemperatureEnergyExchange : EnergyExchangeSystem {
public:
    this(string fname, GasModel gmodel)
    {
        // For 2-T model, one entry in T_modes, so index is 0.
        int mode = 0;
        mGmodel = gmodel;
        mGsEq = new GasState(gmodel);
        mMolef.length = gmodel.n_species;
        mNumden.length = gmodel.n_species;

        // Load in table of energy exchange mechanisms from the verbose lua file 
        auto L = init_lua_State();
        doLuaFile(L, fname);

        lua_getglobal(L, "mechanism");
        lua_pushnil(L); // dummy first key
        while (lua_next(L, -2) != 0) { // -1 is the dummy key, -2 is the mechanism table
            mEEM ~= createEnergyExchangeMechanism(L, mode, gmodel);
            lua_pop(L, 1); // discard value but keep key so that lua_next can remove it (?!)
        }
        lua_pop(L, 1); // remove mechanisms table
        lua_close(L);
    }

    @nogc
    void evalRelaxationTimes(in GasState gs)
    {
        mGmodel.massf2molef(gs, mMolef);
        mGmodel.massf2numden(gs, mNumden);
        foreach (mech; mEEM) {
            mech.evalRelaxationTime(gs, mMolef, mNumden);
        }
    }

    @nogc
    void evalRates(in GasState gs, in ReactionMechanism rmech, ref number[] rates)
    {
        // Compute a star state at transrotational equilibrium.
        mGsEq.copy_values_from(gs);
        mGsEq.T_modes[0] = gs.T;
        mGmodel.update_thermo_from_pT(mGsEq);
        number rate = 0.0;
        // Compute rate change from VT exchange
        foreach (mech; mEEM) {
            rate += mech.rate(gs, mGsEq, mMolef, mNumden, rmech);
        }
        rates[0] = rate;
    }

private:
    int[] mVibRelaxers;
    int[][] mHeavyParticles;
    number[] mMolef;
    number[] mNumden;
    GasModel mGmodel;
    GasState mGsEq;
    EnergyExchangeMechanism[] mEEM;
    //EnergyExchangeMechanism[][] mVT;
    // EnergyExchangeMechanism[] mET;
    // EnergyExchangeMechanism[] mChemCoupling;
}

