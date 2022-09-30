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

class EnergyExchangeSystem {
public:
    this (GasModel gmodel)
    {
        mGmodel = gmodel;
        mGsEq = GasState(gmodel);
        mMolef.length = gmodel.n_species;
        mNumden.length = gmodel.n_species;
    }

    @nogc void evalRelaxationTimes(in GasState gs)
    {
        mGmodel.massf2molef(gs, mMolef);
        mGmodel.massf2numden(gs, mNumden);
        foreach (mech; mEEM) {
            mech.evalRelaxationTime(gs, mMolef, mNumden);
        }
    }

    @nogc void evalRates(in GasState gs, in ReactionMechanism rmech, ref number[] rates)
    {
        // Compute a star state at transrotational equilibrium.
        mGsEq.copy_values_from(gs);
        mGsEq.T_modes[] = gs.T;
        mGmodel.update_thermo_from_pT(mGsEq);
        rates[] = to!number(0.0);
        // Compute rate change from VT exchange
        foreach (mech; mEEM) {
            rates[mech.mode_i] += mech.rate(gs, mGsEq, mMolef, mNumden, rmech);
        }
    }

    @nogc
    void eval_source_terms(in ReactionMechanism rmech, ref GasState Q, ref number[] rates, ref number[] source)
    {
        evalRelaxationTimes(Q);
        evalRates(Q, rmech, rates);
        rates2source(Q, rates, source);
    }

    @nogc
    void rates2source(in GasState Q, in number[] rates, ref number[] source)
    {
        foreach(i, rate; rates) source[i] = rate*Q.rho;
    }

private:
    GasModel mGmodel;
    GasState mGsEq;
    number[] mMolef;
    number[] mNumden;
    EnergyExchangeMechanism[] mEEM;
}

class TwoTemperatureEnergyExchange : EnergyExchangeSystem {
public:
    this(string fname, GasModel gmodel)
    {
        super(gmodel);

        // For 2-T model, one entry in T_modes, so index is 0.
        int mode_i = 0;

        // For 2-T, the other mode is translation/rotation. We represent this with -1
        int mode_j = -1;

        // Load in table of energy exchange mechanisms from the verbose lua file
        auto L = init_lua_State();
        doLuaFile(L, fname);

        lua_getglobal(L, "mechanism");
        lua_pushnil(L); // dummy first key
        while (lua_next(L, -2) != 0) { // -1 is the dummy key, -2 is the mechanism table
            mEEM ~= createEnergyExchangeMechanism(L, mode_i, mode_j, gmodel);
            lua_pop(L, 1); // discard value but keep key so that lua_next can remove it (?!)
        }
        lua_pop(L, 1); // remove mechanisms table
        lua_close(L);
    }
}

class ThreeTemperatureEnergyExchange : EnergyExchangeSystem {
public:
    this(string fname, GasModel gmodel)
    {
        super(gmodel);
        mGsEqE = GasState(gmodel);

        // Load in table of energy exchange mechanisms from the verbose lua file
        auto L = init_lua_State();
        doLuaFile(L, fname);

        lua_getglobal(L, "mechanism");
        lua_pushnil(L); // dummy first key
        while (lua_next(L, -2) != 0) { // -1 is the dummy key, -2 is the mechanism table
            int mode_i, mode_j;
            string type = getString(L, -1, "type");
            switch (type) {
                case "V-T":
                    mode_i = 0;
                    mode_j = -1;
                    break;
                case "E-T":
                    mode_i = 1;
                    mode_j = -1;
                    break;
                case "C-E":
                    mode_i = 1;
                    mode_j = -1;
                    break;
                case "V-E":
                    mode_i = 0;
                    mode_j = 1;
                    break;
                case "E-V":
                    mode_i = 0;
                    mode_j = 1;
                    break;
                case "C-V":
                    mode_i = 0;
                    mode_j = -1;
                    break;
                default:
                    throw new Error("Unknown type of energy exchange");
            }
            mEEM ~= createEnergyExchangeMechanism(L, mode_i, mode_j, gmodel);
            lua_pop(L, 1); // discard value but keep key so that lua_next can remove it (?!)
        }
        lua_pop(L, 1); // remove mechanisms table
        lua_close(L);
    }

    @nogc override void evalRates(in GasState gs, in ReactionMechanism rmech, ref number[] rates)
    {
        // Compute a star state at transrotational equilibrium.
        mGsEq.copy_values_from(gs);
        mGsEq.T_modes[] = gs.T;
        mGmodel.update_thermo_from_pT(mGsEq);

        // Compute a star state at vibration/electron equilibrium.
        mGsEqE.copy_values_from(gs);
        mGsEqE.T_modes[] = gs.T_modes[1];
        mGsEqE.T = gs.T_modes[1];
        mGmodel.update_thermo_from_pT(mGsEqE);

        // For each mechanism, we need to decide which equilibrium state with which to
        // compute the rate. That decision is stored here
        GasState * gs_eq;

        // compute the rates
        rates[] = to!number(0.0);
        foreach (mech; mEEM) {
            // If mode_j is less than zero, this means energy is being exchanged with the
            // transrotation mode. Therefore we use transrotational equilibrium. Otherwise
            // we'll use vibration/electron equilibrium
            gs_eq = (mech.mode_j < 0) ? &mGsEq : &mGsEqE;
            number rate = mech.rate(gs, *gs_eq, mMolef, mNumden, rmech);
            rates[mech.mode_i] += rate;

            // If the energy is exchanged between modes, subtract the energy
            // from the other mode
            if (mech.mode_j >= 0) rates[mech.mode_j] -= rate;
        }
        // if the electronic energy wants to go negative, stop it
        // assumes electrons are present, and in the last position
        // if ((gs.massf[$-1] < 1e-13 || gs.u_modes[$-1] < 1e-3) && rates[1] < 0.0) rates[1] = 0.0;
    }

private:
    GasState mGsEqE;
}
