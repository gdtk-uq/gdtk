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

interface EnergyExchangeSystem {
    @nogc void evalRelaxationTimes(in GasState gs);
    @nogc void evalRates(in GasState gs, ref number[] rates);
}

class TwoTemperatureEnergyExchange : EnergyExchangeSystem {
public:

    this(string fname, GasModel gmodel)
    {
	// For 2-T model, one entry in T_modes, so index is 0.
	int mode = 0;
	
        mGmodel = gmodel;
        mGsEq = new GasState(gmodel);

	// Configure other parameters via a Lua state
	auto L = init_lua_State();
	doLuaFile(L, fname);

	// Check on species order before proceeding.
        lua_getglobal(L, "species");
	if (lua_isnil(L, -1)) {
            string errMsg = format("There is no species listing in your kinetics input file: '%s'\n", fname);
            throw new Error(errMsg);
        }
        foreach (isp; 0 .. gmodel.n_species) {
            lua_rawgeti(L, -1, isp);
            auto sp = to!string(luaL_checkstring(L, -1));
            if (sp != gmodel.species_name(isp)) {
                string errMsg = "Species order is incompatible between gas model and kinetics input.\n";
                errMsg ~= format("Kinetics input file is: '%s'.\n", fname);
                errMsg ~= format("In gas model: index %d ==> species '%s'\n", isp, gmodel.species_name(isp));
                errMsg ~= format("In kinetics: index %d ==> species '%s'\n", isp, sp);
                throw new Error(errMsg);
            }
            lua_pop(L, 1);
        }
        lua_pop(L, 1);

	lua_getglobal(L, "vibrational_relaxers");
	auto nVib = to!int(lua_objlen(L, -1));
	mVibRelaxers.length = nVib;
	foreach (int isp; 0 .. nVib) {
	    lua_rawgeti(L, -1, isp+1);
	    auto sp = to!string(luaL_checkstring(L, -1));
	    mVibRelaxers[isp] = gmodel.species_index(sp);
	    lua_pop(L, 1);
	}
	lua_pop(L, 1);

	mMolef.length = gmodel.n_species;
	// We make mVT of length n_species for convenience, even though
	// we'll only fill array entries associated vibRelaxers.
	// This just makes indexing consistent, and avoid having
	// a separate indexing just for the vibRelaxers.
        mVT.length = gmodel.n_species;

	lua_getglobal(L, "mechanism");
	
	foreach (ip; mVibRelaxers) {
	    mVT[ip].length = gmodel.n_species;
	    foreach (iq; 0 .. gmodel.n_species) {
		auto p = gmodel.species_name(ip);
		auto q = gmodel.species_name(iq);
		string key = p ~ ":" ~ q ~ "|VT";
		lua_getfield(L, -1, key.toStringz);
                if (!lua_isnil(L, -1)) {
                    mVT[ip][iq] = createVTMechanism(L, ip, iq, mode, gmodel);
                    mHeavyParticles ~= iq;
                }
		lua_pop(L, 1);
	    }
	}
	lua_pop(L, 1);
    }
    
    @nogc
    void evalRelaxationTimes(in GasState gs)
    {
        mGmodel.massf2molef(gs, mMolef);
        foreach (p; mVibRelaxers) {
            foreach (q; mHeavyParticles) {
                mVT[p][q].evalRelaxationTime(gs, mMolef);
            }
        }
    }
    
    @nogc
    void evalRates(in GasState gs, ref number[] rates)
    {
        // Compute a star state at transrotational equilibrium.
        mGsEq.copy_values_from(gs);
        mGsEq.T_modes[0] = gs.T;
        mGmodel.update_thermo_from_pT(mGsEq);
        number rate = 0.0;
        // Compute rate change from VT exchange
        foreach (p; mVibRelaxers) {
            foreach (q; mHeavyParticles) {
                rate += mVT[p][q].rate(gs, mGsEq, mMolef);
            }
        }
        rates[0] = rate;
        // TODO: Compute rate change from ET exchange
        // TODO: Compute rate change from chemistry coupling.
    }

private:
    int[] mVibRelaxers;
    int[] mHeavyParticles;
    number[] mMolef;
    GasModel mGmodel;
    GasState mGsEq;
    EnergyExchangeMechanism[][] mVT;
    // EnergyExchangeMechanism[] mET;
    // EnergyExchangeMechanism[] mChemCoupling;
    
}

