/**
 * equilibrium_gas.d
 *
 * A thermally-perfect gas mix with equilibrium thermochemistry.
 *
 * Authors: Hugh McDougall, Peter J. and Rowan G.
 * Version: 2018-05-24: initial shell copied from fuel_air_mix.d.
 */

module gas.equilibrium_gas;

import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;
import gas.therm_perf_gas;
import gas.thermo.cea_thermo_curves;
import gas.thermo.perf_gas_mix_eos;
import gas.thermo.therm_perf_gas_mix_eos;
import gas.diffusion.viscosity;
import gas.diffusion.therm_cond;
import gas.diffusion.cea_viscosity;
import gas.diffusion.cea_therm_cond;
import gas.diffusion.wilke_mixing_viscosity;
import gas.diffusion.wilke_mixing_therm_cond;
import std.math;
import std.stdio;
import std.string;
import std.file;
import std.json;
import std.conv;
import util.lua;
import util.lua_service;
import core.stdc.stdlib : exit;
import nm.bracketing;
import nm.brent; 
import kinetics.reaction_mechanism;
import kinetics.reaction;

// First, the basic gas model.

class EquilibriumGas: ThermallyPerfectGas {
public:
    this(lua_State *L) {
        super(L); 
    } // end constructor

    this(in string fname)
    {
        auto L = init_lua_State();
        doLuaFile(L, fname);
        this(L);
        lua_close(L);
        create_species_reverse_lookup();
    }

    override string toString() const
    {
        char[] repr;
        repr ~= "EquilibriumGas =(species=[TODO]";
        repr ~= ")";
        return to!string(repr);
    }

    override void update_thermo_from_pT(GasState Q) 
    {
        super.update_thermo_from_pT(Q); // [TODO] replace with Hugh's functions
    }
    override void update_thermo_from_rhou(GasState Q)
    {
        super.update_thermo_from_rhou(Q); // [TODO] replace with Hugh's functions
    }
    override void update_thermo_from_rhoT(GasState Q)
    {
        super.update_thermo_from_rhoT(Q); // [TODO] replace with Hugh's functions
    }
    override void update_thermo_from_rhop(GasState Q)
    {
        super.update_thermo_from_rhop(Q); // [TODO] replace with Hugh's functions
    }
    override void update_thermo_from_ps(GasState Q, double s)
    {
        super.update_thermo_from_ps(Q, s);// [TODO] replace with Hugh's functions
    }

    // [TODO] other functions from therm_perf_gas.d, as needed.
    
} // end class EquilibriumGas

// Unit test of the basic gas model...

version(equilibrium_gas_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
        auto gm = new EquilibriumGas("sample-data/therm-perf-5-species-air.lua");
        auto gd = new GasState(5, 0);
        assert(approxEqual(3.621, gm.LJ_sigmas[0]), failedUnitTest());
        assert(approxEqual(97.530, gm.LJ_epsilons[0]), failedUnitTest());
               
        gd.p = 1.0e6;
        gd.T = 2000.0;
        gd.massf = [0.2, 0.2, 0.2, 0.2, 0.2];
        gm.update_thermo_from_pT(gd);
        assert(approxEqual(11801825.6, gd.u, 1.0e-6), failedUnitTest());
        assert(approxEqual(1.2840117, gd.rho, 1.0e-6), failedUnitTest());
        return 0;
    }
}
