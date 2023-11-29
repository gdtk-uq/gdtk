/**
 * fuel_air_mix.d
 *
 * Fuel+Air->Products reacting gas as used by JJ.
 *
 * Authors: JJ Hoste, Peter J. and Rowan G.
 * Version: 2017-04-22: initial shell copied from powers-aslam-gas.
 */

module gas.fuel_air_mix;

import std.math;
import std.stdio;
import std.string;
import std.file;
import std.json;
import std.conv;
import util.lua;
import util.lua_service;
import ntypes.complex;
import nm.number;
import nm.bracketing;
import nm.brent;

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
import kinetics.reaction_mechanism;
import kinetics.reaction;

// First, the basic gas model.

class FuelAirMix: ThermallyPerfectGas {
public:
    this(lua_State *L) {
        super(L);
        type_str = "FuelAirMix";
        lua_getglobal(L, "FuelAirMix");
        // Now, pull out the numeric value parameters.
        _A_edm = getDouble(L, -1, "Aedm");
        _B_edm = getDouble(L, -1, "Bedm");
        _laminar_limit = getBool(L, -1, "lamLimit");
        lua_pop(L, 1); // dispose of the table
    } // end constructor

    override string toString() const
    {
        char[] repr;
        repr ~= "FuelAirMix =(";
        repr ~= " Aedm=" ~ to!string(_A_edm);
        repr ~= ", Bedm=" ~ to!string(_B_edm);
        repr ~= ")";
        return to!string(repr);
    }

private:
    //settings specific to EDM model
    double _A_edm;
    double _B_edm;
    bool _laminar_limit;
} // end class FuelAirMix

// Unit test of the basic gas model...

version(fuel_air_mix_test) {
    import std.stdio;
    import util.msg_service;

    int main() {
        lua_State* L = init_lua_State();
        doLuaFile(L, "sample-data/fuel-air-mix-model.lua");
        auto gm = new FuelAirMix(L);
        lua_close(L);
        auto gd = GasState(2, 0);
        gd.p = 1.0e5;
        gd.T = 300.0;
        gd.massf[0] = 0.75; gd.massf[1] = 0.25;
        /+
         [FIX-ME]
        assert(isClose(gm.R(gd), 287.0, 1.0e-4), failedUnitTest());
        assert(gm.n_modes == 0, failedUnitTest());
        assert(gm.n_species == 2, failedUnitTest());
        assert(isClose(gd.p, 1.0e5, 1.0e-6), failedUnitTest());
        assert(isClose(gd.T, 300.0, 1.0e-6), failedUnitTest());
        assert(isClose(gd.massf[0], 0.75, 1.0e-6), failedUnitTest());
        assert(isClose(gd.massf[1], 0.25, 1.0e-6), failedUnitTest());

        gm.update_thermo_from_pT(gd);
        gm.update_sound_speed(gd);
        double my_rho = 1.0e5 / (287.0 * 300.0);
        assert(isClose(gd.rho, my_rho, 1.0e-4), failedUnitTest());
        double my_Cv = gm.dudT_const_v(gd);
        double my_u = my_Cv*300.0 - 0.25*300000.0;
        assert(isClose(gd.u, my_u, 1.0e-3), failedUnitTest());
        double my_Cp = gm.dhdT_const_p(gd);
        double my_a = sqrt(my_Cp/my_Cv*287.0*300.0);
        assert(isClose(gd.a, my_a, 1.0e-3), failedUnitTest());
        gm.update_trans_coeffs(gd);
        assert(isClose(gd.mu, 0.0, 1.0e-6), failedUnitTest());
        assert(isClose(gd.k, 0.0, 1.0e-6), failedUnitTest());
        +/
        return 0;
    }
}
