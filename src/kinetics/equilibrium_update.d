/**
 * equilibrium_update.d
 *
 * Thermally-perfect gas mix in equilibrium.
 *
 * This kinetics file accompanies the gas model in gas/therm_perf_gas_equil.d
 *
 * Authors: Nick Gibbons.
 * Version: 2020-xx-xx
 */

module kinetics.equilibrium_update;

import std.math;
import nm.complex;
import nm.number;

import gas;
import util.lua;
import util.lua_service;
import kinetics.thermochemical_reactor;

final class EquilibriumUpdate : ThermochemicalReactor {
    
    this(string fname, GasModel gmodel)
    {
        super(gmodel); // hang on to a reference to the gas model
        // We may need to pick a number of pieces out of the file.
        auto L = init_lua_State();
        doLuaFile(L, fname);
        lua_getglobal(L, "config");
        // Now, pull out the numeric value parameters.
        // _alpha = getDouble(L, -1, "alpha");
        lua_pop(L, 1); // dispose of the table
        lua_close(L);
    }
    
    override void opCall(GasState Q, double tInterval,
                         ref double dtChemSuggest, ref double dtThermSuggest, 
                         ref number[maxParams] params)
    {
        // Nick, the simplest code to look at, as an example,
        // is in the powers_aslam_kinetics.d file.

        // Since the internal energy and density in the (isolated) reactor is fixed,
        // we need to evaluate the new temperature, pressure, etc.
        _gmodel.update_thermo_from_rhou(Q);
        _gmodel.update_sound_speed(Q);
    }

private:
    // Maybe some parameters stored.
} // end class EquilibriumUpdate


version(equilibrium_update_test) {
    import std.stdio;
    import util.msg_service;
    import gas.therm_perf_gas_equil;

    int main() {
        auto gm = new ThermallyPerfectGasEquilibrium("../gas/sample-data/therm-perf-equil-5-species-air.lua");
        auto gd = new GasState(5, 0);

        gd.p = 0.1*101.35e3;
        gd.T = 2500.0;
        gd.massf = [0.74311527, 0.25688473, 0.0, 0.0, 0.0];
        gm.update_thermo_from_pT(gd);
        assert(approxEqual(0.7321963 , gd.massf[0], 1.0e-6)); 
        assert(approxEqual(0.23281198, gd.massf[1], 1.0e-6));
        assert(approxEqual(0.0, gd.massf[2], 1.0e-6));
        assert(approxEqual(0.01160037, gd.massf[3], 1.0e-6));
        assert(approxEqual(0.02339135, gd.massf[4], 1.0e-6));

        return 0;
    }
}
