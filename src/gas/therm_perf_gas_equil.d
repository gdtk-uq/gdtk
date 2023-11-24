// therm_perf_gas_equil.d

/*
    This is a totally identical subclass of ThermallyPerfectGas designed to invoke
    the EquilibriumUpdate reactor class present in the kinetics folder. Since reactors
    in eilmer 4 are tied to different gas models, we need a different gas model even
    though ThermallyPerfectGas would work just fine with EquilibriumUpdate.

    @author: Nick Gibbons/PJ
*/

module gas.therm_perf_gas_equil;

import std.math;
import std.stdio;
import std.string;
import std.conv : to;
import util.lua;
import util.lua_service;
import ntypes.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.therm_perf_gas;


class ThermallyPerfectGasEquilibrium: ThermallyPerfectGas {
public:
    this(lua_State* L)
    {
        super(L);
        type_str = "ThermallyPerfectGasEquilibrium";
    } // end constructor using Lua interpreter

    this(in string fname)
    {
        auto L = init_lua_State();
        doLuaFile(L, fname);
        this(L);
        lua_close(L); // We no longer need the Lua interpreter.
    } // end constructor from a Lua file

} // end class ThermallyPerfectGasEquilibrium

version(therm_perf_gas_equil_test) {
    int main() {
        import util.msg_service;

        FloatingPointControl fpctrl;
        // Enable hardware exceptions for division by zero, overflow to infinity,
        // invalid operations, and uninitialized floating-point variables.
        // Copied from https://dlang.org/library/std/math/floating_point_control.html
        fpctrl.enableExceptions(FloatingPointControl.severeExceptions);
        //
        auto gm = new ThermallyPerfectGasEquilibrium("sample-data/therm-perf-equil-5-species-air.lua");

        auto gd = GasState(5, 0);
        assert(isClose(3.621, gm.LJ_sigmas[0]), failedUnitTest());
        assert(isClose(97.530, gm.LJ_epsilons[0]), failedUnitTest());

        gd.p = 1.0e6;
        gd.T = 2000.0;
        gd.massf = [0.2, 0.2, 0.2, 0.2, 0.2];
        gm.update_thermo_from_pT(gd);
        assert(isClose(11801825.6, gd.u, 1.0e-6), failedUnitTest());
        assert(isClose(1.2840117, gd.rho, 1.0e-6), failedUnitTest());

        gd.rho = 2.0;
        gd.u = 14.0e6;
        gm.update_thermo_from_rhou(gd);
        assert(isClose(3373757.4, gd.p, 1.0e-6), failedUnitTest());
        assert(isClose(4331.944, gd.T, 1.0e-6), failedUnitTest());

        gd.T = 10000.0;
        gd.rho = 1.5;
        gm.update_thermo_from_rhoT(gd);
        assert(isClose(5841068.3, gd.p, 1.0e-6), failedUnitTest());
        assert(isClose(20340105.9, gd.u, 1.0e-6), failedUnitTest());

        gd.rho = 10.0;
        gd.p = 5.0e6;
        gm.update_thermo_from_rhop(gd);
        assert(isClose(11164648.5, gd.u, 1.0e-6), failedUnitTest());
        assert(isClose(1284.012, gd.T, 1.0e-6), failedUnitTest());

        gd.p = 1.0e6;
        number s = 10000.0;
        gm.update_thermo_from_ps(gd, s);
        assert(isClose(2560.118, gd.T, 1.0e-6), failedUnitTest());
        assert(isClose(12313952.52, gd.u, 1.0e-6), failedUnitTest());
        assert(isClose(1.00309, gd.rho, 1.0e-5), failedUnitTest());

        s = 11000.0;
        number h = 17.0e6;
        gm.update_thermo_from_hs(gd, h, s);
        assert(isClose(5273.103, gd.T, 1.0e-5), failedUnitTest());
        assert(isClose(14946629.7, gd.u, 1.0e-5), failedUnitTest());
        assert(isClose(0.4603513, gd.rho, 1.0e-5), failedUnitTest());
        assert(isClose(945271.84, gd.p, 1.0e-4), failedUnitTest());

        gd.T = 4000.0;
        gm.update_trans_coeffs(gd);
        assert(isClose(0.00012591, gd.mu, 1.0e-5), failedUnitTest());
        assert(isClose(0.2448263, gd.k, 1.0e-5), failedUnitTest());

        // [TODO]
        // entropy, enthalpy and sound speed tests
        return 0;
    } // end main()
}
