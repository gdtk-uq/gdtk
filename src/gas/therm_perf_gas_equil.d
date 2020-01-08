// therm_perf_gas_equil.d
// Thermally-perfect gas with equilibrium chemistry.

module gas.therm_perf_gas_equil;

import std.math;
import std.stdio;
import std.string;
import std.conv : to;
import std.algorithm;
import util.lua;
import util.lua_service;
import nm.complex;
import nm.number;

import gas.gas_model;
import gas.gas_state;
import gas.therm_perf_gas;

class ThermallyPerfectGasEquilibrium: ThermallyPerfectGas {
public:
    this(lua_State* L)
    // Construct the model from parameters that are contained in a Lua interpreter,
    // but delegate all of the hard work to Rowan's ThermallyPerfectGas.
    {
        super(L);
    } // end constructor using Lua interpreter

    this(in string fname)
    {
        super(fname);
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
        auto gm = new ThermallyPerfectGasEquilibrium("sample-data/therm-perf-5-species-air.lua");
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
    } // end main()
}
