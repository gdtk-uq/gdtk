/*
 * gas_mixtures.d
 *
 * Author: Rowan G.
 * Version: 2021-03-01
 */

module gas.diffusion.gas_mixtures;

import std.math;
import std.string;
import std.conv;

import util.lua;
import util.lua_service;
import nm.number;

import gas;
import gas.diffusion.transport_properties_model;
import gas.diffusion.viscosity;
import gas.diffusion.therm_cond;
import gas.diffusion.wilke_mixing_viscosity;
import gas.diffusion.wilke_mixing_therm_cond;

class GasMixtureTransProps : TransportPropertiesModel {
public:

    this(lua_State *L, string[] speciesNames)
    {
        Viscosity[] vms;
        ThermalConductivity[] tcms;
        double[] molMasses;
        foreach (isp, spName; speciesNames) {
            lua_getglobal(L, "db");
            lua_getfield(L, -1, spName.toStringz);

            molMasses ~= getDouble(L, -1, "M");

            lua_getfield(L, -1, "viscosity");
            vms ~= createViscosityModel(L);
            lua_pop(L, 1);

            lua_getfield(L, -1, "thermal_conductivity");
            tcms ~=  createThermalConductivityModel(L);
            lua_pop(L, 1);

            lua_pop(L, 1);
            lua_pop(L, 1);
        }
        mVisc = new WilkeMixingViscosity(vms, molMasses);
        mThermCond = new WilkeMixingThermCond(tcms, molMasses);
    }
    
    @nogc
    override void updateTransProps(GasState gs)
    {
        mVisc.update_viscosity(gs);
        mThermCond.update_thermal_conductivity(gs);
    }

private:
    WilkeMixingViscosity mVisc;
    WilkeMixingThermCond mThermCond;
}


version(gas_mixtures_test) {
    int main() {
        import util.msg_service;

        FloatingPointControl fpctrl;
        // Enable hardware exceptions for division by zero, overflow to infinity,
        // invalid operations, and uninitialized floating-point variables.
        // Copied from https://dlang.org/library/std/math/floating_point_control.html
        fpctrl.enableExceptions(FloatingPointControl.severeExceptions);

        auto L = init_lua_State();
        doLuaFile(L, "sample-data/therm-perf-5-species-air.lua");
        string[] speciesNames;
        getArrayOfStrings(L, LUA_GLOBALSINDEX, "species", speciesNames);
        auto tp = new GasMixtureTransProps(L, speciesNames);
        lua_close(L);

        auto gs = new GasState(5, 0);
        gs.p = 1.0e6;
        gs.T = 4000.0;
        gs.massf = [to!number(0.2), to!number(0.2), to!number(0.2), to!number(0.2), to!number(0.2)];
        tp.updateTransProps(gs);
        assert(approxEqualNumbers(to!number(0.00012591), gs.mu, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(to!number(0.2448263), gs.k, 1.0e-6), failedUnitTest());

        return 0;
    }
}
