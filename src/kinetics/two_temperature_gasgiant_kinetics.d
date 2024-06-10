/**
 * Two-temperature hydrogen, helium gas model for Gas-Giant entry simulations.
 * Authors: Yu Liu, RG and PJ.
 * Date: 2018-10-18
 *
 * Using description and data from:
 *
 */

module kinetics.two_temperature_gasgiant_kinetics;

import std.math;
import ntypes.complex;
import nm.number;

import gas;
import util.lua;
import util.lua_service;
import kinetics.thermochemical_reactor;

final class UpdateGasGiant : ThermochemicalReactor {

    this(GasModel gmodel)
    {
        super(gmodel); // hang on to a reference to the gas model
    }

    override void opCall(ref GasState Q, double tInterval, ref double dtSuggest,
                         ref number[maxParams] params)
    {
        // Do something to update the state of the gas over time interval tInterval.
        //
        // Since the internal energy and density in the (isolated) reactor is fixed,
        // we need to evaluate the new temperature, pressure, etc.
        _gmodel.update_thermo_from_rhou(Q);
        _gmodel.update_sound_speed(Q);
    }

    @nogc override void eval_source_terms(GasModel gmodel, ref GasState Q, ref number[] source) {
        string errMsg = "eval_source_terms not implemented for two_temperature_gasgiant_kinetics.";
        throw new ThermochemicalReactorUpdateException(errMsg);
    }

private:
    // Put specific data or helper functions here.
} // end class UpdateGasGiant


version(two_temperature_gasgiant_kinetics_test) {
    import std.stdio;
    import util.msg_service;
    import std.math : isClose;
    import gas.two_temperature_gasgiant;
    void main() {
        auto gm = new TwoTemperatureGasGiant();
        auto gd = GasState(6, 1);
        gd.p = 1.0e5;
        gd.T = 310.0;
        gd.T_modes[0] = 310;
        gd.massf[gas.two_temperature_gasgiant.Species.H2] = 1.0;
        gd.massf[gas.two_temperature_gasgiant.Species.H] = 0.0;
        gd.massf[gas.two_temperature_gasgiant.Species.Hplus] = 0.0;
        gd.massf[gas.two_temperature_gasgiant.Species.eminus] = 0.0;
        gd.massf[gas.two_temperature_gasgiant.Species.He] = 0.0;
        gd.massf[gas.two_temperature_gasgiant.Species.Heplus] = 0.0;
        assert(isClose(gm.R(gd), 4124.506, 1.0e-4), failedUnitTest());
        assert(gm.n_modes == 1, failedUnitTest());
        assert(gm.n_species == 6, failedUnitTest());
        assert(isClose(gd.p, 1.0e5, 1.0e-6), failedUnitTest());
        assert(isClose(gd.T, 310.0, 1.0e-6), failedUnitTest());
        assert(isClose(gd.massf[gas.two_temperature_gasgiant.Species.H2], 1.0, 1.0e-6), failedUnitTest());
        assert(isClose(gd.massf[gas.two_temperature_gasgiant.Species.H], 0.0, 1.0e-6), failedUnitTest());
        assert(isClose(gd.massf[gas.two_temperature_gasgiant.Species.Hplus], 0.0, 1.0e-6), failedUnitTest());
        assert(isClose(gd.massf[gas.two_temperature_gasgiant.Species.eminus], 0.0, 1.0e-6), failedUnitTest());
        assert(isClose(gd.massf[gas.two_temperature_gasgiant.Species.He], 0.0, 1.0e-6), failedUnitTest());
        assert(isClose(gd.massf[gas.two_temperature_gasgiant.Species.Heplus], 0.0, 1.0e-6), failedUnitTest());

        gm.update_thermo_from_pT(gd);
        double my_rho = 1.0e5 / (4124.506 * 310.0);
        assert(isClose(gd.rho, my_rho, 1.0e-4), failedUnitTest());
        // Put some discerning tests here.
        assert(isClose(1.0, 1.0, 1.0e2), failedUnitTest());
    }
}
