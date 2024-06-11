// init_thermochemical_reactor.d
//
// Utility function to make a specific thermochemical updater.
// Each reactor is associated with a particular gas model,
// so we use the gas model as the key to decide which reactor
// to construct.
//

module kinetics.init_thermochemical_reactor;

import ntypes.complex;
import nm.number;

import gas;
import kinetics.thermochemical_reactor;
import std.stdio;

// We need to know about the schemes that are available.
import gas.composite_gas;
import gas.therm_perf_gas;
import gas.therm_perf_gas_equil;
import gas.ideal_gas_ab;
import gas.two_temperature_reacting_argon;
import gas.two_temperature_argon_plus_ideal;
import gas.ideal_dissociating_gas;
import gas.fuel_air_mix;
import gas.two_temperature_nitrogen;
import gas.two_temperature_dissociating_nitrogen;
import gas.vib_specific_nitrogen;
import gas.vib_specific_co;
import gas.two_temperature_air;
import gas.two_temperature_gasgiant;

import kinetics.chemistry_update;
import kinetics.equilibrium_update;
import kinetics.powers_aslam_kinetics;
import kinetics.yee_kotov_kinetics;
import kinetics.multi_temperature_thermochemical_reactor;
import kinetics.two_temperature_argon_kinetics;
import kinetics.ideal_dissociating_gas_kinetics;
import kinetics.fuel_air_mix_kinetics;
import kinetics.two_temperature_nitrogen_kinetics;
import kinetics.two_temperature_dissociating_nitrogen_kinetics;
import kinetics.two_temperature_argon_with_ideal_gas;
import kinetics.vib_specific_nitrogen_kinetics;
import kinetics.vib_specific_co_kinetics;
import kinetics.two_temperature_air_kinetics;
import kinetics.two_temperature_gasgiant_kinetics;


ThermochemicalReactor init_thermochemical_reactor(GasModel gmodel, string fileName1="", string fileName2="")
{
    /*
    Construct a new ThermochemicalReactor object that is appropriate for use with GasModel gmodel.
    Since thermochemical reactors are tied to specific gas models, this routine needs to test
    the actual type of gasmodel that we have, even though OOP normally avoids this.

    In d, this can be done in one of two ways:

    (cast(ObjectType) object) !is null

    or

    typeid(object) is typeid(ObjectType)

    The first is useful for testing for a family of objects, it returns true if object is an ObjectType
    or anything that inherits from ObjectType. The second is useful for testing whether something is
    an actual specific type, excluding parent objects it might be related to. In this routine,
    ThermallyPerfectGas and ThermallyPerfectGasEquilibrium both have different reactor types,
    in spite of being related by inheritance, so the typeid method is used. The other ones
    all use the cast method, since it allows new gas models to be created by inheritance without
    changes being made to this routine.
    */
    ThermochemicalReactor reactor; // start with a null reference

    if (typeid(gmodel) is typeid(ThermallyPerfectGas)) {
        reactor = new ChemistryUpdate(fileName1, gmodel);
    }
    if (typeid(gmodel) is typeid(CompositeGas)) {
        auto cg = cast(CompositeGas) gmodel;
        if (cg.physicalModel == "thermally-perfect-gas") {
            reactor = new ChemistryUpdate(fileName1, gmodel);
        }
        else {
            reactor = new MultiTemperatureThermochemicalReactor(fileName1, fileName2, gmodel);
        }
    }
    if (typeid(gmodel) is typeid(ThermallyPerfectGasEquilibrium)) {
        reactor = new EquilibriumUpdate(fileName1, gmodel);
    }
    if ((cast(IdealGasAB) gmodel) !is null) {
        auto ig = cast(IdealGasAB) gmodel;
        if (ig.modelType == "Yee-Kotov") {
            reactor = new UpdateAB_YeeKotov(fileName1, gmodel);
        }
        else {
            reactor = new UpdateAB(fileName1, gmodel);
        }
    }
    if ((cast(TwoTemperatureReactingArgon) gmodel) !is null) {
        reactor = new UpdateArgonFrac(fileName1, gmodel);
    }
    if ((cast(TwoTemperatureArgonPlusIdealGas) gmodel) !is null) {
        reactor = new UpdateArgonFracWithIdeal(fileName1, gmodel);
    }
    if ((cast(IdealDissociatingGas) gmodel) !is null) {
        reactor = new UpdateIDG(fileName1, gmodel);
    }
    if ((cast(FuelAirMix) gmodel) !is null) {
        reactor = new MixingLimitedUpdate(fileName1, gmodel);
    }
    if ((cast(TwoTemperatureNitrogen) gmodel) !is null) {
        reactor = new VibRelaxNitrogen(fileName1, gmodel);
    }
    if ((cast(TwoTemperatureDissociatingNitrogen) gmodel) !is null) {
        reactor = new TwoTemperatureDissociatingNitrogenKinetics(fileName1, fileName2, gmodel);
    }
    if ((cast(VibSpecificNitrogen) gmodel) !is null) {
        reactor = new VibSpecificNitrogenRelaxation(fileName1, gmodel);
    }
    if (typeid(gmodel) is typeid(VibSpecificCO)) {
        reactor = new VibSpecificCORelaxation(fileName1, gmodel);
    }
    if (typeid(gmodel) is typeid(VibSpecificCOMixture)) {
        reactor = new VibSpecificCOMixtureRelaxation(fileName1, fileName2, gmodel);
    }
    if ((cast(TwoTemperatureAir) gmodel) !is null) {
        reactor = new TwoTemperatureAirKinetics(fileName1, fileName2, gmodel);
    }
    if ((cast(TwoTemperatureGasGiant) gmodel) !is null) {
        reactor = new UpdateGasGiant(gmodel);
    }
    if (reactor is null) {
        throw new ThermochemicalReactorUpdateException("Oops, failed to set up a ThermochemicalReactor.");
    }
    return reactor;
} // end init_thermochemical_reactor()

