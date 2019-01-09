// init_thermochemical_reactor.d
//
// Utility function to make a specific thermochemical updater.
// Each reactor is associated with a particular gas model,
// so we use the gas model as the key to decide which reactor
// to construct.
//

module kinetics.init_thermochemical_reactor;

import nm.complex;
import nm.number;

import gas;
import kinetics.thermochemical_reactor;

// We need to know about the schemes that are available.
import gas.therm_perf_gas;
import gas.powers_aslam_gas;
import gas.two_temperature_reacting_argon;
import gas.ideal_dissociating_gas;
import gas.fuel_air_mix;
import gas.two_temperature_nitrogen;
import gas.vib_specific_nitrogen;
import gas.two_temperature_air;
import gas.pseudo_species_gas;
import gas.electronically_specific_gas;
import gas.two_temperature_gasgiant;

import kinetics.chemistry_update;
import kinetics.powers_aslam_kinetics;
import kinetics.two_temperature_argon_kinetics;
import kinetics.ideal_dissociating_gas_kinetics;
import kinetics.fuel_air_mix_kinetics;
import kinetics.two_temperature_nitrogen_kinetics;
import kinetics.vib_specific_nitrogen_kinetics;
import kinetics.two_temperature_air_kinetics;
import kinetics.electronically_specific_kinetics;
version (with_dvode)
{
    import kinetics.pseudo_species_kinetics;
}
import kinetics.two_temperature_gasgiant_kinetics;


ThermochemicalReactor init_thermochemical_reactor(GasModel gmodel, string fileName1="", string fileName2="")
{
    ThermochemicalReactor reactor; // start with a null reference
    if ((cast(ThermallyPerfectGas) gmodel) !is null) {
        reactor = new ChemistryUpdate(fileName1, gmodel);
    }
    if ((cast(PowersAslamGas) gmodel) !is null) {
        reactor = new UpdateAB(fileName1, gmodel);
    }
    if ((cast(TwoTemperatureReactingArgon) gmodel) !is null) {
        reactor = new UpdateArgonFrac(fileName1, gmodel);
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
    if ((cast(VibSpecificNitrogen) gmodel) !is null) {
        reactor = new VibSpecificNitrogenRelaxation(fileName1, gmodel);
    }
    if ((cast(TwoTemperatureAir) gmodel) !is null) {
        reactor = new TwoTemperatureAirKinetics(fileName1, fileName2, gmodel);
    }
    if ((cast(ElectronicallySpecificGas) gmodel) !is null) {
        reactor = new ElectronicallySpecificKinetics(fileName1, gmodel);
    }
    version (with_dvode)
    {
        if ((cast(PseudoSpeciesGas) gmodel) !is null) {
            reactor = new PseudoSpeciesKinetics(gmodel);
        }
    }
    if ((cast(TwoTemperatureGasGiant) gmodel) !is null) {
        reactor = new UpdateGasGiant(gmodel);
    }
    if (reactor is null) {
        throw new ThermochemicalReactorUpdateException("Oops, failed to set up a ThermochemicalReactor.");
    }
    return reactor;
} // end init_thermochemical_reactor()

