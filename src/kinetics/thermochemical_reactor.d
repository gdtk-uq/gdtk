/**
 * thermochemical_reactor.d
 *
 * Author: Peter J
 */

module kinetics.thermochemical_reactor;

import nm.complex;
import nm.number;

import gas;
import gas.therm_perf_gas;
import gas.powers_aslam_gas;
import gas.two_temperature_reacting_argon;
import gas.ideal_dissociating_gas;
import gas.fuel_air_mix;
import gas.two_temperature_nitrogen;
import gas.vib_specific_nitrogen;
import gas.two_temperature_air;

import kinetics.chemistry_update;
import kinetics.powers_aslam_kinetics;
import kinetics.two_temperature_argon_kinetics;
import kinetics.ideal_dissociating_gas_kinetics;
import kinetics.fuel_air_mix_kinetics;
import kinetics.two_temperature_nitrogen_kinetics;
import kinetics.vib_specific_nitrogen_kinetics;
import kinetics.two_temperature_air_kinetics;


class ThermochemicalReactorUpdateException : Exception {
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}

class ThermochemicalReactor {
public:
    this(GasModel gmodel)
    {
        // We need a reference to the original gas model object
        // to update the GasState data at a later time.
        _gmodel = gmodel;
    }

    // All the work happens when calling the concrete object
    // which updates the GasState over the (small) time, tInterval.
    //
    // The array params is there to allow extra information to be passed in.
    // For example, the mixing-limited combustion model by JJ Hoste needs
    // some information about the local flow state beyond the usual gas state.
    abstract void opCall(GasState Q, double tInterval,
                         ref double dtChemSuggest, ref double dtThermSuggest, 
                         ref number[] params);
    
    // We will need to access this referenced model from the Lua functions
    // so it needs to be public.
    GasModel _gmodel;
} // end class ThermochemicalReactor


// The utility function to make new ThermochemicalReactor objects needs to know about
// the schemes that are available.

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
    if (reactor is null) {
        throw new ThermochemicalReactorUpdateException("Oops, failed to set up a ThermochemicalReactor.");
    }
    return reactor;
} // end init_thermochemical_reactor()

