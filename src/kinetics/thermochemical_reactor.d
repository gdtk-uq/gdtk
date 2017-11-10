/**
 * thermochemical_reactor.d
 *
 * Author: Peter J
 */

module kinetics.thermochemical_reactor;

import gas;
import gas.therm_perf_gas;
import gas.powers_aslam_gas;
import gas.two_temperature_reacting_argon;
import gas.ideal_dissociating_gas;
import gas.fuel_air_mix;
import gas.two_temperature_nitrogen;
import gas.vib_specific_nitrogen;

import kinetics.chemistry_update;
import kinetics.powers_aslam_kinetics;
import kinetics.two_temperature_argon_kinetics;
import kinetics.ideal_dissociating_gas_kinetics;
import kinetics.fuel_air_mix_kinetics;
import kinetics.two_temperature_nitrogen_kinetics;
import kinetics.vib_specific_nitrogen_kinetics;


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
    abstract void opCall(GasState Q, double tInterval, ref double dtSuggest,
			 ref double[] params);
    
    // We will need to access this referenced model from the Lua functions
    // so it needs to be public.
    GasModel _gmodel;
} // end class ThermochemicalReactor


// The utility function to make new ThermochemicalReactor objects needs to know about
// the schemes that are available.


ThermochemicalReactor init_thermochemical_reactor(GasModel gmodel, string fileName="")
{
    ThermochemicalReactor reactor; // start with a null reference
    if ((cast(ThermallyPerfectGas) gmodel) !is null) {
	reactor = new ChemistryUpdate(fileName, gmodel);
    }
    if ((cast(PowersAslamGas) gmodel) !is null) {
	reactor = new UpdateAB(fileName, gmodel);
    }
    if ((cast(TwoTemperatureReactingArgon) gmodel) !is null) {
	reactor = new UpdateArgonFrac(fileName, gmodel);
    }
    if ((cast(IdealDissociatingGas) gmodel) !is null) {
	reactor = new UpdateIDG(fileName, gmodel);
    }
    if ((cast(FuelAirMix) gmodel) !is null) {
	reactor = new MixingLimitedUpdate(fileName, gmodel);
    }
    if ((cast(TwoTemperatureNitrogen) gmodel) !is null) {
	reactor = new VibRelaxNitrogen(fileName, gmodel);
    }
    if ((cast(VibSpecificNitrogen) gmodel) !is null) {
	reactor = new VibSpecificNitrogenRelaxtion(fileName, gmodel);
    }
    if (reactor is null) {
	throw new ThermochemicalReactorUpdateException("Oops, failed to set up a ThermochemicalReactor.");
    }
    return reactor;
} // end init_thermochemical_reactor()

