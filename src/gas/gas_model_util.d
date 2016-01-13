/**
 * gasmodelutil.d
 * Utility functions that make use of the gasmodel class and its derived classes.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-06-22, first cut, exploring the options.
 */

module gas.gas_model_util;

import gas.gas_model;
import gas.ideal_gas;
import gas.therm_perf_gas;
import gas.very_viscous_air;
import gas.co2gas;
import gas.co2gas_sw;
import gas.sf6virial;
import std.file;
import std.stdio;
import util.lua;
import util.lua_service;
import core.stdc.stdlib : exit;

/**
 * We get the instructions for setting up the GasModel object
 * from a Lua file.  The first item should be a model name
 * which we use to select the specific GasModel class.
 * We should then find a table in the file that corresponds
 * to the model name. We will get this table and construct
 * a specific object based on that table. An init() function
 * associated with each class will know how to pick out the
 * specific parameters of interest.
 * As new GasModel classes are added to the collection, just 
 * add a new case to the switch statement below.
 */
GasModel init_gas_model(in string file_name="gas-model.lua") {
    lua_State* L;
    try { 
        L = init_lua_State(file_name);
    } catch (Exception e) {
        writeln("ERROR: in function init_gas_model() in gasmodelutil.d");
        writeln("ERROR: There was a problem parsing the input file: ", file_name);
	writeln("ERROR: There could be a Lua syntax error OR the file might not exist.");
	writeln("ERROR: Quitting at this point.");
 	exit(1);
    }
    string gas_model_name;
    try {
    	gas_model_name = getString(L, LUA_GLOBALSINDEX, "model");
    } catch (Exception e) {
        writeln("ERROR: in function init_gas_model() in gas_model_util.d");
        writeln("ERROR: There was a problem reading the 'model' name" );
	writeln("ERROR: in the gas model input Lua file.");
	writeln("ERROR: Quitting at this point.");
        exit(1);
    }
    GasModel gm;
    switch ( gas_model_name ) {
    case "IdealGas":
	gm = new IdealGas(L);
		break;
    case "ThermallyPerfectGas":
	gm = new ThermallyPerfectGas(L);
	break;
    case "VeryViscousAir":
	gm = new VeryViscousAir(L);
	break;
    case "CO2Gas":
	gm = new CO2Gas(L);
	break;
    case "CO2GasSW":
	gm = new CO2GasSW(L);
	break;
    case "SF6Virial":
    	gm = new SF6Virial(L);
    	break;
    default:
	gm = new IdealGas();
    }
    lua_close(L);
    return gm;
}


version(gas_model_util_test) {
    int main() {
	import std.math;
	auto gm = init_gas_model("sample-data/ideal-air-gas-model.lua");
	auto gd = new GasState(gm, 100.0e3, 300.0);
	assert(approxEqual(gm.R(gd), 287.086, 1.0e-4), "gas constant");
	assert(gm.n_modes == 1, "number of energy modes");
	assert(gm.n_species == 1, "number of species");
	assert(approxEqual(gd.p, 1.0e5), "pressure");
	assert(approxEqual(gd.T[0], 300.0, 1.0e-6), "static temperature");
	assert(approxEqual(gd.massf[0], 1.0, 1.0e-6), "massf[0]");

	gm.update_thermo_from_pT(gd);
	gm.update_sound_speed(gd);
	assert(approxEqual(gd.rho, 1.16109, 1.0e-4), "density");
	assert(approxEqual(gd.e[0], 215314.0, 1.0e-4), "internal energy");
	assert(approxEqual(gd.a, 347.241, 1.0e-4), "sound speed");
	gm.update_trans_coeffs(gd);
	assert(approxEqual(gd.mu, 1.84691e-05, 1.0e-6), "viscosity");
	assert(approxEqual(gd.k[0], 0.0262449, 1.0e-6), "conductivity");

	return 0;
    }
}
