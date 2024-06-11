// init_gas_model.d
// Utility function to make a specific GasModel
// from information contained in a Lua script.
//

module gas.init_gas_model;

import std.format;
import util.lua;
import util.lua_service;
import util.msg_service;

import gas.gas_model;
// Needs to know about all of the gas-model modules that are in play.
import gas.ideal_gas;
import gas.ideal_helium;
import gas.cubic_gas;
import gas.cea_gas;
import gas.composite_gas;
import gas.therm_perf_gas;
import gas.therm_perf_gas_equil;
import gas.very_viscous_air;
import gas.uniform_lut;
import gas.uniform_lut_plus_ideal;
import gas.adaptive_lut_CEA;
import gas.ideal_gas_ab;
import gas.two_temperature_reacting_argon;
import gas.two_temperature_argon_plus_ideal;
import gas.ideal_dissociating_gas;
import gas.two_temperature_air;
import gas.two_temperature_nitrogen;
import gas.two_temperature_dissociating_nitrogen : TwoTemperatureDissociatingNitrogen;
import gas.vib_specific_nitrogen;
import gas.vib_specific_co;
import gas.fuel_air_mix;
import gas.equilibrium_gas;
// RJG, disabled while I figure out the cyclic dependency in static constructors
// import gas.electronically_specific_gas: ElectronicallySpecificGas;
import gas.two_temperature_gasgiant: TwoTemperatureGasGiant;

GasModel init_gas_model(string file_name="gas-model.lua")
/**
 * Get the instructions for setting up the GasModel object from a Lua file.
 * The first item in the file should be a model name which we use to select
 * the specific GasModel class.
 * The constructor for each specific gas model will know how to pick out the
 * specific parameters of interest.
 * As new GasModel classes are added to the collection, just
 * add a new case to the switch statement below.
 */
{
    lua_State* L;

    try {
        L = init_lua_State();
        doLuaFile(L, file_name);
    } catch (LuaInputException e) {
        string msg = "In function init_gas_model() in gas_model.d";
        msg ~= format(" there was a problem parsing the input file: %s\n", file_name);
        msg ~= "There could be a Lua syntax error OR the file might not exist.\n";
        msg ~= "Lua error message follows:\n";
        msg ~= e.msg;
        throw new GasModelException(msg);
    }
    string gas_model_name;
    try {
        gas_model_name = getString(L, "model");
    } catch (Exception e) {
        string msg = "In function init_gas_model() in gas_model.d, ";
        msg ~= "there was a problem reading the 'model' name";
        msg ~= format(" from the gas model input Lua file: %s", file_name);
        throw new GasModelException(msg);
    }
    GasModel gm;
    try {
        switch (gas_model_name) {
        case "IdealGas":
            gm = new IdealGas(L);
            break;
        case "IdealHelium":
            gm = new IdealHelium();
            break;
        case "CubicGas":
            gm = new CubicGas(L);
            break;
        case "CEAGas":
            gm = new CEAGas(L);
            break;
        case "CompositeGas":
            gm = new CompositeGas(L);
            break;
        case "ThermallyPerfectGas":
            gm = new ThermallyPerfectGas(L);
            // There is a special step if we're using the Gas Giant transport properties model.
            lua_getglobal(L, "use_gas_giant_transport_properties");
            if (!lua_isnil(L, -1)) {
                ThermallyPerfectGas tpg = cast(ThermallyPerfectGas) gm;
                tpg.attachGasModelToGasGiantModel();
            }
            lua_pop(L, 1);
            break;
        case "ThermallyPerfectGasEquilibrium":
            gm = new ThermallyPerfectGasEquilibrium(L);
            break;
        case "VeryViscousAir":
            gm = new VeryViscousAir(L);
            break;
        case "look-up table":
            gm = new UniformLUT(L);
            break;
        case "UniformLUTPlusIdealGas":
            gm = new UniformLUTPlusIdealGas(L);
            break;
        case "CEA adaptive look-up table":
            gm = new AdaptiveLUT(L);
            break;
        case "IdealGasAB":
            gm = new IdealGasAB(L);
            break;
        case "TwoTemperatureReactingArgon":
            gm = new TwoTemperatureReactingArgon(L);
            break;
        case "TwoTemperatureArgonPlusIdealGas":
            gm = new TwoTemperatureArgonPlusIdealGas(L);
            break;
        case "IdealDissociatingGas":
            gm = new IdealDissociatingGas(L);
            break;
        case "TwoTemperatureAir":
            gm = new TwoTemperatureAir(L);
            break;
        case "TwoTemperatureNitrogen":
            gm = new TwoTemperatureNitrogen();
            break;
        case "TwoTemperatureDissociatingNitrogen":
            gm = new TwoTemperatureDissociatingNitrogen(L);
            break;
        case "VibSpecificNitrogen":
            gm = new VibSpecificNitrogen(L);
            break;
        case "VibSpecificCO":
            gm = new VibSpecificCO(L);
            break;
        case "VibSpecificCOMixture":
            gm = new VibSpecificCOMixture(L);
            break;
        case "FuelAirMix":
            gm = new FuelAirMix(L);
            break;
        case "EquilibriumGas":
            gm = new EquilibriumGas(L);
            break;
        case "TwoTemperatureGasGiant":
            gm = new TwoTemperatureGasGiant();
            break;
        default:
            string errMsg = format("The gas model '%s' is not available.", gas_model_name);
            throw new Error(errMsg);
        }
    } catch (Error e) {
        string msg = "Error while trying to initialize gas model " ~ gas_model_name ~ "\n";
        msg ~= e.msg;
        throw new Error(msg);
    }
    lua_close(L);
    return gm;
} // end init_gas_model()


version(init_gas_model_test) {
    int main() {
        import std.math;
        import std.stdio;
        import std.conv;
        import ntypes.complex;
        import nm.number;

        // Methods for testing gas state class
        import gas.gas_state;
        GasState* gd = new GasState(2, 1); // Use a pointer because we will change the GasModel later.
        gd.massf[0] = 0.8;
        gd.massf[1] = 0.2;
        number[] phi = [to!number(9.0), to!number(16.0)];
        assert(approxEqualNumbers(to!number(10.4), mass_average(*gd, phi), 1.0e-6));

        // Iterative methods test using idealgas single species model
        // These assume IdealGas class is working properly
        import gas.init_gas_model: init_gas_model;
        GasModel gm;
        try {
            gm = init_gas_model("sample-data/ideal-air-gas-model.lua");
        }
        catch (Exception e) {
            writeln(e.msg);
            string msg;
            msg ~= "Test of iterative methods in gas_model.d require file:";
            msg ~= " ideal-air-gas-model.lua in directory: gas/sample_data";
            throw new Exception(msg);
        }

        gd = new GasState(gm, 100.0e3, 300.0);
        assert(approxEqualNumbers(gm.R(*gd), to!number(287.086), 1.0e-4), "gas constant");
        assert(gm.n_modes == 0, "number of energy modes");
        assert(gm.n_species == 1, "number of species");
        assert(approxEqualNumbers(gd.p, to!number(1.0e5)), "pressure");
        assert(approxEqualNumbers(gd.T, to!number(300.0), 1.0e-6), "static temperature");
        assert(approxEqualNumbers(gd.massf[0], to!number(1.0), 1.0e-6), "massf[0]");

        gm.update_thermo_from_pT(*gd);
        gm.update_sound_speed(*gd);
        assert(approxEqualNumbers(gd.rho, to!number(1.16109), 1.0e-4), "density");
        assert(approxEqualNumbers(gd.u, to!number(215314.0), 1.0e-4), "internal energy");
        assert(approxEqualNumbers(gd.a, to!number(347.241), 1.0e-4), "sound speed");
        gm.update_trans_coeffs(*gd);
        assert(approxEqualNumbers(gd.mu, to!number(1.84691e-05), 1.0e-6), "viscosity");
        assert(approxEqualNumbers(gd.k, to!number(0.0262449), 1.0e-6), "conductivity");

        // Select arbitrary energy and density and establish a set of
        // variables that are thermodynamically consistent
        number e_given = 1.0e7;
        number rho_given = 2.0;
        GasState Q = GasState(gm); // Use the struct directly.
        Q.u = e_given;
        Q.rho = rho_given;
        gm.update_thermo_from_rhou(Q);
        number p_given = Q.p;
        number T_given = Q.T;

        // Initialise the same state from the different property combinations
        // Test pT iterative update
        Q.p = p_given;
        Q.T = T_given;
        update_thermo_state_pT(gm, Q);
        // Determine correct entropy/enthalpy for updates that use them
        number s_given = gm.entropy(Q);
        number h_given = gm.enthalpy(Q);
        assert(approxEqualNumbers(Q.rho, rho_given, 1.0e-6),  failedUnitTest());
        assert(approxEqualNumbers(Q.u, e_given, 1.0e-6), failedUnitTest());
        // Test rhoT iterative update
        Q.rho = rho_given;
        Q.T = T_given;
        update_thermo_state_rhoT(gm, Q);
        assert(approxEqualNumbers(Q.u, e_given, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(Q.p, p_given, 1.0e-6),  failedUnitTest());
        // Test rhop iterative update
        Q.rho = rho_given;
        Q.p = p_given;
        assert(approxEqualNumbers(Q.T, T_given, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(Q.u, e_given, 1.0e-6), failedUnitTest());
        // Test  ps iterative update
        Q.p = p_given;
        update_thermo_state_ps(gm, Q, s_given);
        assert(approxEqualNumbers(Q.T, T_given, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(Q.u, e_given, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(Q.rho, rho_given, 1.0e-6), failedUnitTest());
        // Test hs iterative update
        assert(approxEqualNumbers(Q.T, T_given, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(Q.u, e_given, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(Q.rho, rho_given, 1.0e-6), failedUnitTest());
        assert(approxEqualNumbers(Q.p, p_given, 1.0e-6), failedUnitTest());

        version(complex_numbers) {
            // Check du/dT = Cv
            number u0 = Q.u;
            double h = 1.0e-20;
            Q.T += complex(0.0,h);
            update_thermo_state_rhoT(gm, Q);
            double myCv = Q.u.im/h;
            assert(isClose(myCv, gm.dudT_const_v(Q).re), failedUnitTest());
        }
        return 0;
    }
}
