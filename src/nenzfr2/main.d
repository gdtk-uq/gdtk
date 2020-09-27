// main.d for nenzfr2
// PJ 2020-09-26 Initial code built from Python prototype.

import std.stdio;
import std.array;
import std.string;
import std.conv;
import dyaml;
import nm.secant;
import gas;
import kinetics;
import gasflow;


int main()
{
    writeln("NENZFR2: shock-tunnel with nonequilibrium nozzle flow.");
    writeln("Revision: PUT_REVISION_STRING_HERE");
    writeln("Compiler-name: PUT_COMPILER_NAME_HERE");
    //
    // Read the input.
    auto config = dyaml.Loader.fromFile("input.yaml").load();
    // Extract our parameters and write them to stdout.
    writeln("Input data.");
    writeln("  "~config["title"].as!string);
    //
    // The first gas model is for the shock-tube analysis,
    // assuming frozen reactions or full equilibrium.
    string gm1_filename = config["gas-model-1"].as!string;
    writeln("  gas-model-1= ", gm1_filename);
    auto gm1 = init_gas_model(gm1_filename);
    //
    // The second gas model is for the expansion process that has finite-rate 1T chemistry.
    string gm2_filename = config["gas-model-2"].as!string;
    writeln("  gas-model-2= ", gm2_filename);
    string reactions_filename = config["reactions"].as!string;
    string reactions_filename2 = "";
    try {
        reactions_filename2 = config["reactions-file2"].as!string;
    } catch (YAMLException e) {
        // Do nothing, but assume 1T chemistry,  if we cannot set the second reactions file.
    }
    writeln("  reactions-file= ", reactions_filename);
    writeln("  reactions_file2= ", reactions_filename2);
    auto gm2 = init_gas_model(gm2_filename);
    auto reactor = init_thermochemical_reactor(gm2, reactions_filename, reactions_filename2);
    //
    string[] species;
    foreach(string name; config["species"]) { species ~= name; }
    writeln("  species= ", species);
    double[] molef;
    foreach(name; species) {
        double mf = 0.0;
        try {
            mf = to!double(config["molef"][name].as!string);
        } catch (YAMLException e) {
            // Assume 0.0.
        }
        molef ~= mf;
    }
    writeln("  molef= ", molef);
    //
    // Initial gas state in shock tube.
    double T1 = to!double(config["T1"].as!string);
    double p1 = to!double(config["p1"].as!string);
    double Vs = to!double(config["Vs"].as!string);
    double pe = 0.0;
    try {
        pe = to!double(config["pe"].as!string);
    } catch (YAMLException e) {
        // Assume 0.0.
    }
    double ar = 1.0;
    try {
        ar = to!double(config["ar"].as!string);
    } catch (YAMLException e) {
        // Assume 1.0.
    }
    writeln("  T1= ", T1);
    writeln("  p1= ", p1);
    writeln("  Vs= ", Vs);
    writeln("  pe= ", pe);
    writeln("  ar= ", ar);
    //
    // Set up equilibrium-gas flow analysis of shock tube.
    // Let's assume a cea2 gas model.
    writeln("Initial gas state.");
    GasState state1 = new GasState(gm1);
    state1.p = p1; state1.T = T1; state1.massf = [1.0,];
    gm1.update_thermo_from_pT(state1);
    gm1.update_sound_speed(state1);
    writeln("  state1: ", state1);
    double H1 = gm1.internal_energy(state1) + state1.p/state1.rho;
    writeln("  H1= ", H1);
    //
    writeln("Start incident-shock calculation.");
    GasState state2 = new GasState(gm1);
    double[2] velocities = normal_shock(state1, Vs, state2, gm1);
    double V2 = velocities[0];
    double Vg = velocities[1];
    writeln("  V2= ", V2, " Vg= ", Vg);
    writeln("  state2: ", state2);
    //
    writeln("Start reflected-shock calculation.");
    GasState state5 = new GasState(gm1);
    double Vr = reflected_shock(state2, Vg, state5, gm1);
    //
    writeln("Start calculation of isentropic relaxation.");
    GasState state5s = new GasState(gm1);
    state5s.copy_values_from(state5);
    // Entropy is set, then pressure is relaxed via an isentropic process.
    state5s.p = (pe > 0.0) ? pe : state5.p;
    double entropy5 = gm1.entropy(state5);
    writeln("  state5.entropy= ", entropy5);
    gm1.update_thermo_from_ps(state5s, entropy5);
    writeln("  state5s: ", state5s);
    double H5s = gm1.internal_energy(state5s) + state5s.p/state5s.rho; // stagnation enthalpy
    writeln("  H5s= ", H5s);
    //
    writeln("Start isentropic relaxation to throat (Mach 1)");
    double error_at_throat(double x)
    {
        // Returns Mach number error as pressure is changed.
        GasState state = new GasState(gm1);
        double V = expand_from_stagnation(state5s, x, state, gm1);
        gm1.update_sound_speed(state);
        return (V/state.a) - 1.0;
    }
    double x6 = 1.0;
    try {
        x6 = nm.secant.solve!(error_at_throat, double)(0.95, 0.90, 1.0e-4);
    } catch (Exception e) {
        writeln("Failed to find throat conditions iteratively.");
    }
    GasState state6 = new GasState(gm1);
    double V6 = expand_from_stagnation(state5s, x6, state6, gm1);
    double mflux6 = state6.rho * V6;  // mass flux per unit area, at throat
    writeln("  state6= ", state6);
    writeln("  V6= ", V6);
    writeln("  mflux6= ", mflux6);
    //
    writeln("Start isentropic relaxation to nozzle exit of given area.");
    // The mass flux going through the nozzle exit has to be the same
    // as that going through the nozzle throat.
    double error_at_exit(double x)
    {
        // Returns mass_flux error as for a given exit pressure."
        GasState state = new GasState(gm1);
        double V = expand_from_stagnation(state5s, x, state, gm1);
        double mflux = state.rho * V * ar;
        return (mflux-mflux6)/mflux6;
    }
    // It appears that we need a pretty good starting guess for the pressure ratio.
    // Maybe a low value is OK.
    double x7 = x6;
    try {
        x7 = nm.secant.solve!(error_at_exit, double)(0.001*x6, 0.00005*x6, 1.0e-4, 1.0/state5s.p, 1.0);
    } catch (Exception e) {
        writeln("Failed to find exit conditions iteratively.");
        x7 = x6;
    }
    GasState state7 = new GasState(gm1);
    double V7 = expand_from_stagnation(state5s, x7, state7, gm1);
    double mflux7 = state7.rho * V7 * ar;
    writeln("  area_ratio= ", ar);
    GasState state7_pitot = new GasState(gm1);
    pitot_condition(state7, V7, state7_pitot, gm1);
    writeln("  state7= ", state7);
    writeln("  V7= ", V7);
    writeln("  mflux7= ", mflux7);
    writeln("  pitot7= ", state7_pitot.p);
    //
    return 0;
} // end main
