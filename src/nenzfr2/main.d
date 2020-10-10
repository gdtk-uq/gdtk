// main.d for nenzfr2
// PJ 2020-09-26 Initial code built from Python prototype.

import std.stdio;
import std.array;
import std.string;
import std.conv;
import std.getopt;
import std.file;
import dyaml;
import nm.secant;
import nm.schedule;
import gas;
import gas.cea_gas;
import gas.therm_perf_gas;
import kinetics;
import gasflow;


int main(string[] args)
{
    int exitFlag = 0; // Presume OK in the beginning.
    // Be careful with the usageMsg string; it has embedded newline characters.
    string usageMsg = "Usage: nenzfr2 <input-file>
Options:
   --verbosity=<int>   0 == very terse output
                       1 == key results printed
                       2 == echo input as well as printing more detailed results
   --help              Print this help message.";
    if (args.length < 2) {
        writeln("Too few arguments. You need to specify the input file.");
        writeln(usageMsg);
        stdout.flush();
        exitFlag = 1;
        return exitFlag;
    }
    int verbosityLevel = 1; // default to commenting on major steps
    bool helpWanted = false;
    try {
        getopt(args,
               "verbosity", &verbosityLevel,
               "help", &helpWanted
               );
    } catch (Exception e) {
        writeln("Problem parsing command-line options.");
        writeln("Arguments not processed: ");
        args = args[1 .. $]; // Dispose of program name in first argument.
        foreach (myarg; args) writeln("    arg: ", myarg);
        writeln(usageMsg);
        stdout.flush();
        exitFlag = 1;
        return exitFlag;
    }
    if (verbosityLevel >= 1) {
        writeln("NENZFR2: shock-tunnel with nonequilibrium nozzle flow.");
    }
    if (verbosityLevel >= 2) {
        writeln("Revision: PUT_REVISION_STRING_HERE");
        writeln("Compiler-name: PUT_COMPILER_NAME_HERE");
    }
    if (helpWanted) {
        writeln(usageMsg);
        stdout.flush();
        exitFlag = 0;
        return exitFlag;
    }
    // Read the input file.
    string inputFile = args[1].strip();
    if (!exists(inputFile)) {
        writeln("Did not find input-file: ", inputFile);
        exitFlag = 2;
        return exitFlag;
    }
    // Extract our job parameters from the YAML input file.
    auto config = dyaml.Loader.fromFile(inputFile).load();
    //
    // The first gas model is for the shock-tube analysis,
    // assuming frozen reactions or full equilibrium.
    string gm1_filename = config["gas-model-1"].as!string;
    auto gm1 = init_gas_model(gm1_filename);
    //
    // The second gas model is for the expansion process that has finite-rate 1T chemistry.
    string gm2_filename = config["gas-model-2"].as!string;
    string reactions_filename = config["reactions"].as!string;
    string reactions_filename2 = "";
    try {
        reactions_filename2 = config["reactions-file2"].as!string;
    } catch (YAMLException e) {
        // Do nothing, but assume 1T chemistry,  if we cannot set the second reactions file.
    }
    auto gm2 = init_gas_model(gm2_filename);
    auto reactor = init_thermochemical_reactor(gm2, reactions_filename, reactions_filename2);
    //
    string[] species;
    foreach(string name; config["species"]) { species ~= name; }
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
    // Nozzle-exit area ratio for terminating expansion.
    double ar = 1.0;
    try {
        ar = to!double(config["ar"].as!string);
    } catch (YAMLException e) {
        // Assume 1.0.
    }
    // Nozzle area-ratio schedule.
    double[] xi; foreach(string val; config["xi"]) { xi ~= to!double(val); }
    double[] ai; foreach(string val; config["ai"]) { ai ~= to!double(val); }
    if (verbosityLevel >= 2) {
        writeln("Input data.");
        writeln("  "~config["title"].as!string);
        writeln("  gas-model-1= ", gm1_filename);
        writeln("  gas-model-2= ", gm2_filename);
        writeln("  reactions-file= ", reactions_filename);
        writeln("  reactions_file2= ", reactions_filename2);
        writeln("  species= ", species);
        writeln("  molef= ", molef);
        writeln("  T1= ", T1);
        writeln("  p1= ", p1);
        writeln("  Vs= ", Vs);
        writeln("  pe= ", pe);
        writeln("  ar= ", ar);
        writeln("  xi= ", xi);
        writeln("  ai= ", ai);
    }
    // Set up equilibrium-gas flow analysis of shock tube.
    // Let's assume a cea2 gas model.
    if (verbosityLevel >= 1) { writeln("Initial gas state."); }
    GasState state1 = new GasState(gm1);
    state1.p = p1; state1.T = T1; state1.massf = [1.0,];
    gm1.update_thermo_from_pT(state1);
    gm1.update_sound_speed(state1);
    double H1 = gm1.internal_energy(state1) + state1.p/state1.rho;
    if (verbosityLevel >= 1) {
        writeln("  state1: ", state1);
        writeln("  H1= ", H1);
    }
    //
    if (verbosityLevel >= 1) { writeln("Start incident-shock calculation."); }
    GasState state2 = new GasState(gm1);
    double[2] velocities = normal_shock(state1, Vs, state2, gm1);
    double V2 = velocities[0];
    double Vg = velocities[1];
    if (verbosityLevel >= 1) {
        writeln("  V2= ", V2, " Vg= ", Vg);
        writeln("  state2: ", state2);
    }
    //
    if (verbosityLevel >= 1) { writeln("Start reflected-shock calculation."); }
    GasState state5 = new GasState(gm1);
    double Vr = reflected_shock(state2, Vg, state5, gm1);
    //
    if (verbosityLevel >= 1) { writeln("Start calculation of isentropic relaxation."); }
    GasState state5s = new GasState(gm1);
    state5s.copy_values_from(state5);
    // Entropy is set, then pressure is relaxed via an isentropic process.
    state5s.p = (pe > 0.0) ? pe : state5.p;
    double entropy5 = gm1.entropy(state5);
    gm1.update_thermo_from_ps(state5s, entropy5);
    double H5s = gm1.internal_energy(state5s) + state5s.p/state5s.rho; // stagnation enthalpy
    if (verbosityLevel >= 1) {
        writeln("  state5.entropy= ", entropy5);
        writeln("  state5s= ", state5s);
        writeln("  H5s= ", H5s, " H5s-H1=", H5s-H1);
    }
    //
    if (verbosityLevel >= 1) { writeln("Start isentropic relaxation to throat (Mach 1)"); }
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
        writeln("  Exception message: %s", e.msg);
        exitFlag = 2;
        return exitFlag;
    }
    GasState state6 = new GasState(gm1);
    double V6 = expand_from_stagnation(state5s, x6, state6, gm1);
    double mflux6 = state6.rho * V6;  // mass flux per unit area, at throat
    if (verbosityLevel >= 1) {
        writeln("  state6= ", state6);
        writeln("  V6= ", V6);
        writeln("  mflux6= ", mflux6);
    }
    //
    if (verbosityLevel >= 1) { writeln("Start isentropic relaxation to nozzle exit of given area."); }
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
        // Note that we have the throat conditions and may proceed
        // with a nonequilibrium expansion calculation.
        x7 = x6;
    }
    GasState state7 = new GasState(gm1);
    double V7 = expand_from_stagnation(state5s, x7, state7, gm1);
    double mflux7 = state7.rho * V7 * ar;
    GasState state7_pitot = new GasState(gm1);
    pitot_condition(state7, V7, state7_pitot, gm1);
    if (verbosityLevel >= 1) {
        writeln("  area_ratio= ", ar);
        writeln("  state7= ", state7);
        writeln("  V7= ", V7);
        writeln("  mflux7= ", mflux7);
        writeln("  pitot7= ", state7_pitot.p);
        writeln("End of stage 1: shock-tube and frozen/eq nozzle analysis.");
    }
    // We will continue with a non-equilibrium chemistry expansion
    // only if we have the correct gas models in play.
    auto gm_cea = cast(CEAGas) gm1;
    if (gm_cea is null) {
        writeln("Cannot continue with nonequilibrium expansion.");
        writeln("  Gas model 1 is not of class CEAGas.");
        exitFlag = 2;
        return exitFlag;
    }
    auto gm_tp = cast(ThermallyPerfectGas) gm2;
    if (gm_tp is null) {
        writeln("Cannot continue with nonequilibrium expansion.");
        writeln("  Gas model 2 is not of class ThermallyPerfectGas.");
        exitFlag = 2;
        return exitFlag;
    }
    //
    if (state6.ceaSavedData is null) {
        exitFlag = 3;
        return exitFlag;
    }
    if (verbosityLevel >= 2) {
        writeln("Throat state mass fractions from CEA.");
        writeln("massf=", state6.ceaSavedData.massf);
    }
    //
    // Supersonic expansion.
    auto ar_schedule = new Schedule(xi, ai);
    int n = 20;
    double dx = (xi[$-1] - xi[0])/n;
    foreach (i; 0..n+1) {
        double x = xi[0] + i*dx;
        double a = ar_schedule.interpolate_value(x);
        writeln("x= ", x, " a= ", a);
    }
    return exitFlag;
} // end main
