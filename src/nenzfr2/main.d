// main.d for nenzfr2
// PJ 2020-09-26 Initial code built from Python prototype.

import std.stdio;
import std.array;
import std.string;
import std.conv;
import dyaml;

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
    double pe = to!double(config["pe"].as!string);
    double ar = to!double(config["ar"].as!string);
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
    return 0;
} // end main
