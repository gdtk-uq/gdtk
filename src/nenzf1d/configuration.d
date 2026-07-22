// configuration.d
// Store the output of a yaml configuration file in a struct of primitive variables.
// This will help with calling nenzf1d from external programs, who may want to modify
// these directly.
//
// Authors: Nick Gibbons, Peter J., Rowan Gollan

module nenzf1d.configuration;

import std.conv;
import std.stdio;
import std.string;

import dyaml;


struct Config{
    string gm1_filename;
    string gm2_filename;
    string reactions_filename;
    string reactions_filename2="";
    string[] species;
    double[] molef;
    double T1=300.0;
    double p1=1.0e5;
    double Vs=0.0;
    double pe=0.0;
    double Te=0.0;
    double meq_throat=1.0;
    double ar=1.0;
    double pp_ps=0.0;
    double C=1.0;
    double[] xi;
    double[] di;
    double x_end;
    double t_final=1.0e-2;  // fixme, this used to default to 2.0 * x_end / v
    double t_inc=1.0e-10;
    double t_inc_factor = 1.0001;
    double t_inc_max= 1.0e-7;
}

// Read the D-YAML Node data and convert it into primitive data types.
// Note that some parameters have defaults, which are simply left in
// place if they are not specified in the input file.
Config process_config(Node configdata) {
    Config config;
    config.gm1_filename = configdata["gas-model-1"].as!string;
    //
    // The second gas model is for the expansion process that has finite-rate 1T chemistry.
    config.gm2_filename = configdata["gas-model-2"].as!string;
    config.reactions_filename = configdata["reactions"].as!string;
    if ("reactions-file2" in configdata) config.reactions_filename2 = configdata["reactions-file2"].as!string;
    //
    foreach(string name; configdata["species"]) { config.species ~= name; }
    foreach(name; config.species) {
        double mf = 0.0;
        if ("molef" in configdata) {
            if (name in configdata["molef"]) mf = configdata["molef"][name].as!double;
        }
        config.molef ~= mf;
    }
    //
    // Initial gas state in shock tube.
    // These initial-state and shock speed data are needed for the shock-tube analysis
    // but not for the case of skipping to the direct setting of the nozzle-supply condition.
    if ("T1" in configdata) config.T1 = configdata["T1"].as!double;
    if ("p1" in configdata) config.p1 = configdata["p1"].as!double;
    if ("Vs" in configdata) config.Vs = configdata["Vs"].as!double;
    //
    // Observed relaxation pressure for reflected-shock, nozzle-supply region.
    // A value of 0.0, or not present, indicates that we should use the ideal shock-reflection pressure.
    if ("pe" in configdata) config.pe = configdata["pe"].as!double;
    // Specified temperature for nozzle-supply region.
    // A value of 0.0, or not present, indicates that we should use do the analysis of the shock-tube processes,
    // but a nonzero value (together with a nonzero pe) skips directly to setting the
    // stagnation condition for the nozzle-supply region.
    if ("Te" in configdata) {
        config.Te = configdata["Te"].as!double;
        if ("pe" !in configdata) throw new Error("Missing stagnation temperature pe. When using Te, this parameter is required.");
    } else {
        if ("T1" !in configdata) throw new Error("Missing shock tube fill temperature T1. Without Te this parameter is required.");
        if ("p1" !in configdata) throw new Error("Missing shock tube fill pressure p1. Without Te this parameter is required.");
        if ("Vs" !in configdata) throw new Error("Missing incident shock velocity Vs. Without Te this parameter is required.");
    }
    //
    // Mach number (in equilibrium gas) at nozzle throat.
    // Nominally, it would be 1.0 for sonic flow.
    // It may be good to expand a little more so that,
    // on changing to the frozen-gas sound speed in the nonequilibrium gas model,
    // the flow remains slightly supersonic.
    if ("meq_throat" in configdata) config.meq_throat = configdata["meq_throat"].as!double;
    //
    // Nozzle-exit area ratio for terminating expansion.
    if ("ar" in configdata) config.ar = configdata["ar"].as!double;
    //
    // Alternatively, we might stop on pPitot/pSupply becoming less than pp_ps.
    // A value of zero, or not present, removes this stopping criterion.
    if ("pp_ps" in configdata) config.pp_ps = configdata["pp_ps"].as!double;
    // pPitot = C * rho*V^^2
    // In a number of sphere simulations, for flows representative of T4 flow conditions,
    // values of pPitot/(rho*v^^2) appeared to be in the range 0.96 to 1.0.
    // TODO: Per discussion with PJ in July 2026, we might want to change the
    // default value here to 0.94, as recommended by Sopek et al.
    // (doi.org/10.1016/j.actaastro.2024.07.008).
    if ("C" in configdata) config.C = configdata["C"].as!double;
    //
    // Nozzle x,diameter schedule. These parameters are required
    foreach(double val; configdata["xi"]) { config.xi ~= val; }
    foreach(double val; configdata["di"]) { config.di ~= val; }

    // Part 2 of config, all optional parameters
    config.x_end = config.xi[$-1];  // Nozzle-exit position for terminating expansion.
    if ("x_end" in configdata)        config.x_end        = configdata["x_end"].as!double;
    if ("t_final" in configdata)      config.t_final      = configdata["t_final"].as!double;
    if ("t_inc" in configdata)        config.t_inc        = configdata["t_inc"].as!double;
    if ("t_inc_factor" in configdata) config.t_inc_factor = configdata["t_inc_factor"].as!double;
    if ("t_inc_max" in configdata)    config.t_inc_max    = configdata["t_inc_max"].as!double;
    return config;
} // end process_config()
