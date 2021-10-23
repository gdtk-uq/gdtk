// configuration.d
// Store the output of a yaml configuration file in a struct of primitive variables.
// This will help with calling nenzf1d from external programs, who may want to modify
// these directly.
//
// Authors: Nick Gibbons, Peter J., Rowan Gollan

module configuration;

import std.stdio;
import std.string;
import std.conv;
import dyaml;


struct Config{
    string gm1_filename;
    string gm2_filename;
    string reactions_filename;
    string reactions_filename2="";
    string[] species;
    double[] molef;
    double T1;
    double p1;
    double Vs;
    double pe=0.0;
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
    try {
        config.reactions_filename2 = configdata["reactions-file2"].as!string;
    } catch (YAMLException e) {
        // We cannot set the second reactions file so assume 1T chemistry.
    }
    //
    foreach(string name; configdata["species"]) { config.species ~= name; }
    foreach(name; config.species) {
        double mf = 0.0;
        try { mf = to!double(configdata["molef"][name].as!string); } catch (YAMLException e) {}
        config.molef ~= mf;
    }
    //
    // Initial gas state in shock tube.
    config.T1 = to!double(configdata["T1"].as!string);
    config.p1 = to!double(configdata["p1"].as!string);
    config.Vs = to!double(configdata["Vs"].as!string);
    // Observed relaxation pressure for reflected-shock, nozzle-supply region.
    // A value of 0.0 indicates that we should use the ideal shock-reflection pressure.
    try { config.pe = to!double(configdata["pe"].as!string); } catch (YAMLException e) {}
    // Mach number (in equilibrium gas) at nozzle throat.
    // Nominally, it would be 1.0 for sonic flow.
    // It may be good to expand a little more so that,
    // on changing to the frozen-gas sound speed in the nonequilibrium gas model,
    // the flow remains slightly supersonic.
    try { config.meq_throat = to!double(configdata["meq_throat"].as!string); } catch (YAMLException e) {}
    // Nozzle-exit area ratio for terminating expansion.
    try { config.ar = to!double(configdata["ar"].as!string); } catch (YAMLException e) {}
    // Alternatively, we might stop on pPitot/pSupply becoming less than pp_ps.
    try { config.pp_ps = to!double(configdata["pp_ps"].as!string); } catch (YAMLException e) {}
    // pPitot = C * rho*V^^2
    // A value of C=1.0 is a good default.
    // In a number of sphere simulations, for flows representative of T4 flow conditions,
    // values of pPitot/(rho*v^^2) appeared to be in the range 0.96 to 1.0.
    try { config.C = to!double(configdata["C"].as!string); } catch (YAMLException e) {}
    // Nozzle x,diameter schedule.
    foreach(string val; configdata["xi"]) { config.xi ~= to!double(val); }
    foreach(string val; configdata["di"]) { config.di ~= to!double(val); }

    // Part 2 of config
    config.x_end = config.xi[$-1];  // Nozzle-exit position for terminating expansion.
    try { config.x_end = to!double(configdata["x_end"].as!string); } catch (YAMLException e) {}
    try { config.t_final = to!double(configdata["t_final"].as!string); } catch (YAMLException e) {}
    try { config.t_inc = to!double(configdata["t_inc"].as!string); } catch (YAMLException e) {}
    try { config.t_inc_factor = to!double(configdata["t_inc_factor"].as!string); } catch (YAMLException e) {}
    try { config.t_inc_max = to!double(configdata["t_inc_max"].as!string); } catch (YAMLException e) {}
    return config;
} // end process_config()
