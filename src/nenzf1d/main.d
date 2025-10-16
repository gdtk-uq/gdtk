// main.d: Top level function for directly called nenzf1d,
// a quasi-one-dimensional calculator for shock tunnel conditions.

import std.file;
import std.getopt;
import std.math;
import std.stdio;
import std.string;

import dyaml;

import nenzf1d.configuration;
import nenzf1d.shock_tube_nozzle;

import util.buildinfo : buildCfg;

int main(string[] args)
{
    int exitFlag = 0; // Presume OK in the beginning.
    version(enable_fp_exceptions) {
        FloatingPointControl fpctrl;
        // Enable hardware exceptions for division by zero, overflow to infinity,
        // invalid operations, and uninitialized floating-point variables.
        // Copied from https://dlang.org/library/std/math/floating_point_control.html
        fpctrl.enableExceptions(FloatingPointControl.severeExceptions);
    }
    //
    // Be careful with the usageMsg string; it has embedded newline characters.
    string usageMsg = "Usage: nenzf1d <input-file>
Options:
   --verbosity=<int>   0 == very terse output
                       1 == key results printed (default)
                       2 == echo input as well as printing more detailed results
                       3 == debug printing as well
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
        writeln("NENZF1D: 1D analysis of shock-tunnel with nonequilibrium nozzle flow.");
        writeln("  Revision-id: ", buildCfg.revisionId);
        writeln("  Revision-date: ", buildCfg.revisionDate);
        writeln("  Compiler-name: ", buildCfg.compilerName);
        writeln("  Build-date: ", buildCfg.buildDate);
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
    Node configdata = dyaml.Loader.fromFile(inputFile).load();
    if (verbosityLevel >= 1) {
        writeln("  "~configdata["title"].as!string);
    }
    // Convert raw configdata into primitive data types
    Config config = process_config(configdata);
    //
    // Run nenzf1d simulation and output results
    try{
        auto result = analyse(verbosityLevel, config);
        writeln("Exit condition:");
        writefln("  x           %g m", result.x);
        writefln("  area-ratio  %g", result.area_ratio);
        writefln("  velocity    %g m/s", result.velocity);
        writefln("  Mach        %g", result.Mach_number);
        writefln("  p_pitot     %g kPa (C.rho.V^2)", result.p_pitot/1000);
        if (result.rayleigh_pitot > 0.0) {
            writefln("  p_pitot     %g kPa (Rayleigh-Pitot, frozen kinetics)",
                    result.rayleigh_pitot/1000);
        }
        writefln("  pressure    %g kPa", result.pressure/1000.0);
        writefln("  density     %g kg/m^3", result.density);
        writefln("  temperature %g K", result.temperature);
        foreach (i; 0 .. result.T_modes.length) {
            string label = format("T_modes[%d]", i);
            writefln("%s%-12s%g K", "  ", label, result.T_modes[i]);
        }
        foreach (i, name; config.species) {
            string label = format("massf[%s]", name);
            writefln("%s%-12s%g", "  ", label, result.massf[i]);
        }
        writefln("  viscosity   %g Pa.s", result.viscosity);
        //
        writeln("Expansion error-indicators:");
        writefln("  relerr-mass %g", result.massflux_rel_err);
        writefln("  relerr-H    %g", result.enthalpy_rel_err);
        if (result.rayleigh_pitot > 0.0) {
            writefln("  relerr-pitot %g", result.pitot_rel_err);
        }
    } catch (Exception e) {
        writeln("Exception thrown in nenzf1d.run!");
        writefln("  Exception message: %s", e.msg);
        exitFlag=3;
    }
    return exitFlag;
} // end main()
