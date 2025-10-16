// puffin_main.d
// Puffin, a calculator for steady-state 2D supersonic flow.
//
// PA Jacobs
// 2022-01-21 Start code.
//

import std.algorithm;
import std.conv;
import std.file;
import std.getopt;
import std.math;
import std.parallelism;
import std.path;
import std.stdio;
import std.string;

import puffin.config;
import puffin.marching_calc;
import puffin.streamtube;

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
    // We assemble the usage messages as multi-line strings.
    // Be careful when editing them and try to limit the line length
    // to something that is likely to easily fit on a console,
    // say 80 characters.
    string usageMsg = "Usage: puffin ... [OPTION]...
Top-level arguments include the following.
Argument:                            Comment:
--------------------------------------------------------------------------------
Actions:
  Running a space-marching calculation is the only real action.
  --help                             writes this help message

Parameters:
  --job=<string>                     file names built from this string
  --max-cpus=<int>                   to process blocks in parallel, set to more than 1
  --verbosity=<int>                  level of commentary as the work is done
                                       0=very little written to console
                                       1=major steps commentary (default)
                                       2=minor steps commentary
--------------------------------------------------------------------------------
";
    //
    if ( args.length < 1 ) {
        writeln("Too few arguments.");
        writeln(usageMsg);
        stdout.flush();
        exitFlag = 1;
        return exitFlag;
    }
    //
    string jobName = "";
    int maxCPUs = 1; // default to single-threaded simulation
    int verbosityLevel = 1; // default to commenting on major steps
    bool helpWanted = false;
    try {
        getopt(args,
               "job", &jobName,
               "verbosity", &verbosityLevel,
               "max-cpus", &maxCPUs,
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
    if (verbosityLevel > 0) {
        writeln("Puffin 2D supersonic steady-flow calculator.");
        writeln("Revision: ", buildCfg.revisionId);
        writeln("Compiler-name: ", buildCfg.compilerName);
        writeln("Build-flavour: ", buildCfg.buildFlavour);
    }
    if (helpWanted) {
        writeln(usageMsg);
        stdout.flush();
        exitFlag = 0;
        return exitFlag;
    }
    //
    if (jobName.length > 0) {
        // Clean up the jobName, by removing any extension or path details, if necessary.
        string ext = extension(jobName);
        if (!ext.empty && ext != ".py") {
            writeln("If you are going to supply an extension for your job name, it needs to be \".py\"");
            exitFlag = 1;
            return exitFlag;
        }
        string dir = dirName(jobName);
        if (dir != ".") {
            writeln("You are expected to start with your job script in your working directory.");
            exitFlag = 1;
            return exitFlag;
        }
        string bn = baseName(jobName);
        if (ext.length > 0) {
            jobName = bn.replace(ext, "");
        } else {
            jobName = bn;
        }
    }
    if (jobName.length == 0) {
        writeln("Need to specify a job name.");
        writeln(usageMsg);
        exitFlag = 1;
        return exitFlag;
    }
    Config.job_name = jobName;
    Config.verbosity_level = verbosityLevel;
    Config.maxCPUs = min(max(maxCPUs, 1), totalCPUs);
    //
    if (verbosityLevel > 0) { writeln("Do a Calculation."); }
    init_calculation();
    do_space_marching_calculation();
    if (verbosityLevel > 0) { writeln("Done."); }
    //
    return exitFlag;
} // end main
