// lorikeet_main.d
// Lorikeet, a calculator for transient 2D compressible flow.
//
// PA Jacobs
// 2022-12-12 Adapt the Puffin and Chicken codes into Lorikeet.
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
import puffin.fluidblock;
import puffin.transient_calc;

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
    string usageMsg = "Usage: lrtk-run ... [OPTION]...
Top-level arguments include the following.
Argument:                            Comment:
--------------------------------------------------------------------------------
Actions:
  Running a transient-flow simulation is the only real action.
  --help                             writes this help message

Parameters:
  --job=<string>                     file names built from this string
  --tindx=<int>                      starting index for flow-field data (default 0)
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
    int tindx = 0;
    try {
        getopt(args,
               "job", &jobName,
               "verbosity", &verbosityLevel,
               "tindx", &tindx,
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
        writeln("Lorikeet 2D compressible transient-flow simulator.");
        writeln("Revision: PUT_REVISION_STRING_HERE");
        writeln("Compiler-name: PUT_COMPILER_NAME_HERE");
        //
        write("Build-flavour: ");
        version(flavour_debug) { writeln("debug"); }
        version(flavour_profile) { writeln("profile"); }
        version(flavour_fast) { writeln("fast"); }
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
    if (verbosityLevel > 0) { writeln("Do a simulation."); }
    init_simulation(tindx);
    do_time_integration();
    if (verbosityLevel > 0) { writeln("Done."); }
    //
    return exitFlag;
} // end main
