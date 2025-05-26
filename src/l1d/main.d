// main.d for the Lagrangian 1D Gas Dynamics, also known as L1d4.
//
// PA Jacobs
// 2020-04-05 Minimal code.
// The intention today is that we build just enough to run
// the Sod shock tube and David Gildfind's expanding test gas.
//

import std.stdio;
import std.string;
import std.file;
import std.path;
import std.getopt;
import std.conv;
import std.math;

import config;
import simcore;
import postprocess;

int main(string[] args)
{
    int exitFlag = 0; // Presume OK in the beginning.
    version(enable_fp_exceptions) {
        FloatingPointControl fpctrl;
        // Enable hardware exceptions for division by zero, overflow to infinity,
        // invalid operations, and uninitialized floating-point variables.
        // See https://dlang.org/phobos/std_math_hardware.html#.FloatingPointControl
        // fpctrl.enableExceptions(FloatingPointControl.severeExceptions);
        // More fine-grained control, as explored in Eilmer.
        // fpctrl.enableExceptions(FloatingPointControl.invalidException);
        fpctrl.enableExceptions(FloatingPointControl.divByZeroException);
        fpctrl.enableExceptions(FloatingPointControl.overflowException);
    }

    // We assemble the usage messages as multi-line strings.
    // Be careful when editing them and try to limit the line length
    // to something that is likely to easily fit on a console,
    // say 80 characters.
    string usageMsg = "Usage: l1d4 ... [OPTION]...
Top-level arguments include the following.
Argument:                            Comment:
--------------------------------------------------------------------------------
Actions:
  --run-simulation                   run the simulation
  --time-slice                       extract a single-time slice for gas slugs
  --piston-history                   assemble the history dataset for a piston
  --cell-history                     Assemble the history dataset for a cell
  --xt-data                          generate an xt-dataset for a flow variable (GNUPlot)
  --xt-data-vtk                      generate an xt-dataset for all flow variables (VTK)
  --xt-data-json                     generate an xt-dataset for all flow variables (JSON)
  --trim-solution-files              remove solution file content after tindx-end
  --help                             writes this help message

Parameters:
  --job=<string>                     file names built from this string
  --tindx=<int>                      time index of the starting data (default 0)
  --tindx-end=<int>                  time index of the final data for xt-dataset
                                       (default 9999)
  --pindx=<int>                      piston index for history generation
  --x0=<double>                      cell initial position for history generation
  --var-name=<name>                  flow variable name
  --log10                            take logarithm of flow variable
                                       Useful to get a better range for pressure.
  --millisec                         write time values in milliseconds

  --verbosity=<int>                  level of commentary as the work is done
                                       0=very little written to console
                                       1=major steps commentary (default)
                                       2=minor steps commentary
--------------------------------------------------------------------------------
";
    //
    if ( args.length < 2 ) {
        writeln("Too few arguments.");
        writeln(usageMsg);
        stdout.flush();
        exitFlag = 1;
        return exitFlag;
    }
    //
    string jobName = "";
    bool runSimulation = false;
    bool timeSlice = false;
    bool pistonHistory = false;
    bool cellHistory = false;
    bool xtData = false;
    bool xtDataVTK = false;
    bool xtDataJSON = false;
    bool trimSolutionFiles = false;
    int tindx = 0;
    int tindxEnd = 9999;
    int pistonIndx = 0;
    double x0 = 0.0;
    string varName = "";
    bool takeLog = false;
    bool milliSec = false;
    int verbosityLevel = 1; // default to commenting on major steps
    bool helpWanted = false;
    try {
        getopt(args,
               "job", &jobName,
               "run-simulation", &runSimulation,
               "time-slice", &timeSlice,
               "piston-history", &pistonHistory,
               "cell-history", &cellHistory,
               "xt-data", &xtData,
               "xt-data-vtk", &xtDataVTK,
               "xt-data-json", &xtDataJSON,
               "trim-solution-files", &trimSolutionFiles,
               "tindx", &tindx,
               "tindx-end", &tindxEnd,
               "pindx", &pistonIndx,
               "x0", &x0,
               "var-name", &varName,
               "log10", &takeLog,
               "millisec", &milliSec,
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
    if (verbosityLevel > 0) {
        writeln("L1d 4.0 compressible-flow 1D simulation code.");
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
    L1dConfig.job_name = jobName;
    L1dConfig.verbosity_level = verbosityLevel;

    // Get to work to do one task...
    if (runSimulation) {
        writeln("Run a simulation.");
        init_simulation(tindx);
        integrate_in_time();
    } else if (timeSlice) {
        extract_time_slice(tindx);
    } else if (pistonHistory) {
        assemble_piston_history(pistonIndx);
    } else if (cellHistory) {
        generate_cell_history(x0, tindx, tindxEnd, milliSec);
    } else if (xtData) {
        generate_xt_dataset_gnuplot(varName, tindx, tindxEnd, takeLog);
    } else if (xtDataVTK) {
        generate_xt_dataset(tindx, tindxEnd, milliSec, "VTK");
    } else if (xtDataJSON) {
        generate_xt_dataset(tindx, tindxEnd, milliSec, "JSON");
    } else if (trimSolutionFiles) {
        trim_solution_files(tindxEnd);
    } else {
        writeln("You did not ask for anything to be done.");
    }
    if (verbosityLevel > 0) { writeln("Done."); }
    return exitFlag;
} // end main
