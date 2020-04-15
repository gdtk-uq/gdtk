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

import config;
import simcore;

int main(string[] args)
{
    int exitFlag = 0; // Presume OK in the beginning.

    // We assemble the usage messages as multi-line strings.
    // Be careful when editing them and try to limit the line length
    // to something that is likely to easily fit on a console,
    // say 80 characters.
    string briefUsageMsg = "Usage: l1d4-run... [OPTION]...
Top-level arguments include the following.
Argument:                            Comment:
--------------------------------------------------------------------------------
  --job=<string>                     file names built from this string
  --tindx=<int>                      time index of the starting data (default 0)
  --help                             writes this help message
--------------------------------------------------------------------------------
";
    //
    if ( args.length < 2 ) {
        writeln("Too few arguments.");
        writeln(briefUsageMsg);
        stdout.flush();
        exitFlag = 1;
        return exitFlag;
    }
    //
    string jobName = "";
    int verbosityLevel = 1; // default to having a little information
    int tindxStart = 0;
    bool helpWanted = false;
    try {
        getopt(args,
               "job", &jobName,
               "verbosity", &verbosityLevel,
               "tindx", &tindxStart,
               "help", &helpWanted
               );
    } catch (Exception e) {
        writeln("Problem parsing command-line options.");
        writeln("Arguments not processed: ");
        args = args[1 .. $]; // Dispose of program name in first argument.
        foreach (myarg; args) writeln("    arg: ", myarg);
        writeln(briefUsageMsg);
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
        writeln(briefUsageMsg);
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
        writeln(briefUsageMsg);
        exitFlag = 1;
        return exitFlag;
    }
    L1dConfig.job_name = jobName;
    L1dConfig.verbosity_level = verbosityLevel;

    // Get to work.
    init_simulation(tindxStart);

    if (verbosityLevel > 0) { writeln("Done simulation."); }

    return exitFlag;
} // end main
