// main.d of
// Puffin, a calculator for steady-state 2D supersonic flow.
//
// PA Jacobs
// 2022-01-21 Start code.
//

import std.stdio;
import std.string;
import std.file;
import std.path;
import std.getopt;
import std.conv;

import config;
import marching_calc;
import streamtube;

int main(string[] args)
{
    int exitFlag = 0; // Presume OK in the beginning.

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
    int verbosityLevel = 1; // default to commenting on major steps
    bool helpWanted = false;
    try {
        getopt(args,
               "job", &jobName,
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
        writeln("Puffin 2D supersonic flow calculator.");
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
    //
    if (verbosityLevel > 0) { writeln("Do a Calculation."); }
    init_calculation();
    do_space_marching_calculation();
    if (verbosityLevel > 0) { writeln("Done."); }
    //
    return exitFlag;
} // end main
