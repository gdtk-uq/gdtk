import std.getopt;
import std.stdio;

import command;
import prepgrids;
import prepflow;
import runsteady;
import snapshot2vtk;
// Eilmer4 imports
import globalconfig : GlobalConfig;

Command[string] commands;
Command helpCmd;

static this()
{
    // Initialise helpCmd in this module so that it has access to
    // all commands.
    helpCmd.main = &printHelp;
    helpCmd.description = "Display help information for a topic/command or a general overview.";
    helpCmd.shortDescription = "Display help about using Eilmer.";
    helpCmd.helpMsg =
`lmr help [TOPIC]

Show help for a given Eilmer command or topic.

   With no arguments, print a list of commonly used commands.

   Given a command name or topic, print specific help.
`;

    // Add commands here.
    commands["help"] = helpCmd;

    commands["prep-grids"] = prepGridCmd;
    {
        // add alias
        commands["prep-grid"] = commands["prep-grids"];
    }
    commands["prep-flow"] = prepFlowCmd;
    commands["run-steady"] = runSteadyCmd;
    commands["snapshot2vtk"] = snapshot2vtkCmd;
}

void main(string[] args)
{
    bool helpWanted = false;
    bool versionWanted = false;
    bool versionLongWanted = false;
    getopt(args,
           std.getopt.config.stopOnFirstNonOption,
           "help|h", &helpWanted,
           "version|v", &versionWanted,
           "version-full", &versionLongWanted, 
    );

    if (versionLongWanted) {
        printVersion(false);
    }
    else if (versionWanted) {
        printVersion();
    }
    if (helpWanted) printHelp(args);

    if (args.length < 2) {
        // Nothing asked for. Print help and exit.
        printHelp(args);
        return;
    }
    auto cmd = args[1];

    if (cmd in commands) {
        return (*commands[cmd].main)(args);
    }
    // If we've made it here, then we've chosen a bad command.
    writefln("lmr: '%s' is not an lmr command. See 'lmr help'.", cmd);
    return;
}

void printHelp(string[] args)
{
    if (args.length >= 3) {
        auto cmd = args[2];
        if (cmd in commands) {
            writeln(commands[cmd].helpMsg);
            return;
        }
        // If we've made it here, then we've chosen a bad command.
        writefln("lmr: '%s' is not an lmr command. See 'lmr help'.", cmd);
    }
    // else just print general help
    string generalHelp =
`usage: lmr [-v |--version] [-h | --help] [--version-long]
            <command> [<args>]

== Eilmer simulation program ==

List of commonly used commands:

at preparation stage
   prep-grids      build grids for simulation
   prep-flow       build initial flow field for simulation

at simulation stage
   run-steady      run steady-state solver
   ...

at post-processing stage
   snapshot2vtk    convert a snapshot to VTK format for visualisation

`;
    write(generalHelp);
    return;
    
}


void printVersion(bool shortVersion=true)
{
    if (GlobalConfig.is_master_task) {
        writeln("Eilmer 4.0 compressible-flow simulation code.");
        writeln("Revision-id: PUT_REVISION_STRING_HERE");
        writeln("Revision-date: PUT_REVISION_DATE_HERE");
        writeln("Compiler-name: PUT_COMPILER_NAME_HERE");
        writeln("Build-date: PUT_BUILD_DATE_HERE");

        if (shortVersion) return;

        write("Build-flavour: ");
        version(flavour_debug) { writeln("debug"); }
        version(flavour_profile) { writeln("profile"); }
        version(flavour_fast) { writeln("fast"); }
        write("Profiling: ");
        version(diagnostics) { writeln("included"); } else { writeln("omitted"); }
        //
        write("Capabilities:");
        version(multi_species_gas) {
            write(" multi-species-gas");
        } else {
            write(" single-species-gas");
        }
        version(multi_T_gas) {
            write(" multi-temperature-gas");
        } else {
            write(" single-temperature-gas");
        }
        version(MHD) {
            write(" MHD");
        } else {
            write(" no-MHD");
        }
        version(turbulence) {
            write(" turbulence");
        } else {
            write(" no-turbulence-modelling");
        }
    }
}


