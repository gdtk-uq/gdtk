import std.getopt;
import std.stdio;
import std.algorithm;

import lmrconfig;
import command;
import checkjacobian;
import computenorms;
import probeflow;
import limiter2vtk;
import prepgrids;
import prepsim;
import prepmappedcells;
import runsim;
import snapshot2vtk;
import structured2unstructured;
import revisionid;

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
`lmr help [-a] <TOPIC>

Show help for a given Eilmer command or topic.

   With no arguments, print a list of commonly used commands.

   With option "-a", print all commands with short descriptions.

   Given a command name or topic, print specific help.
`;

    // Add commands here.
    commands["help"] = helpCmd;
    // Try to add commands in alphabetical order from here down.
    // 1. Add user commands
    commands["compute-norms"] = compNormsCmd;
    commands["probe-flow"] = probeFlowCmd;
    commands["limiter2vtk"] = limiter2vtkCmd;
    commands["prep-grids"] = prepGridCmd;
    commands["prep-grid"] = commands["prep-grids"]; // alias for prep-grids
    commands["prep-sim"] = prepSimCmd;
    commands["prep-flow"] = commands["prep-sim"]; // alias for prep-sim
    commands["prep-mapped-cells"] = prepMappedCellsCmd;
    commands["revision-id"] = revisionIdCmd;
    commands["run"] = runCmd;
    commands["snapshot2vtk"] = snapshot2vtkCmd;
    commands["structured2unstructured"] = structured2unstructuredCmd;
    // add alias for structured2unstructured
    commands["sgrid2ugrid"] = commands["structured2unstructured"];
    // 2. Add dev/diag commands
    commands["check-jacobian"] = checkJacCmd;
}

void main(string[] args)
{
    bool helpWanted = false;
    bool versionWanted = false;
    bool versionLongWanted = false;
    NumberType numberType;
    getopt(args,
           std.getopt.config.stopOnFirstNonOption,
           "h|help", &helpWanted,
           "v|version", &versionWanted,
           "version-long", &versionLongWanted,
           "number-type", &numberType
    );

    if (args.length < 2) {
        // Nothing asked for. Print help and exit.
        printHelp(args);
        return;
    }

    if (versionLongWanted) {
        printVersion(false);
        return;
    }
    else if (versionWanted) {
        printVersion();
        return;
    }
    if (helpWanted) printHelp(args);

    // Special cases for version options written as commands.
    if (args[1] == "version-long") {
        printVersion(false);
        return;
    }
    else if (args[1] == "version") {
        printVersion();
        return;
    }

    if (args[1].startsWith("--number-type")) {
        args.remove(1);
    }

    auto cmd = args[1];

    if (cmd == "run") {
        // We need to treat this one specially because of how delegation
        // is made based on number type (real or complex)
        return runsim.delegateAndExecute(args, numberType);
    }

    if (cmd in commands) {
        return (*commands[cmd].main)(args);
    }
    // If we've made it here, then we've chosen a bad command.
    writefln("lmr: '%s' is not an lmr command. See 'lmr help'.", cmd);
    return;
}

void listHelpForAllCommands()
{
    writeln("See 'lmr help <command>' to read about a specific subcommand.");
    writeln("");
    auto cmds = commands.keys;
    cmds.sort;
    // First print regular user commands
    writeln("Available commands");
    foreach (cmd; cmds) {
        if (commands[cmd].type == LmrCmdType.user)
            writefln("   %-24s %s", cmd, commands[cmd].shortDescription);
    }
    writeln("");
    writeln("Developer/diagnostics commands");
    foreach (cmd; cmds) {
        if (commands[cmd].type == LmrCmdType.dev)
            writefln("   %-24s %s", cmd, commands[cmd].shortDescription);
    }
    writeln("");
    writeln("Meta commands");
    writefln("   %-24s Print condensed version information about lmr program.", "version");
    writefln("   %-24s Print full version information about lmr program.", "version-long");
}

void printHelp(string[] args)
{
    if (args.length >= 3) {
        auto arg = args[2];
        if (arg == "-a" || arg == "--all") {
            // list all commands with short description.
            listHelpForAllCommands();
            return;
        }
        // Next look to see if a command has been supplied.
        if (arg in commands) {
            writeln(commands[arg].helpMsg);
            return;
        }
        // If we've made it here, then we've chosen a bad command.
        writefln("lmr: '%s' is not an lmr command. See 'lmr help'.", arg);
    }
    // else just print general help
    string generalHelp =
`usage: lmr [-h | --help] [help -a]
           [-v | --version] [--version-long]
           [--number-type=real_values|complex_values]
            <command> [<args>]

== Eilmer simulation program ==

List of commonly used commands:

at preparation stage
   prep-grids      build grids for simulation
   prep-sim        build initial flow field for simulation

at simulation stage
   run             run flow solver

at post-processing stage
   snapshot2vtk    convert a snapshot to VTK format for visualisation
   probe-flow      reports the flow-field data at specified location(s)

== Notes ==
--number-type option, if used, must appear before "run" command.
It is advanced usage to control delegation of "run" to lmr-run
(which uses real numbers) or lmrZ-run (which uses complex numbers).

`;
    write(generalHelp);
    return;

}


void printVersion(bool shortVersion=true)
{
    if (GlobalConfig.is_master_task) {
        writeln("Eilmer 4.0 compressible-flow simulation code.");
        writeln("Revision-id: ", lmrCfg.revisionId);
        writeln("Revision-date: ", lmrCfg.revisionDate);
        writeln("Compiler-name: ", lmrCfg.compilerName);
        writeln("Build-date: ", lmrCfg.buildDate);

        if (shortVersion) return;

        write("Build-flavour: ");
        version(flavour_debug) { writeln("debug"); }
        version(flavour_profile) { writeln("profile"); }
        version(flavour_fast) { writeln("fast"); }
        write("Profiling: ");
        version(diagnostics) { writeln("included"); } else { writeln("omitted"); }
        //
        writeln("Capabilities:");
        version(multi_species_gas) {
            writeln("   multi-species-gas");
        } else {
            writeln("   single-species-gas");
        }
        version(multi_T_gas) {
            writeln("   multi-temperature-gas");
        } else {
            writeln("   single-temperature-gas");
        }
        version(MHD) {
            writeln("   MHD");
        } else {
            writeln("   no-MHD");
        }
        version(turbulence) {
            writeln("   turbulence");
        } else {
            writeln("   no-turbulence-modelling");
        }
    }
}


