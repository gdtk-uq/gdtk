import std.getopt;
import std.stdio;
import std.algorithm;

import lmrconfig;
import command;
import computenorms;
import customscript;
import probeflow;
import sliceflow;
import extractline;
import limiter2vtk;
import residual2vtk;
import lmr.commands.listspecies;
import lmr.commands.plotdiagnostics;
import lmr.commands.prepenergyexchange;
import lmr.commands.prepgas;
import lmr.commands.prepreactions;
import lmr.commands.slicesolid;
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
    commands["custom-script"] = customScriptCmd;
    commands["extract-line"] = extractLineCmd;
    commands["limiter2vtk"] = limiter2vtkCmd;
    commands["list-species"] = listSpeciesCmd;
    // alias for list-species, provided for consistency with prep-gas tool
    commands["list-available-species"] = commands["list-species"];
    commands["residual2vtk"] = residual2vtkCmd;
    commands["plot-diagnostics"] = plotDiagnosticsCmd;
    commands["prep-energy-exchange"] = prepExchCmd;
    // alias for 'prep-energy-exchange' provided for consistency with prep-kinetics tool
    commands["prep-kinetics"] = commands["prep-energy-exchange"];
    commands["prep-gas"] = prepGasCmd;
    commands["prep-grids"] = prepGridCmd;
    commands["prep-grid"] = commands["prep-grids"]; // alias for prep-grids
    commands["prep-reactions"] = prepReacCmd;
    // alias for 'prep-reactions' provided for consistency with prep-chem tool
    commands["prep-chem"] = commands["prep-reactions"];
    commands["prep-sim"] = prepSimCmd;
    commands["prep-flow"] = commands["prep-sim"]; // alias for prep-sim
    commands["prep-mapped-cells"] = prepMappedCellsCmd;
    commands["probe-flow"] = probeFlowCmd;
    commands["revision-id"] = revisionIdCmd;
    commands["run"] = runCmd;
    commands["slice-flow"] = sliceFlowCmd;
    // alias for slice-flow because fluid comes naturally when also working with solids
    commands["slice-fluid"] = commands["slice-flow"];
    commands["slice-solid"] = sliceSolidCmd;
    commands["snapshot2vtk"] = snapshot2vtkCmd;
    commands["structured2unstructured"] = structured2unstructuredCmd;
    // add alias for structured2unstructured
    commands["sgrid2ugrid"] = commands["structured2unstructured"];
    // 2. Add dev/diag commands
}

int main(string[] args)
{
    bool helpWanted = false;
    bool versionWanted = false;
    bool versionLongWanted = false;
    NumberType numberType;
    try {
        getopt(args,
               std.getopt.config.stopOnFirstNonOption,
               "h|help", &helpWanted,
               "v|version", &versionWanted,
               "version-long", &versionLongWanted,
               "number-type", &numberType
               );
    } catch (Exception e) {
        writeln("Eilmer top-level program quitting.");
        writeln("There is something wrong with the command-line arguments/options.");
        writeln(e.msg);
        writeln("");
        printHelp([""]);
        return 1;
    }

    if (versionLongWanted) {
        printVersion(false);
        if (helpWanted || canFind(args, "help")) {
            writeln("");
            printHelp([""]); // general help
        }
        return 0;
    }
    else if (versionWanted) {
        printVersion();
        if (helpWanted || canFind(args, "help")) {
            writeln("");
            printHelp([""]); // general help
        }
        return 0;
    }

    if (args.length < 2) {
        writeln("Eilmer top-level program quitting.");
        writeln("No subcommand supplied as the first command-line argument.");
        writeln("");
        printHelp([""]); // general help
        return 1;
    }

    // Special cases for version options written as commands.
    if (args[1] == "version-long") {
        printVersion(false);
        return 0;
    }
    else if (args[1] == "version") {
        printVersion();
        return 0;
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
    return 1;
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

int printHelp(string[] args)
{
    if (args.length >= 3) {
        auto arg = args[2];
        if (arg == "-a" || arg == "--all") {
            // list all commands with short description.
            listHelpForAllCommands();
            return 0;
        }
        // Next look to see if a command has been supplied.
        if (arg in commands) {
            writeln(commands[arg].helpMsg);
            return 0;
        }
        // If we've made it here, then we've chosen a bad command.
        writefln("lmr: '%s' is not an lmr command. See 'lmr help'.", arg);
    }
    // else just print general help
    string generalHelp =
`== Eilmer simulation program ==

Usage: lmr [-h | --help] [help -a]
           [-v | --version] [--version-long]
           [--number-type=real_values|complex_values]
            <command> [<args>]

Examples:
lmr help             ==> prints this general help
lmr help -a          ==> lists all commands
lmr help <command>   ==> prints help for <command>

== Commonly used commands ==

at preparation stage
   prep-grids      build grids for simulation
   prep-sim        build initial flow field for simulation

at simulation stage
   run             run flow solver

at post-processing stage
   snapshot2vtk    convert a snapshot to VTK format for visualisation
   probe-flow      reports the flow-field data at specified location(s)
   slice-flow      reports the flow-field data along slices, in index directions
   extract-line    reports the flow-field data along lines in 3D space

== Notes ==
--number-type option, if used, must appear before "run" command.
It is advanced usage to control delegation of "run" to lmr-run
(which uses real numbers) or lmrZ-run (which uses complex numbers).

`;
    write(generalHelp);
    return 0;
}


void printVersion(bool shortVersion=true)
{
    if (GlobalConfig.is_master_task) {
        writeln("Eilmer 5.0 compressible-flow simulation code.");
        writeln("Revision-id: ", lmrCfg.revisionId);
        writeln("Revision-date: ", lmrCfg.revisionDate);
        writeln("Compiler-name: ", lmrCfg.compilerName);
        writeln("Build-date: ", lmrCfg.buildDate);
        write("Build-flavour: ");
        version(flavour_debug) { writeln("debug"); }
        version(flavour_profile) { writeln("profile"); }
        version(flavour_fast) { writeln("fast"); }
        write("Profiling: ");
        version(diagnostics) { writeln("included"); } else { writeln("omitted"); }

        if (shortVersion) return;

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


