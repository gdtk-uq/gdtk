/**
 * Command to restore a default commands file (related to steady mode run-time).
 *
 * Author: RJG
 * Date: 2025-04-05
 */

module lmr.commands.restore_default_run_time_commands;

import std.file : copy, FileException;
import std.stdio : writeln, writefln;
import std.string : format;
import std.getopt;
import std.process : environment;

import lmr.commands.command;
import lmr.lmrconfig : lmrCfg;
import lmr.lmrerrors: LmrError, lmrErrorExit;

Command restoreDfltCmd;

string cmdName = "restore-default-run-time-commands";

static this()
{
    restoreDfltCmd.main = &main_;
    restoreDfltCmd.description = "Restore defaults for steady-mode run-time commands.";
    restoreDfltCmd.shortDescription = restoreDfltCmd.description;
    restoreDfltCmd.helpMsg = format(
`lmr %s [option]

Restore the default file for steady solver mode run-time commands.

This might be useful when preparing for a restarted solution from
a snapshot > 0, particularly if you have hand-edited the commands
file and want a clean restart with no actions.

option:

 -v, --verbose
     Explicitly print what is being copied.
    
`, cmdName);
}

int main_(string[] args)
{
    int verbosity = 0;
    try {
        getopt(args,
               config.bundling,
               "v|verbose+", &verbosity
               );
    }
    catch (Exception e) {
        string exitMsg = format("Eilmer %s command quitting.", cmdName);
        exitMsg ~= ("There is something wrong with the command-line arguments/options.");
        exitMsg ~= format("%s", e.msg);
        lmrErrorExit(LmrError.commandLineError, exitMsg);
    }

    string fromFile = format("%s/share/default-commands-to-steady-mode", environment.get("DGD"));
    string toFile = lmrCfg.nkCmdsFile;

    if (verbosity > 0) {
        writeln("Copying defaults file.");
        writefln("  FROM: %s", fromFile);
        writefln("  TO:   %s", toFile);
    }

    try {
        copy(fromFile, toFile);
    }
    catch (FileException e) {
        string exitMsg = format("ERROR: lmr %s did not complete succesfully.\n", cmdName);
        exitMsg ~= "Are you calling this command AFTER you have prepared a simulation?\n";
        exitMsg ~= "This command relies on lmrsim/ directory being present.\n";
        exitMsg ~= format("%s", e.msg);
        lmrErrorExit(LmrError.inputOutput, exitMsg);
    }

    return 0;
}
