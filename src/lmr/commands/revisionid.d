/**
 * Command to display revision id to stdou.
 *
 * Authors: RJG, PJ, KAD, NNG
 * Date: 2022-08-08
 */

module lmr.commands.revisionid;

import std.getopt;
import std.stdio;
import std.string;

import lmr.commands.command;
import lmr.lmrconfig;

Command revisionIdCmd;
string cmdName = "revision-id";

static this()
{
    revisionIdCmd.main = &main_;
    revisionIdCmd.description = "Print version control revision ID for build of Eilmer.";
    revisionIdCmd.shortDescription = "Print version control revision ID.";
    revisionIdCmd.helpMsg =
`lmr revision-id [options]

Print the git revision id of the source code for this build of the
Eilmer simulation program.

If no options are provided, a short form revision ID is given.

option:

 -f, --full
     print the revision ID as its full 40-character SHA-1 hash.
     default: false

`;

}

int main_(string[] args)
{
    bool showFull = false;
    try {
        getopt(args,
               config.bundling,
               "f|full", &showFull);
    } catch (Exception e) {
        writefln("Eilmer %s program quitting.", cmdName);
        writeln("There is something wrong with the command-line arguments/options.");
        writeln(e.msg);
        return 1;
    }

    if (showFull) {
        writeln(lmrCfg.fullRevisionId);
        return 0;
    }

    writeln(lmrCfg.revisionId);
    return 0;
}


