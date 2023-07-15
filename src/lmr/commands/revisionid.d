/**
 * Command to display revision id to stdou.
 *
 * Authors: RJG, PJ, KAD, NNG
 * Date: 2022-08-08
 */

module revisionid;

import std.stdio;
import std.getopt;
import std.string;

import lmrconfig;
import command;

Command revisionIdCmd;

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

void main_(string[] args)
{
    bool showFull = false;

    getopt(args,
           config.bundling,
           "f|full", &showFull);

    if (showFull) {
        writeln(lmrCfg.fullRevisionId);
        return;
    }

    writeln(lmrCfg.revisionId);
    return;
}


