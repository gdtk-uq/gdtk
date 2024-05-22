/**
 * Module for slicing solid-domain snapshots, picking out specific data.
 *
 * Authors: RJG, PJ, KAD, NNG
 * Date: 2024-05-22, adapted from slideflow.d
 */

module lmr.commands.slicesolid;

import std.stdio;
import std.file;
import std.string;
import std.regex : regex, replaceAll;
import std.array;
import std.format : format;
import std.getopt;
import std.conv : to;
import std.algorithm;

import globalconfig;
import fileutil;
import solidsolution;
import lmrconfig : lmrCfg;
import init : initConfiguration;
import cmdhelper;

import command;

Command sliceSolidCmd;

string cmdName = "slice-solid";

static this()
{
    sliceSolidCmd.main = &main_;
    sliceSolidCmd.description = "Slice solid-field snapshots (for structured grids).";
    sliceSolidCmd.shortDescription = sliceSolidCmd.description;
    sliceSolidCmd.helpMsg = format(
`lmr %s [options]

Slice solid-field snapshots across index-directions,
for a selection of field variables.
The really only makes sense for structured-grid blocks.

If no selection of variable names is supplied, then the default action is
to report all field variables (--names=all).

If no options related to snapshot selection are given,
then the default is to process the final snapshot.

options ([+] can be repeated):

 -n, --names
     comma separated list of variable names for reporting
     examples:
       --names=T,e,k
       --names=T
     default:
       --names=all
     The output will always start with pos.x,pos.y and, for 3D, pos.z.

 -l, --slice-list
     slices the solid field in a range of blocks by accepting a string of the form
     "blk-range,i-range,j-range:k-range;"

     examples:
       --slice-list=0:2,:,$,0
       will select blocks 0 and 1, writing out the top strip (j=njc-1) of cells, stepping in i
       Note that you can specify the last value of any index with $.

       --slice-list="0:2,:,0,0;2,$,:,0"
       will select blocks 0 and 1, writing out a strip of cells stepping in i, keeping j=0, k=0, plus
       from block 2, write out a strip of cells stepping in j, keeping i=nic-1 and k=0
       Note that you need the quotes to stop the shell from cutting your command at the semicolon.

     default: none (report nothing)

 -s, --snapshot
     an integer of the snapshot to process for slicing
     examples:
       --snapshot=2 : processes snapshot 2
       -s 4 : processes snapshot 4
     default: none (empty array)


 -f, --final
     process the final snapshot
     If given, this overrides the "-s" option (above).
     default: false

 -o, --output
     write output to a file
     example:
       --output=solid-profile-data.txt
     default: none (just write to STDOUT)

 -v, --verbose [+]
     Increase verbosity.

`, cmdName);

}

int main_(string[] args)
{
    int verbosity = 0;
    int snapshot = -1; // -1 signals no selection set
    bool finalSnapshot = false;
    string namesStr;
    string outFilename;
    string sliceListStr;
    try {
        getopt(args,
               config.bundling,
               "v|verbose+", &verbosity,
               "n|names", &namesStr,
               "s|snapshot", &snapshot,
               "f|final", &finalSnapshot,
               "o|output", &outFilename,
               "l|slice-list", &sliceListStr
               );
    } catch (Exception e) {
        writefln("Eilmer %s program quitting.", cmdName);
        writeln("There is something wrong with the command-line arguments/options.");
        writeln(e.msg);
        return 1;
    }

    if (verbosity > 0) {
        writefln("lmr %s: Begin program.", cmdName);
    }

    if (namesStr.empty) {
        // add default, for when nothing supplied
        namesStr ~= "all";
    }
    // Use stdout if no output filename is supplied,
    // or open a file ready for use if one is.
    File outfile = outFilename.empty() ? stdout : File(outFilename, "w");
    //
    initConfiguration(); // To read in GlobalConfig
    auto availSnapshots = determineAvailableSnapshots();
    // Set to final snapshot if given no other options
    auto snap = to!int(availSnapshots[$-1]);
    // Test if a specific snapshot is given
    if (snapshot >= 0) snap = snapshot;
    //Override other options if final requested
    if (finalSnapshot) snap = to!int(availSnapshots[$-1]);
    if (verbosity > 0) {
        writefln("lmr %s: Slicing solid field for snapshot %s.", cmdName, snap);
    }
    sliceListStr = sliceListStr.strip();
    sliceListStr = sliceListStr.replaceAll(regex("\""), "");
    //
    auto soln = new SolidSolution(to!int(snap), GlobalConfig.nSolidBlocks);
    bool threeD =  canFind(soln.solidBlocks[0].variableNames, "pos.z");
    //
    string[] namesList;
    auto namesVariables = namesStr.split(",");
    foreach (var; namesVariables) {
        var = strip(var);
        if (var.toLower() == "all") {
            foreach (name; soln.solidBlocks[0].variableNames) {
                if (!canFind(["pos.x","pos.y","pos.z"], name)) { namesList ~= name; }
            }
        } else {
            if (canFind(soln.solidBlocks[0].variableNames, var)) {
                namesList ~= var;
            } else {
                writefln("Ignoring requested variable: %s", var);
                writeln("This does not appear in list of solid solution variables.");
            }
        }
    }
    string headerStr = "pos.x pos.y";
    if (threeD) headerStr ~= " pos.z";
    foreach (var; namesList) { headerStr ~= " " ~ var; }
    outfile.writeln(headerStr);
    //
    size_t nFluidBlocks = GlobalConfig.nFluidBlocks;
    size_t nSolidBlocks = GlobalConfig.nSolidBlocks;
    foreach (sliceStr; sliceListStr.split(";")) {
        auto rangeStrings = sliceStr.split(",");
        auto blk_range = decode_range_indices(rangeStrings[0], nFluidBlocks, nFluidBlocks+nSolidBlocks);
        foreach (ib; blk_range[0] .. blk_range[1]) {
            auto sib = ib - nFluidBlocks; // internally in a SolidSolution, solid blocks number from 0
            auto blk = soln.solidBlocks[sib];
            // We need to do the decode in the context of each block because
            // the upper limits to the indices are specific to the block.
            auto i_range = decode_range_indices(rangeStrings[1], 0, blk.nic);
            auto j_range = decode_range_indices(rangeStrings[2], 0, blk.njc);
            auto k_range = decode_range_indices(rangeStrings[3], 0, blk.nkc);
            foreach (k; k_range[0] .. k_range[1]) {
                foreach (j; j_range[0] .. j_range[1]) {
                    foreach (i; i_range[0] .. i_range[1]) {
                        outfile.write(format(" %g %g", soln.solidBlocks[sib]["pos.x", i, j, k],
                                             soln.solidBlocks[sib]["pos.y", i, j, k]));
                        if (threeD) {
                            outfile.write(format(" %g", soln.solidBlocks[sib]["pos.z", i, j, k]));
                        }
                        foreach (var; namesList) {
                            outfile.write(format(" %g", soln.solidBlocks[sib][var, i, j, k]));
                        }
                        outfile.write("\n");
                    }
                }
            }
        } // end foreach ib
    } // end foreach sliceStr

    if (verbosity > 0) {
        writefln("lmr %s: Done.", cmdName);
    }
    outfile.close();
    return 0;
}
