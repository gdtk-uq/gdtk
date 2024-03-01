/**
 * Module for slicing flow-field snapshots, picking out specific data.
 *
 * Authors: RJG, PJ, KAD, NNG
 * Date: 2024-03-01, adapted from probeflow.d and src/eilmer/postprocess.d
 */

module sliceflow;

import std.stdio;
import std.file;
import std.string;
import std.regex : regex, replaceAll;
import std.array;
import std.format : format;
import std.getopt;
import std.conv : to;
import std.range;
import std.algorithm;

import globalconfig;
import fileutil;
import flowsolution;
import lmrconfig : lmrCfg;
import init : initConfiguration;
import cmdhelper;

import command;

Command sliceFlowCmd;

string cmdName = "slice-flow";

static this()
{
    sliceFlowCmd.main = &main_;
    sliceFlowCmd.description = "Slice flow-field snapshots (for structured grids).";
    sliceFlowCmd.shortDescription = sliceFlowCmd.description;
    sliceFlowCmd.helpMsg = format(
`lmr %s [options]

Slice flow-field snapshots across index-directions,
for a selection of flow field variables.
The really only makes sense for structured-grid blocks.

If no selection of variable names is supplied, then the default
is to use --names="rho,p,T,vel.x,vel.y"

If no options related to snapshot selection are given,
then the default is to process the final snapshot.

options ([+] can be repeated):

 -n, --names
     comma separated list of variable names for reporting
     examples:
       --names="rho"
     default:
       --names="rho,p,T,vel.x,vel.y"
     The output will always start with pos.x,pos.y,pos.z

 -l, --slice-list
     slices the flow field in a range of blocks by accepting a string of the form
     "blk-range,i-range,j-range:k-range;"

     example:
       --slice-list="0:2,:,0,0;2,$,:,0"
       will select blocks 0 and 1, writing out a strip of cells stepping in i, keeping j=0, k=0, plus
       from block 2, write out a strip of cells stepping in j, keeping i=nic-1 and k=0
       Note that you can specify the last value of an index with $.

     default: none (report nothing)

 -s, --snapshot[s]+
     comma separated array of snapshots to convert
     examples:
       --snapshots=0,1,5 : processes snapshots 0, 1 and 5
       --snapshot=2 : processes snapshot 2 only
       --snapshot 1  --snapshot 4 : processes snapshots 1 and 4
     default: none (empty array)


 -f, --final
     process the final snapshot
     default: false

 -a, --all
     process all snapshots
     default: false

 -o, --output
     write output to a file
     example:
       --output=norms-data.txt
     default: none (just write to STDOUT)

 -v, --verbose [+]
     Increase verbosity.

`, cmdName);

}

int main_(string[] args)
{
    int verbosity = 0;
    int[] snapshots;
    bool finalSnapshot = false;
    bool allSnapshots = false;
    bool binaryFormat = false;
    string namesStr;
    string outFilename;
    string sliceListStr;
    string luaRefSoln;
    try {
        getopt(args,
               config.bundling,
               "v|verbose+", &verbosity,
               "n|names", &namesStr,
               "s|snapshots|snapshot", &snapshots,
               "f|final", &finalSnapshot,
               "a|all", &allSnapshots,
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

    if (namesStr.empty) { // add default, for when nothing supplied
        namesStr ~= "rho,p,T,vel.x,vel.y";
    }

    auto namesVariables = namesStr.split(",");
    foreach (var; namesVariables) var = strip(var);

    // Use stdout if no output filename is supplied,
    // or open a file ready for use if one is.
    File outfile = outFilename.empty() ? stdout : File(outFilename, "w");

    initConfiguration(); // To read in GlobalConfig
    auto availSnapshots = determineAvailableSnapshots();
    auto snaps2process = determineSnapshotsToProcess(availSnapshots, snapshots, allSnapshots, finalSnapshot);

    if (verbosity > 0) {
        writefln("lmr %s: Slicing flow field.", cmdName);
    }

    sliceListStr = sliceListStr.strip();
    sliceListStr = sliceListStr.replaceAll(regex("\""), "");

    foreach (snap; snaps2process) {
        if (verbosity > 1) {
            writefln("lmr %s: Slicing flow field for snapshot %s.", cmdName, snap);
        }
        outfile.writefln("# snapshot: %s", snap);
        auto soln = new FlowSolution(to!int(snap), GlobalConfig.nFluidBlocks);
        bool is3dimensional =  canFind(soln.flowBlocks[0].variableNames, "pos.z");
        outfile.write("# pos.x pos.y");
        if (is3dimensional) outfile.write(" pos.z");
        foreach (var; namesVariables) {
            if (!canFind(soln.flowBlocks[0].variableNames, var)) {
                writefln("Ignoring requested variable %q.", var);
                writeln("This does not appear in list of flow solution variables.");
                continue;
            }
            outfile.write(" ", var);
        }
        outfile.write("\n");
        foreach (sliceStr; sliceListStr.split(";")) {
            auto rangeStrings = sliceStr.split(",");
            auto blk_range = decode_range_indices(rangeStrings[0], 0, soln.nBlocks);
            foreach (ib; blk_range[0] .. blk_range[1]) {
                auto blk = soln.flowBlocks[ib];
                // We need to do the decode in the context of each block because
                // the upper limits to the indices are specific to the block.
                auto i_range = decode_range_indices(rangeStrings[1], 0, blk.nic);
                auto j_range = decode_range_indices(rangeStrings[2], 0, blk.njc);
                auto k_range = decode_range_indices(rangeStrings[3], 0, blk.nkc);
                foreach (k; k_range[0] .. k_range[1]) {
                    foreach (j; j_range[0] .. j_range[1]) {
                        foreach (i; i_range[0] .. i_range[1]) {
                            outfile.write(format(" %g %g", soln.flowBlocks[ib]["pos.x", i, j, k],
                                                 soln.flowBlocks[ib]["pos.y", i, j, k]));
                            if (is3dimensional) {
                                outfile.write(format(" %g", soln.flowBlocks[ib]["pos.z", i, j, k]));
                            }
                            foreach (var; namesVariables) {
                                if (canFind(soln.flowBlocks[0].variableNames, var)) {
                                    outfile.write(format(" %g", soln.flowBlocks[ib][var, i, j, k]));
                                }
                            }
                            outfile.write("\n");
                        }
                    }
                }
            } // end foreach ib
        } // end foreach sliceStr
        outfile.writeln(""); // A blank line between snapshots
    }

    if (verbosity > 0) {
        writefln("lmr %s: Done.", cmdName);
    }
    outfile.close();
    return 0;
}


size_t[] decode_range_indices(string rangeStr, size_t first, size_t endplus1)
// Decode strings such as "0:$", ":", "0:3", "$"
// On input, first and endplus1 represent the largest, available range.
// Return the pair of numbers that can be used in a foreach loop range.
{
    if (rangeStr == ":") {
        return [first, endplus1];
    }
    if (canFind(rangeStr, ":")) {
        // We have a range specification to pull apart.
        auto items = rangeStr.split(":");
        first = to!size_t(items[0]);
        if (items.length > 1 && items[1] != "$") {
            // Presume that we have a second integer.
            size_t new_endplus1 = to!size_t(items[1]);
            if (new_endplus1 < endplus1) endplus1 = new_endplus1;
        }
    } else if (rangeStr == "$") {
        // With just a single "$" specified, we want only the last index.
        first = endplus1 - 1;
    }else {
        // Presume that we have a single integer.
        first = to!size_t(rangeStr);
        if (first < endplus1) endplus1 = first+1;
    }
    return [first, endplus1];
} // end decode_range_indices()

