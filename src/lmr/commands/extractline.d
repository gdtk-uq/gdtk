/**
 * Module for postprocessing flow-field snapshots, picking out lines of data.
 *
 * Authors: RJG, PJ, KAD, NNG
 * Date: 2024-03-04, adapted from sliceflow.d and src/eilmer/postprocess.d
 */

module extractline;

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

import geom : Vector3;
import globalconfig;
import fileutil;
import flowsolution;
import lmrconfig : lmrCfg;
import init : initConfiguration;
import cmdhelper;

import command;

Command extractLineCmd;

string cmdName = "extract-line";

static this()
{
    extractLineCmd.main = &main_;
    extractLineCmd.description = "Extract lines of data from flow-field snapshots.";
    extractLineCmd.shortDescription = extractLineCmd.description;
    extractLineCmd.helpMsg = format(
`lmr %s [options]

Extract lines of data from flow-field snapshots,
for a selection of flow field variables.
Works on structured-grid and unstructured-grid blocks.

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

 -l, --line-list
     lines from point 0 to point 1 and sample numbers are specified as a string of the form
     "x0,y0,z0,x1,y1,z1,n;"
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
    string lineListStr;
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
               "l|line-list", &lineListStr
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

    lineListStr = lineListStr.strip();
    lineListStr = lineListStr.replaceAll(regex("\""), "");

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
        size_t[2][] cells_found; // accumulate the identies of the cells found here
        foreach (lineStr; lineListStr.split(";")) {
            auto items = lineStr.split(",");
            if (items.length != 7) {
                string errMsg = "The 'extract-line' string requires exactly 7 values.\n";
                errMsg ~= format("You have provided %d items.\n", items.length);
                errMsg ~= format("The problematic string is: %s\n", lineStr);
                throw new Error(errMsg);
            }
            Vector3 p0 = Vector3(to!double(items[0]), to!double(items[1]), to!double(items[2]));
            Vector3 p1 = Vector3(to!double(items[3]), to!double(items[4]), to!double(items[5]));
            size_t n = to!size_t(items[6]);
            auto count = soln.find_enclosing_cells_along_line(p0, p1, n, cells_found);
            writeln("# Info: Found ", count, " cells from point ", p0, " to point ", p1);
        } // end foreach lineStr
        foreach(i; 0 .. cells_found.length) {
            size_t ib = cells_found[i][0]; size_t idx = cells_found[i][1];
            outfile.write(format(" %g %g", soln.flowBlocks[ib]["pos.x", idx],
                                 soln.flowBlocks[ib]["pos.y", idx]));
            if (is3dimensional) {
                outfile.write(format(" %g", soln.flowBlocks[ib]["pos.z", idx]));
            }
            foreach (var; namesVariables) {
                if (canFind(soln.flowBlocks[0].variableNames, var)) {
                    outfile.write(format(" %g", soln.flowBlocks[ib][var, idx]));
                }
            }
            outfile.write("\n");
        }
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

