/**
 * Module for probing flow-field snapshots, picking out specific data.
 *
 * Authors: RJG, PJ, KAD, NNG
 * Date: 2024-02-18, adapted from computenorms.d and src/eilmer/postprocess.d
 */

module probeflow;

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
import std.string;
import std.regex;
import std.uni;

import globalconfig;
import fileutil;
import flowsolution;
import lmrconfig : lmrCfg;
import init : initConfiguration;
import cmdhelper;

import command;

Command probeFlowCmd;

string cmdName = "probe-flow";

static this()
{
    probeFlowCmd.main = &main_;
    probeFlowCmd.description = "Probe flow-field snapshots.";
    probeFlowCmd.shortDescription = probeFlowCmd.description;
    probeFlowCmd.helpMsg = format(
`lmr %s [options]

Probe flow-field snapshots based on a selection of field variables.

If no selection of variable names is supplied, then the default action is
to probe all flow field variables (--names=all).

If no options related to snapshot selection are given,
then the default is to process the final snapshot.

options ([+] can be repeated):

 -n, --names
     comma separated list of variable names for reporting
     examples:
       --names=rho,vel.x,vel.y
       --names=rho
     default:
       --names=all

 --add-vars
     comma separated array of auxiliary variables to add to the flow solution
     eg. --add-vars=mach,pitot
     Other variables include:
         total-h, total-p, total-T,
         enthalpy, entropy, molef, conc,
         Tvib (for some gas models)
         nrf (non-rotating-frame velocities)
         cyl (cylindrical coordinates: r, theta)

 -l, --location
     probes the flow field at locations 0, 1 and 2 by accepting a string of the form
     "x0,y0,z0;x1,y1,z1;x2,y2,z2"

     For 2D, simply supply 0 for z, or omit.
     To avoid the shell cutting your command at the semicolon, use quotes.
     example:
       --location="0.25,0.25,0.0;0.75,0.75"

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
    string locationStr;
    string luaRefSoln;
    string addVarsStr;
    try {
        getopt(args,
               config.bundling,
               "v|verbose+", &verbosity,
               "n|names", &namesStr,
               "s|snapshots|snapshot", &snapshots,
               "f|final", &finalSnapshot,
               "a|all", &allSnapshots,
               "o|output", &outFilename,
               "l|location", &locationStr,
               "add-vars", &addVarsStr
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

    if (namesStr.empty) { // add default of "all" when nothing supplied
        namesStr ~= "all";
    }

    // Use stdout if no output filename is supplied,
    // or open a file ready for use if one is.
    File outfile = outFilename.empty() ? stdout : File(outFilename, "w");

    initConfiguration(); // To read in GlobalConfig
    auto availSnapshots = determineAvailableSnapshots();
    auto snaps2process = determineSnapshotsToProcess(availSnapshots, snapshots, allSnapshots, finalSnapshot);

    string[] addVarsList;
    addVarsStr = addVarsStr.strip();
    addVarsStr = addVarsStr.replaceAll(regex("\""), "");
    if (addVarsStr.length > 0) {
        addVarsList = addVarsStr.split(",");
    }

    if (verbosity > 0) {
        writefln("lmr %s: Probing flow field.", cmdName);
    }

    locationStr = locationStr.strip();
    locationStr = locationStr.replaceAll(regex("\""), "");
    double[] xp, yp, zp;
    foreach(triple; locationStr.split(";")) {
        auto items = triple.split(",");
        xp ~= to!double(items[0]);
        yp ~= to!double(items[1]);
        zp ~= (items.length > 2) ? to!double(items[2]) : 0.0;
    }
    if (xp.length == 0) {
        writeln("No probe locations given.");
    }

    foreach (snap; snaps2process) {
        if (verbosity > 1) {
            writefln("lmr %s: Probing flow field for snapshot %s.", cmdName, snap);
        }
        auto soln = new FlowSolution(to!int(snap), GlobalConfig.nFluidBlocks);
        soln.add_aux_variables(addVarsList);
        //
        string[] namesList;
        auto namesVariables = namesStr.split(",");
        foreach (var; namesVariables) {
            var = strip(var);
            if (var.toLower() == "all") {
                foreach (name; soln.flowBlocks[0].variableNames) namesList ~= name;
            } else {
                namesList ~= var;
            }
        }
        outfile.writeln("---"); // YAML document opener
        outfile.writefln("snapshot: %s", snap);
        foreach (ip; 0 .. xp.length) {
            if (verbosity > 2) {
                writefln("lmr %s: Probing flow field at location [%g,%g,%g]",
                         cmdName, xp[ip], yp[ip], zp[ip]);
            }
            auto nearest = soln.find_nearest_cell_centre(xp[ip], yp[ip], zp[ip]);
            size_t iblk = nearest[0]; size_t icell = nearest[1];
            // Each probe point becomes a YAML list item and
            // all of the data for this probe point is in a YAML map.
            outfile.writeln("pointdata:");
            outfile.writefln("  - x: %s\n    y: %s\n    z: %s",
                             soln.get_value_str(iblk, icell, "pos.x"),
                             soln.get_value_str(iblk, icell, "pos.y"),
                             soln.get_value_str(iblk, icell, "pos.z"));
            foreach (var; namesList) {
                if (verbosity > 2) {
                    writefln("lmr %s: Probing flow field for variable= %s", cmdName, var);
                }
                if (!canFind(soln.flowBlocks[0].variableNames, var)) {
                    writefln("Ignoring requested variable: %s", var);
                    writeln("This does not appear in list of flow solution variables:");
                    writeln("  ", soln.flowBlocks[0].variableNames);
                    continue;
                }
                outfile.writefln("    %s: %s", var, soln.get_value_str(iblk, icell, var));
            }
        }
        outfile.writeln("..."); // YAML document closer
    }

    if (verbosity > 0) {
        writefln("lmr %s: Done.", cmdName);
    }
    outfile.close();
    return 0;
}


