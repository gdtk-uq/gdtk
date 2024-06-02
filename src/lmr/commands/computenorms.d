/**
 * Module for computing error norms based on a reference solution.
 *
 * Authors: RJG, PJ, KAD, NNG
 * Date: 2023-07-14
 */

module computenorms;

import std.stdio;
import std.file;
import std.string;
import std.array;
import std.format : format;
import std.getopt;
import std.conv : to;
import std.range;
import std.algorithm;

import globalconfig;
import fileutil;
import flowsolution;
import solidsolution;
import lmrconfig : lmrCfg;
import init : initConfiguration;
import cmdhelper;

import command;

Command compNormsCmd;

string cmdName = "compute-norms";

static this()
{
    compNormsCmd.main = &main_;
    compNormsCmd.description = "Compute field norms (possibly with respect to a reference solution).";
    compNormsCmd.shortDescription = compNormsCmd.description;
    compNormsCmd.helpMsg = format(
`lmr %s [options]

Compute norms for a snapshot based on selection of field variables.

If a reference solution provided (as Lua file), then use that to compute error norms for
selected field variables.

If no selection of norm variables is supplied, then the default action is
to compute a density norm (--norms="rho").

If no options related to snapshot selection are given,
then the default is to process the final snapshot.

options ([+] can be repeated):

 -n, --norms
     comma separated list of variables for norms calculation
     examples:
       --norms="rho,velx,vely"
       --norms="rho"
     default:
       --norms="rho"
       If no norms-list supplied, then just process density.

 --solid-norms
     comma separated list of variables for norms calculation in solid domains
     examples:
       --solid-norms="T,e"
     default: none (don't do anything for solid domain)

 -r, --reference-solution
     a reference solution provided in a Lua file

     the reference solution will be subtracted from the field variables
     depending on which variables are provided in the reference solution

     examples:
       --reference-solution=my_ref_soln.lua
       -r ref-solution.lua
     default: none (just compute field norms)

 --limit-to-region
     limits the norm computation to a box with lower corner (x0,y0,z0)
     and upper corner (x1,y1,z1) by accepting a string of the form
     "x0,y0,z0,x1,y1,z1"

     this is sometimes useful to focus the error estimate on the interior
     away from boundaries

     for 2D, simply supply z = 0.0 for z0 and z1
     examples:
       --limit-to-region="0.25,0.25,0.0,0.75,0.75,0.0"

     default: none (just use entire domain)

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
    string normsStr;
    string solidNormsStr;
    string outFilename;
    string region;
    string luaRefSoln;
    try {
        getopt(args,
               config.bundling,
               "v|verbose+", &verbosity,
               "n|norms", &normsStr,
               "solid-norms", &solidNormsStr,
               "r|reference-solution", &luaRefSoln,
               "s|snapshots|snapshot", &snapshots,
               "f|final", &finalSnapshot,
               "a|all", &allSnapshots,
               "o|output", &outFilename,
               "limit-to-region", &region);
    } catch (Exception e) {
        writefln("Eilmer %s program quitting.", cmdName);
        writeln("There is something wrong with the command-line arguments/options.");
        writeln(e.msg);
        return 1;
    }

    if (verbosity > 0) {
        writefln("lmr %s: Begin program.", cmdName);
    }

    if (normsStr.empty) { // add default of "rho" when nothing supplied
        normsStr ~= "rho";
    }

    auto normsVariables = normsStr.split(",");
    foreach (var; normsVariables) var = strip(var);
    auto solidNormsVariables = solidNormsStr.split(",");
    foreach (var; solidNormsVariables) var = strip(var);

    // Use stdout if no output filename is supplied,
    // or open a file ready for use if one is.
    File outfile = outFilename.empty() ? stdout : File(outFilename, "w");

    initConfiguration(); // To read in GlobalConfig
    auto availSnapshots = determineAvailableSnapshots();
    auto snaps2process = determineSnapshotsToProcess(availSnapshots, snapshots, allSnapshots, finalSnapshot);

    /*
     * When doing time-marching, also add timestep information.
     */
    double[] times;
    if (GlobalConfig.solverMode == SolverMode.transient) {
        times = mapTimesToSnapshots(snaps2process);
    }

    if (verbosity > 0) {
        writefln("lmr %s: Computing norms.", cmdName);
    }

    foreach (isnap, snap; snaps2process) {
        if (verbosity > 1) {
            writefln("lmr %s: Computing norms for fluid snapshot %s.", cmdName, snap);
        }
        double simTime = (GlobalConfig.solverMode == SolverMode.transient) ? times[isnap] : -1.0;
        auto soln = new FlowSolution(to!int(snap), GlobalConfig.nFluidBlocks, simTime);
        outfile.writeln("---"); // YAML document opener
        outfile.writefln("snapshot: %s", snap);
        outfile.writeln("field_type: fluid");
        if (!luaRefSoln.empty) {
            soln.subtract_ref_soln(luaRefSoln);
            outfile.writefln("error-norms-computed: yes");
            outfile.writefln("reference-solution-file: %s", luaRefSoln);
        }
        foreach (var; normsVariables) {
            if (verbosity > 2) {
                writefln("lmr %s: Computing norm for variable= %s", cmdName, var);
            }
            if (!canFind(soln.flowBlocks[0].variableNames, var)) {
                writeln(format("Ignoring requested variable \"%s\". This does not appear in list of flow solution variables.", var));
                continue;
            }
            auto norms = soln.compute_volume_weighted_norms(var, region);
            outfile.writefln("%s:", var);
            outfile.writefln("   L1: " ~ lmrCfg.dblVarFmt, norms[0]);
            outfile.writefln("   L2: " ~ lmrCfg.dblVarFmt, norms[1]);
            outfile.writefln("   Linf: " ~ lmrCfg.dblVarFmt, norms[2]);
            outfile.writefln("   Linf-pos: {x: %.6e, y: %.6e, z: %.6e }", norms[3], norms[4], norms[5]);
        }
        outfile.writeln("..."); // YAML document closer
    }

    if (GlobalConfig.nSolidBlocks > 0 && solidNormsVariables.length > 0) {
        foreach (isnap, snap; snaps2process) {
            if (verbosity > 1) {
                writefln("lmr %s: Computing norms for solid snapshot %s.", cmdName, snap);
            }
            double simTime = (GlobalConfig.solverMode == SolverMode.transient) ? times[isnap] : -1.0;
            auto soln = new SolidSolution(to!int(snap), GlobalConfig.nSolidBlocks, simTime);
            outfile.writeln("---"); // YAML document opener
            outfile.writefln("snapshot: %s", snap);
            outfile.writeln("field_type: solid");
            if (!luaRefSoln.empty) {
                soln.subtract_ref_soln(luaRefSoln);
                outfile.writefln("error-norms-computed: yes");
                outfile.writefln("reference-solution-file: %s", luaRefSoln);
            }
            foreach (var; solidNormsVariables) {
                if (verbosity > 2) {
                    writefln("lmr %s: Computing norm for variable= %s", cmdName, var);
                }
                if (!canFind(soln.solidBlocks[0].variableNames, var)) {
                    writeln(format("Ignoring requested variable \"%s\". This does not appear in list of solid solution variables.", var));
                    continue;
                }
                auto norms = soln.compute_volume_weighted_norms(var, region);
                outfile.writefln("%s:", var);
                outfile.writefln("   L1: " ~ lmrCfg.dblVarFmt, norms[0]);
                outfile.writefln("   L2: " ~ lmrCfg.dblVarFmt, norms[1]);
                outfile.writefln("   Linf: " ~ lmrCfg.dblVarFmt, norms[2]);
                outfile.writefln("   Linf-pos: {x: %.6e, y: %.6e, z: %.6e }", norms[3], norms[4], norms[5]);
            }
            outfile.writeln("..."); // YAML document closer
        }
    }

    if (verbosity > 0) {
        writefln("lmr %s: Done.", cmdName);
    }
    outfile.close();
    return 0;
}


