/**
 * Module for converting Eilmer-native flow fields to VTK format.
 *
 * Authors: RJG, PJ, KAD, NNG
 * Date: 2022-08-08
 */

module snapshot2vtk;

import std.stdio;
import std.file;
import std.format;
import std.getopt;
import std.conv : to;
import std.range;
import std.string;
import std.regex;

import globalconfig;
import fileutil;
import flowsolution;
import solidsolution : SolidSolution;
import lmrconfig : lmrCfg;
import init : initConfiguration;
import vtk_writer;
import cmdhelper;

import command;

Command snapshot2vtkCmd;
string cmdName = "snapshot2vtk";

static this()
{
    snapshot2vtkCmd.main = &main_;
    snapshot2vtkCmd.description = "Convert snapshots of flow fields to VTK format.";
    snapshot2vtkCmd.shortDescription = snapshot2vtkCmd.description;
    snapshot2vtkCmd.helpMsg =
`lmr snapshot2vtk [options]

Convert a flow field using one or more snapshots to VTK format.

If no options related to snapshot selection are given,
then the default is to process the final snapshot.

options ([+] can be repeated):

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

 -b, --binary-format
     selects binary format for output
     default: false

 --add-vars
     comma separated array of auxiliary variables to add in VTK
     eg. --add-vars="mach,pitot"
     Other variables include:
         total-h, total-p, total-T,
         enthalpy, entropy, molef, conc,
         Tvib (for some gas models)
         nrf (non-rotating-frame velocities)
         cyl (cylindrical coordinates: r, theta)

 -r, --subtract-ref-solution
     name of file containing a Lua-format reference solution
     matching field variables in solution files and reference solution
     will be altered as: numerical value - reference value
     examples:
     --subtract-ref-solution=ref-soln.lua
     default: none

 --subtract-solid-ref-solution
     name of file containing a Lua-format reference solution for the solid domains
     matching field variables in solution files and reference solution
     will be altered as: numerical value - reference value
     examples:
     --subtract-solid-ref-solution=solid-ref-soln.lua
     default: none

 -v, --verbose [+]
     Increase verbosity during preparation and writing of VTK files.

`;

}

int main_(string[] args)
{
    int verbosity = 0;
    int[] snapshots;
    bool finalSnapshot = false;
    bool allSnapshots = false;
    bool binaryFormat = false;
    string luaRefSoln;
    string luaSolidRefSoln;
    string addVarsStr;
    try {
        getopt(args,
               config.bundling,
               "v|verbose+", &verbosity,
               "s|snapshots|snapshot", &snapshots,
               "f|final", &finalSnapshot,
               "a|all", &allSnapshots,
               "b|binary-format", &binaryFormat,
               "r|subtract-ref-solution", &luaRefSoln,
               "subtract-solid-ref-solution", &luaSolidRefSoln,
               "add-vars", &addVarsStr);
    } catch (Exception e) {
        writefln("Eilmer %s program quitting.", cmdName);
        writeln("There is something wrong with the command-line arguments/options.");
        writeln(e.msg);
        return 1;
    }

    if (verbosity > 0) {
        writeln("lmr snapshot2vtk: Begin program.");
    }

    initConfiguration(); // To read in GlobalConfig
    size_t nFluidBlocks = GlobalConfig.nFluidBlocks;
    auto availSnapshots = determineAvailableSnapshots();
    auto snaps2process = determineSnapshotsToProcess(availSnapshots, snapshots, allSnapshots, finalSnapshot);

    string[] addVarsList;
    addVarsStr = addVarsStr.strip();
    addVarsStr = addVarsStr.replaceAll(regex("\""), "");
    if (addVarsStr.length > 0) {
        addVarsList = addVarsStr.split(",");
    }

    /*
     * When doing time-marching, also add timestep information.
     */
    double[] times;
    if (GlobalConfig.solverMode == SolverMode.transient) {
        times = mapTimesToSnapshots(snaps2process);
    }

    /*
     * We should clean out if the vtk/ area already has stuff in it.
     * It gets confusing when the meta files don't align with the pieces
     * sitting around in the directory.
     */
    if (lmrCfg.vtkDir.exists) {
        if (verbosity > 1) { writeln("lmr snapshot2vtk: Removing old VTK files."); }
        lmrCfg.vtkDir.rmdirRecurse;
    }

    /*
     * Now write vtk files for each snapshot
     */
    if (verbosity > 0) {
        writeln("lmr snapshot2vtk: Writing VTK files to disk.");
    }

    ensure_directory_is_present(lmrCfg.vtkDir);
    File pvdFile = begin_PVD_file(lmrCfg.vtkDir~"/"~lmrCfg.fluidPrefix~".pvd");
    foreach (isnap, snap; snaps2process) {
        if (verbosity > 1) {
            writefln("lmr snapshot2vtk: Writing fluid snapshot %s to disk.", snap);
        }
        double simTime = (GlobalConfig.solverMode == SolverMode.transient) ? times[isnap] : -1.0;
        auto soln = new FlowSolution(to!int(snap), GlobalConfig.nFluidBlocks, simTime);
        soln.add_aux_variables(addVarsList);
        if (!luaRefSoln.empty) soln.subtract_ref_soln(luaRefSoln);
        // [TODO] add aux variables.
        string pvtuFileName = lmrCfg.fluidPrefix ~ "-s" ~ snap ~ ".pvtu";
        File pvtuFile = begin_PVTU_file(lmrCfg.vtkDir ~ "/" ~ pvtuFileName, soln.flowBlocks[0].variableNames);
        foreach (jb; 0 .. soln.nBlocks) {
            if (verbosity > 2) {
                writefln("lmr snapshot2vtk: Writing fluid block %d for snapshot %s to disk.", jb, snap);
            }
            string vtuFileName = lmrCfg.fluidPrefix ~ "-s" ~ snap ~ "-b" ~ format(lmrCfg.blkIdxFmt, jb) ~ ".vtu";
            add_piece_to_PVTU_file(pvtuFile, vtuFileName);
            write_VTU_file(soln.flowBlocks[jb], soln.gridBlocks[jb], lmrCfg.vtkDir~"/"~vtuFileName, binaryFormat);
            if (GlobalConfig.solverMode == SolverMode.transient) {
                add_dataset_to_PVD_file(pvdFile, times[isnap], vtuFileName);
            }
            else {
                add_dataset_to_PVD_file(pvdFile, to!double(snap), vtuFileName);
            }
        }
        finish_PVTU_file(pvtuFile);
    }
    finish_PVD_file(pvdFile);

    if (GlobalConfig.nSolidBlocks > 0) {
        pvdFile = begin_PVD_file(lmrCfg.vtkDir~"/"~lmrCfg.solidPrefix~".pvd");
        foreach (isnap, snap; snaps2process) {
            if (verbosity > 1) {
                writefln("lmr snapshot2vtk: Writing solid snapshot %s to disk.", snap);
            }
            double simTime = (GlobalConfig.solverMode == SolverMode.transient) ? times[isnap] : -1.0;
            auto soln = new SolidSolution(to!int(snap), GlobalConfig.nSolidBlocks, simTime);
            if (!luaSolidRefSoln.empty) soln.subtract_ref_soln(luaSolidRefSoln);
            string pvtuFileName = lmrCfg.solidPrefix ~ "-s" ~ snap ~ ".pvtu";
            File pvtuFile = begin_PVTU_file(lmrCfg.vtkDir ~ "/" ~ pvtuFileName, soln.solidBlocks[0].variableNames);
            foreach (jb; 0 .. soln.nBlocks) {
                auto id = jb + nFluidBlocks;
                if (verbosity > 2) {
                    writefln("lmr snapshot2vtk: Writing solid block %d for snapshot %s to disk.", id, snap);
                }
                string vtuFileName = lmrCfg.solidPrefix ~ "-s" ~ snap ~ "-b" ~ format(lmrCfg.blkIdxFmt, id) ~ ".vtu";
                add_piece_to_PVTU_file(pvtuFile, vtuFileName);
                write_VTU_file(soln.solidBlocks[jb], soln.gridBlocks[jb], lmrCfg.vtkDir~"/"~vtuFileName, binaryFormat);
                if (GlobalConfig.solverMode == SolverMode.transient) {
                    add_dataset_to_PVD_file(pvdFile, times[isnap], vtuFileName);
                }
                else {
                    add_dataset_to_PVD_file(pvdFile, to!double(snap), vtuFileName);
                }
            }
            finish_PVTU_file(pvtuFile);
        }
        finish_PVD_file(pvdFile);
    }

    if (verbosity > 0) {
        writeln("lmr snapshot2vtk: Done.");
    }
    return 0;
}


