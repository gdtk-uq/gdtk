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

import globalconfig;
import fileutil;
import flowsolution;
import lmrconfig : lmrCfg;
import init : initConfiguration;
import vtk_writer;
import cmdhelper;

import command;

Command snapshot2vtkCmd;

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
     [TODO] not implemented presently

 -r, --subtract-ref-solution
     name of file containing a Lua-format reference solution
     matching field variables in solution files and reference solution
     will be altered as: numerical value - reference value
     examples:
     --subtract-ref-solution=ref-soln.lua
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
    // [TODO] implement --add-vars
    getopt(args,
           config.bundling,
           "v|verbose+", &verbosity,
           "s|snapshots|snapshot", &snapshots,
           "f|final", &finalSnapshot,
           "a|all", &allSnapshots,
           "b|binary-format", &binaryFormat,
           "r|subtract-ref-solution", &luaRefSoln);

    if (verbosity > 0) {
        writeln("lmr snapshot2vtk: Begin program.");
    }

    initConfiguration(); // To read in GlobalConfig
    auto availSnapshots = determineAvailableSnapshots();
    auto snaps2process = determineSnapshotsToProcess(availSnapshots, snapshots, allSnapshots, finalSnapshot);

    /*
     * Now write vtk files for each snapshot
     */
    if (verbosity > 0) {
        writeln("lmr snapshot2vtk: Writing VTK files to disk.");
    }

    ensure_directory_is_present(lmrCfg.vtkDir);
    File pvdFile = begin_PVD_file(lmrCfg.vtkDir~"/"~lmrCfg.flowPrefix~".pvd");
    foreach (snap; snaps2process) {
        if (verbosity > 1) {
            writefln("lmr snapshot2vtk: Writing snapshot %s to disk.", snap);
        }
        auto soln = new FlowSolution(to!int(snap), GlobalConfig.nFluidBlocks);
        if (!luaRefSoln.empty) soln.subtract_ref_soln(luaRefSoln);
        // [TODO] add aux variables.
        string pvtuFileName = lmrCfg.flowPrefix ~ snap ~ ".pvtu";
        File pvtuFile = begin_PVTU_file(lmrCfg.vtkDir ~ "/" ~ pvtuFileName, soln.flowBlocks[0].variableNames);
        foreach (jb; 0 .. soln.nBlocks) {
            if (verbosity > 2) {
                writefln("lmr snapshot2vtk: Writing block %d for snapshot %s to disk.", jb, snap);
            }
            string vtuFileName = lmrCfg.flowPrefix ~ format(lmrCfg.blkIdxFmt, jb) ~ "-" ~ snap ~ ".vtu";
            add_dataset_to_PVD_file(pvdFile, to!double(snap), vtuFileName);
            add_piece_to_PVTU_file(pvtuFile, vtuFileName);
            write_VTU_file(soln.flowBlocks[jb], soln.gridBlocks[jb], lmrCfg.vtkDir~"/"~vtuFileName, binaryFormat);
        }
        finish_PVTU_file(pvtuFile);
    }
    finish_PVD_file(pvdFile);

    if (verbosity > 0) {
        writeln("lmr snapshot2vtk: Done.");
    }
    return 0;
}


