/**
 * Module for converting Eilmer-native residual fields to VTK format.
 *
 * Authors: RJG, PJ, KAD, NNG
 * Date: 2024-03-07
 */

module residual2vtk;

import std.stdio;
import std.file;
import std.format;
import std.getopt;
import std.conv;
import std.range;
import std.bitmanip;
import std.stdint;
import std.algorithm.mutation : remove;

import geom;
import globalconfig;
import fileutil;
import flowsolution;
import lmrconfig;
import init : initConfiguration;
import vtk_writer;
import cmdhelper;
import blockio;

import command;

Command residual2vtkCmd;
string cmdName = "residual2vtk";

static this()
{
    residual2vtkCmd.main = &main_;
    residual2vtkCmd.description = "Convert fields of residual values to VTK format.";
    residual2vtkCmd.shortDescription = residual2vtkCmd.description;
    residual2vtkCmd.helpMsg = format(
`lmr %s [options]

Convert the residual values for flow field using one or more snapshots to VTK format.

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

 -v, --verbose [+]
     Increase verbosity during preparation and writing of VTK files.

`, cmdName);

}

int main_(string[] args)
{
    double[][] data;
    string[] variables;
    string fileFmt;
    int nBlocks;

    int verbosity = 0;
    int[] snapshots;
    bool finalSnapshot = false;
    bool allSnapshots = false;
    bool binaryFormat = false;
    try {
        getopt(args,
               config.bundling,
               "v|verbose+", &verbosity,
               "s|snapshots|snapshot", &snapshots,
               "f|final", &finalSnapshot,
               "a|all", &allSnapshots,
               "b|binary-format", &binaryFormat);
    } catch (Exception e) {
        writefln("Eilmer %s program quitting.", cmdName);
        writeln("There is something wrong with the command-line arguments/options.");
        writeln(e.msg);
        return 1;
    }

    if (verbosity > 0) writefln("lmr %s: Begin program.", cmdName);

    initConfiguration(); // To read in GlobalConfig
    nBlocks = GlobalConfig.nFluidBlocks;
    fileFmt = GlobalConfig.field_format;
    variables = readVariablesFromMetadata(lmrCfg.residualMetadataFile); 
    auto availSnapshots = determineAvailableSnapshots();
    auto snaps2process = determineSnapshotsToProcess(availSnapshots, snapshots, allSnapshots, finalSnapshot);

    /*
     * Now write vtk files for each snapshot
     */
    if (verbosity > 0) writefln("lmr %s: Writing VTK files to disk.", cmdName);

    ensure_directory_is_present(lmrCfg.vtkDir);
    File pvdFile = begin_PVD_file(lmrCfg.vtkDir~"/"~lmrCfg.residualPrefix~".pvd");
    foreach (snap; snaps2process) {
        // We can't process snapshot 0000, no residual values computed
        if (snap == format(lmrCfg.snapshotIdxFmt, 0))
            continue;
        if (verbosity > 1) writefln("lmr %s: Writing snapshot %s to disk.", cmdName, snap);
        // We need to load a flow solution to get access to the grid and number of cells
        auto soln = new FlowSolution(to!int(snap), GlobalConfig.nFluidBlocks);
        string pvtuFileName = lmrCfg.residualPrefix ~ "-" ~ snap ~ ".pvtu";
        File pvtuFile = begin_PVTU_file(lmrCfg.vtkDir ~ "/" ~ pvtuFileName, variables);
        foreach (jb; 0 .. nBlocks) {
            readValuesFromFile(data, residualFilename(to!int(snap), jb), variables, soln.flowBlocks[jb].ncells, fileFmt);
            if (verbosity > 2) writefln("lmr %s: Writing block %d for snapshot %s to disk.", cmdName, jb, snap);
            string vtuFileName = lmrCfg.residualPrefix ~ "-" ~ format(lmrCfg.blkIdxFmt, jb) ~ "-" ~ snap ~ ".vtu";
            add_dataset_to_PVD_file(pvdFile, to!double(snap), vtuFileName);
            add_piece_to_PVTU_file(pvtuFile, vtuFileName);
            writeVTUfile(data, soln.gridBlocks[jb], variables, lmrCfg.vtkDir~"/"~vtuFileName, binaryFormat);
        }
        finish_PVTU_file(pvtuFile);
    }
    finish_PVD_file(pvdFile);

    if (verbosity > 0) writefln("lmr %s: Done.", cmdName);

    return 0;
}
