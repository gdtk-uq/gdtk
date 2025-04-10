/**
 * Module for converting Eilmer surface loads to VTK format.
 *
 * Authors: RBO
 * Date: 2025-04-07
 */

module lmr.commands.loads2vtk;

import std.algorithm;
import std.array;
import std.conv : to;
import std.file;
import std.format : format, formattedRead;
import std.getopt;
import std.range : empty;
import std.stdio;
import std.string;

import geom.grid : Grid_t, gridTypeFromName;
import geom.grid.sgrid : StructuredGrid;
import geom.grid.usgrid : UnstructuredGrid;

import lmr.commands.cmdhelper;
import lmr.commands.command;
import lmr.fileutil;
import lmr.globalconfig;
import lmr.init : initConfiguration;
import lmr.lmrconfig;
import lmr.ufluidblock;
import lmr.vtk_writer;

Command loads2vtkCmd;
string cmdName = "loads2vtk";

static this()
{
    loads2vtkCmd.main = &main_;
    loads2vtkCmd.description = "Convert surface loads to VTK format.";
    loads2vtkCmd.shortDescription = loads2vtkCmd.description;
    loads2vtkCmd.helpMsg = format(
`lmr %s [options]

Convert the surface loads files to VTK format for a three-dimensional grid.

options ([+] can be repeated):

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
    int blkId, bndryId, nBlocks;
    string[] splitHandle;
    UnstructuredGrid grid, bndryGrid;
    StructuredGrid sgrid;
    string gName, fileFmt, vtuFilename;
    size_t nCells;
    Grid_t gridType;

    int verbosity = 0;
    bool binaryFormat = false;
    try {
        getopt(args,
               config.bundling,
               "v|verbose+", &verbosity,
               "b|binary-format", &binaryFormat);
    } catch (Exception e) {
        writefln("Eilmer %s program quitting.", cmdName);
        writeln("There is something wrong with the command-line arguments/options.");
        writeln(e.msg);
        return 1;
    }

    if (verbosity > 0) writefln("lmr %s: Begin program.", cmdName);

    initConfiguration(); // To read in GlobalConfig
    fileFmt = GlobalConfig.field_format;
    nBlocks = GlobalConfig.nFluidBlocks;

    // ensure the grid is three-dimensional
    if (GlobalConfig.dimensions != 3) {
        string errMsg = format("The chosen grid has '%d' dimensions.\n", GlobalConfig.dimensions);
        errMsg ~= "The loads2vtk command only supports three-dimensional grids.";
        throw new Error(errMsg);
    }

    // give warning if grid motion is enabled
    if (GlobalConfig.grid_motion != GridMotion.none) {
        string msg = "Warning: grid motion is enabled but the loads2vtk command does not account for moving grids.\n";
        msg ~= "The loads will be attached to the surface grids from the initial snapshot directory.";
        writeln(msg);
    }

    // find last loads index
    // TODO: add option to process all loads indices
    string loadsMetadataFile = lmrCfg.loadsDir ~ "/" ~ lmrCfg.loadsPrefix ~ "-metadata";
    auto file = File(loadsMetadataFile);
    auto range = file.byLine();
    size_t lastLoadsIdx;
    foreach (line; range) {
        if (!line.empty && line[0] != 'l') {
            auto vars = line.split(" ");
            lastLoadsIdx = to!size_t(vars[0]);
        }
    }

    // generate list of loads files
    string loadsDir = lmrCfg.loadsDir ~ "/" ~ format(lmrCfg.snapshotIdxFmt, lastLoadsIdx);
    auto loadsFiles = dirEntries(loadsDir, SpanMode.shallow).
                 filter!(a => a.isFile()).map!(a => a.name).array();
    sort(loadsFiles);

    // get variable names from header in first loads file
    auto firstLoadsFile = File(loadsFiles[0]);
    auto loadsVarNames = firstLoadsFile.readln().chomp().split();
    size_t nVars = loadsVarNames.length;

    // get grid type for each block
    Grid_t[] gridTypes = new Grid_t[](nBlocks);
    string blkListFile = lmrCfg.blkListFile;
    auto listFile = File(blkListFile);
    auto listFileLine = listFile.readln().chomp(); // only comments on the first line
    foreach (ib; 0 .. nBlocks) {
        listFileLine = listFile.readln().chomp();
        int ib_chk;
        string gridTypeName;
        string label;
        int ncells;
        formattedRead(listFileLine, " %d %s %s %d", &ib_chk, &gridTypeName, &label, &ncells);
        if (ib != ib_chk) {
            string msg = format("Reading %s file ib=%d ib_chk=%d", blkListFile, ib, ib_chk);
            throw new FlowSolverException(msg);
        }
        gridTypes[blkId] = gridTypeFromName(gridTypeName);
    }

    // write VTK files
    if (verbosity > 0) writefln("lmr %s: Writing VTK files to disk.", cmdName);

    ensure_directory_is_present(lmrCfg.vtkDir);
    File pvdFile = begin_PVD_file(lmrCfg.vtkDir ~ "/" ~ lmrCfg.loadsPrefix ~ ".pvd");
    string pvtuFilename = lmrCfg.loadsPrefix ~ "-" ~ format(lmrCfg.snapshotIdxFmt, 0) ~ ".pvtu"; // TODO: update index
    File pvtuFile = begin_PVTU_file(lmrCfg.vtkDir ~ "/" ~ pvtuFilename, loadsVarNames);
    foreach (loadsFile; loadsFiles) {
        // get block and boundary ids
        splitHandle = loadsFile.split("-");
        blkId = to!int(splitHandle[1]);
        bndryId = to!int(splitHandle[3]);

        // get grid type for given blkId
        gridType = gridTypes[blkId];

        // construct boundary grid for given bndryId
        gName = gridFilename(lmrCfg.initialFieldDir, blkId);
        final switch (gridType) {
            case Grid_t.structured_grid:
                /*
                // TODO: add support for structured grids
                sgrid = new StructuredGrid(gName, fileFmt);
                grid = new UnstructuredGrid(sgrid);
                */
                string errMsg = format("Block %d has a structured grid type.\n", blkId);
                errMsg ~= "Currently, the loads2vtk command only supports unstructured grids.";
                throw new Error(errMsg);
                break;
            case Grid_t.unstructured_grid:
                grid = new UnstructuredGrid(gName, fileFmt, true);
                break;
        }
        bndryGrid = cast(UnstructuredGrid) grid.get_boundary_grid(bndryId);
        nCells = bndryGrid.ncells;

        // size data array
        data = new double[][](nCells, nVars);

        // read in data from loads file
        file = File(loadsFile);
        range = file.byLine();
        size_t offset = 0;
        foreach (line; range) {
            if (!line.empty && line[0] != 'p') {
                auto vars = line.split(" ");
                foreach (k, var; vars) {
                    data[offset][k] = to!double(var);
                }
                offset += 1;
            }
        }
        file.close();

        // write VTU file
        vtuFilename = lmrCfg.loadsPrefix 
            ~ "-blk-" 
            ~ format(lmrCfg.blkIdxFmt, blkId) 
            ~ "-bndry-" 
            ~ to!string(bndryId) 
            ~ ".vtu";
        add_dataset_to_PVD_file(pvdFile, 0, vtuFilename); // TODO: update index
        add_piece_to_PVTU_file(pvtuFile, vtuFilename);
        writeVTUfile(data, bndryGrid, loadsVarNames, lmrCfg.vtkDir ~ "/" ~ vtuFilename, binaryFormat);
    }
    finish_PVTU_file(pvtuFile);
    finish_PVD_file(pvdFile);

    if (verbosity > 0) writefln("lmr %s: Done.", cmdName);

    return 0;
}
