/**
 * Cell history writing functions.
 *
 * Author: Rowan G. and Peter J.
 * First code: 2016-01-31
 * History:
 *  2024-03-09 Updated for eilmer 5
 */

module lmr.history;

import std.stdio : File, writeln;
import std.format : formattedWrite;
import std.array : appender, split;
import std.algorithm: find;
import std.range : empty, walkLength;
import std.format : format;
import std.file : append, dirEntries, SpanMode, DirEntry;

import globalconfig : GlobalConfig;
import globaldata : globalBlocks;
import lmrconfig : lmrCfg, historyFilename;
import fileutil : ensure_directory_is_present;
import fluidblock : FluidBlock;
import lmrexceptions : TimeMarchingException;
import lmr.fluidfvcell : FluidFVCell;
import fvcellio : FVCellIO, FieldVarsType, createFVCellIO, buildFluidVariables;
version(mpi_parallel) {
    import mpi : MPI_Barrier, MPI_COMM_WORLD;
}

FVCellIO cio;

/**
 * Initialise directory and files for history cells.
 *
 * Authors: RJG and PAJ
 * Date: 2024-03-09
 */
void initHistoryCells()
{
    // Get non-shared copy for initialisation
    string[] varList = buildFluidVariables();
    cio = createFVCellIO(FieldVarsType.fluid, varList);
    if (GlobalConfig.is_master_task) ensure_directory_is_present(lmrCfg.historyDir);
    version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }

    foreach (ih, hcell; GlobalConfig.hcells) {
        auto blkId = hcell[0];
        auto cellId = hcell[1];
        auto blk = cast(FluidBlock) globalBlocks[blkId];
        assert(blk !is null, "Oops, this should be a FluidBlock object.");
        if (find(GlobalConfig.localFluidBlockIds, blkId).empty) { continue; }
        if ( cellId >= blk.cells.length ) {
            string errMsg = "ERROR: initHistoryCells()\n";
            errMsg ~= format("The requested history cell index %d is not valid for block %d.\n", cellId, blkId);
            errMsg ~= format("This block only has %d cells.\n", blk.cells.length);
            throw new TimeMarchingException(errMsg);
        }
        string basicName = historyFilename(ih, blkId, cellId);
        string nameWithoutDir = basicName.split("/")[2];
        auto numberOfExistingHistFiles = walkLength(dirEntries(lmrCfg.historyDir, nameWithoutDir~".*", SpanMode.shallow));
        string fname = format("%s.%d", basicName, numberOfExistingHistFiles);
        auto f = File(fname, "w");
        writeln("Creating file: ", fname);
        f.write("t");
        foreach (i, var; varList) {
            f.write(format(" %s", var));
        }
        f.write("\n");
        f.close();
    } // end foreach hcell

    // [TODO] RJG, 2024-02-07
    // Add something similar for solid domain cells.

}

/**
 * Write data for history cells to files.
 *
 * Authors: PAJ and RJG
 * Date: 2024-03-09
 */
void writeHistoryCellsToFiles(double sim_time)
{
    foreach (ih, hcell; GlobalConfig.hcells) {
        auto blkId = hcell[0];
        auto cellId = hcell[1];
        if (find(GlobalConfig.localFluidBlockIds, blkId).empty) { continue; }
        auto blk = cast(FluidBlock) globalBlocks[blkId];
        assert(blk !is null, "Oops, this should be a FluidBlock object.");
        auto cell = blk.cells[cellId];
        auto writer = appender!string();
        formattedWrite(writer, "%.18e %s\n", sim_time, cell.historyDataAsString());
        string basicName = historyFilename(ih, blkId, cellId);
        string nameWithoutDir = basicName.split("/")[2];
        auto numberOfExistingHistFiles = walkLength(dirEntries(lmrCfg.historyDir, nameWithoutDir~".*", SpanMode.shallow));
        string fname = format("%s.%d", basicName, numberOfExistingHistFiles-1);
        append(fname, writer.data);
    } // end foreach hcell
}

/**
 * History cell data as string.
 *
 * Authors: RJG and PAJ
 * Date: 2024-03-09
 */
string historyDataAsString(FluidFVCell cell)
{
    auto writer = appender!string();
    foreach (i, var; buildFluidVariables()) {
        formattedWrite(writer, "%.18e ", cio[cell, var]);
    }
    return writer.data;
}
