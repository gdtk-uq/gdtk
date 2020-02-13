/**
 * Eilmer4 cell history writing functions.
 *
 * Author: Rowan G. and Peter J.
 * First code: 2016-01-31
 */

module history;

import std.stdio;
import std.string;
import std.file;
import std.array;
import std.format;
import std.algorithm: find;

import fvcore;
import globalconfig;
import globaldata;
import fvcell;
import fileutil;
import solidfvcell;
import solidblock;
import fluidblock;

string histDir = "hist";

void init_history_cell_files()
{
    foreach (hcell; GlobalConfig.hcells) {
        auto blkId = hcell[0];
        auto cellId = hcell[1];
        auto blk = cast(FluidBlock) globalBlocks[blkId];
        assert(blk !is null, "Oops, this should be a FluidBlock object.");
        if (find(GlobalConfig.localFluidBlockIds, blkId).empty) { continue; }
        if ( cellId >= blk.cells.length ) {
            string errMsg = "ERROR: init_history_cell_files()\n";
            errMsg ~= format("The requested history cell index %d is not valid for block %d.\n", cellId, blkId);
            errMsg ~= format("This block only has %d cells.\n", blk.cells.length);
            throw new FlowSolverException(errMsg);
        }
        string basicName = format("%s-blk-%d-cell-%d.dat", GlobalConfig.base_file_name, blkId, cellId);
        auto foundTheseEntries = dirEntries(histDir, basicName~".*", SpanMode.shallow);
        string[] foundTheseNames;
        foreach (entry; foundTheseEntries) { foundTheseNames ~= entry.name; }
        string fname = format("%s/%s.%d", histDir, basicName, foundTheseNames.length);
        auto f = File(fname, "w");
        f.write("# 1:t ");
        foreach (i, var; GlobalConfig.flow_variable_list) {
            f.write(format("%d:%s ", i+2, var));
        }
        f.write("\n");
        f.close();
    } // end foreach hcell
    //
    // And history cells in solid domain, if present
    foreach (hcell; GlobalConfig.solid_hcells) {
        auto blkId = hcell[0];
        auto cellId = hcell[1];
        if (GlobalConfig.in_mpi_context) { throw new Error("[TODO] not available in MPI context."); }
        if (cellId >= localSolidBlocks[blkId].activeCells.length) {
            string errMsg = "ERROR: init_history_cell_files()\n";
            errMsg ~= format("The requested history cell index %d is not valid for solid block %d.\n", cellId, blkId);
            errMsg ~= format("This solid block only has %d cells.\n", localSolidBlocks[blkId].activeCells.length);
            throw new FlowSolverException(errMsg);
        }
        string fname = format("%s/%s-solid-blk-%d-cell-%d.dat", 
                              histDir, GlobalConfig.base_file_name, blkId, cellId);
        auto f = File(fname, "w");
        auto cellVars = varListForSolidCell();
        f.write("# 1:t ");
        foreach ( i, var; cellVars ) {
            f.write(format("%d:%s ", i+2, var));
        }
        f.write("\n");
        f.close();
    } // end foreach hcell
} // end init_history_cell_files()

void write_history_cells_to_files(double sim_time)
{
    foreach (hcell; GlobalConfig.hcells) {
        auto blkId = hcell[0];
        auto cellId = hcell[1];
        if (find(GlobalConfig.localFluidBlockIds, blkId).empty) { continue; }
        auto blk = cast(FluidBlock) globalBlocks[blkId];
        assert(blk !is null, "Oops, this should be a FluidBlock object.");
        auto cell = blk.cells[cellId];
        auto writer = appender!string();
        formattedWrite(writer, "%.18e %s\n", sim_time, cell.write_values_to_string());
        string basicName = format("%s-blk-%d-cell-%d.dat", GlobalConfig.base_file_name, blkId, cellId);
        auto foundTheseEntries = dirEntries(histDir, basicName~".*", SpanMode.shallow);
        string[] foundTheseNames;
        foreach (entry; foundTheseEntries) { foundTheseNames ~= entry.name; }
        string fname = format("%s/%s.%d", histDir, basicName, foundTheseNames.length-1);
        append(fname, writer.data);
    } // end foreach hcell
    // And history cells in solid domain, if present
    foreach (hcell; GlobalConfig.solid_hcells) {
        auto blkId = hcell[0];
        auto cellId = hcell[1];
        if (GlobalConfig.in_mpi_context) { throw new Error("[TODO] not available in MPI context."); }
        string fname = format("%s/%s-solid-blk-%d-cell-%d.dat", 
                              histDir, GlobalConfig.base_file_name, blkId, cellId);
        auto cell = localSolidBlocks[blkId].activeCells[cellId];
        auto writer = appender!string();
        formattedWrite(writer, "%.18e %s\n", sim_time,
                       cell.writeValuesToString());
        append(fname, writer.data);
    } // end foreach hcell
} // end write_history_cells_to_files()
