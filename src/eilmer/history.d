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
import std.conv;

import globalconfig;
import globaldata;
import fvcell;
import fileutil;
import solidfvcell;
import solidblock;
import fluidblock;
import fluidblockio_old;

string histDir = "hist";

string hist_format_ext(string flow_format)
{
    switch(flow_format) {
        case "gziptext" : return "dat";
        case "rawbinary" : return "bin";
        case "eilmer4text" : return "dat";
        case "eilmer4binary" : return "bin";
        default : throw new Error("How did we get here?");
    }
}

void write_history_header_ascii(string fname){
    auto f = File(fname, "w");
    f.write("# 1:t ");
    foreach (i, var; GlobalConfig.flow_variable_list) {
        f.write(format("%d:%s ", i+2, var));
    }
    f.write("\n");
    f.close();
}

void write_history_header_rawbinary(string fname){
    int[1] integer; char[] chars;
    auto f = File(fname, "wb");

    integer[0] = to!int(GlobalConfig.flow_variable_list.length)+1;
    f.rawWrite(integer);

    integer[0] = 1;
    f.rawWrite(integer);

    chars.length = 1;
    chars[0] = 't';
    f.rawWrite(chars);

    foreach (var; GlobalConfig.flow_variable_list) {
        integer[0] = to!int(var.length);
        f.rawWrite(integer);

        chars.length = var.length;
        chars[] = var[];
        f.rawWrite(chars);
    }
    f.close();
}

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
        string flowFileExt = hist_format_ext(GlobalConfig.flow_format);
        string basicName = format("%s-blk-%d-cell-%d.%s", GlobalConfig.base_file_name, blkId, cellId, flowFileExt);
        auto foundTheseEntries = dirEntries(histDir, basicName~".*", SpanMode.shallow);
        string[] foundTheseNames;
        foreach (entry; foundTheseEntries) { foundTheseNames ~= entry.name; }
        string fname = format("%s/%s.%d", histDir, basicName, foundTheseNames.length);
        if (flowFileExt=="bin"){
            write_history_header_rawbinary(fname);
        } else {
            write_history_header_ascii(fname);
        }
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

void write_history_line_ascii(string fname, double sim_time, FVCell cell){
    auto writer = appender!string();
    formattedWrite(writer, "%.18e %s\n", sim_time, cell.write_values_to_string());
    append(fname, writer.data);
}


void write_history_line_rawbinary(string fname, double sim_time, FVCell cell){
    double[1] dbl;

    auto f = File(fname, "ab");
    dbl[0] = sim_time;
    f.rawWrite(dbl);
    cell.write_values_to_raw_binary(f);
    f.close();
}

void write_history_cells_to_files(double sim_time)
{
    foreach (hcell; GlobalConfig.hcells) {
        auto blkId = hcell[0];
        auto cellId = hcell[1];
        if (find(GlobalConfig.localFluidBlockIds, blkId).empty) { continue; }
        auto blk = cast(FluidBlock) globalBlocks[blkId];
        assert(blk !is null, "Oops, this should be a FluidBlock object.");
        auto cell = blk.cells[cellId];

        string flowFileExt = hist_format_ext(GlobalConfig.flow_format);
        string basicName = format("%s-blk-%d-cell-%d.%s", GlobalConfig.base_file_name, blkId, cellId, flowFileExt);
        auto foundTheseEntries = dirEntries(histDir, basicName~".*", SpanMode.shallow);
        string[] foundTheseNames;
        foreach (entry; foundTheseEntries) { foundTheseNames ~= entry.name; }
        string fname = format("%s/%s.%d", histDir, basicName, foundTheseNames.length-1);
        if (flowFileExt=="bin"){
            write_history_line_rawbinary(fname, sim_time, cell);
        } else {
            write_history_line_ascii(fname, sim_time, cell);
        }
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
