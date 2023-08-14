/**
 * Module for converting Eilmer-native limiter fields to VTK format.
 *
 * Authors: RJG, PJ, KAD, NNG
 * Date: 2023-08-14
 */

module limiter2vtk;

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

Command limiter2vtkCmd;
string cmdName = "limiter2vtk";

static this()
{
    limiter2vtkCmd.main = &main_;
    limiter2vtkCmd.description = "Convert fields of limiter values to VTK format.";
    limiter2vtkCmd.shortDescription = limiter2vtkCmd.description;
    limiter2vtkCmd.helpMsg = format(
`lmr %s [options]

Convert the limiter values for flow field using one or more snapshots to VTK format.

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

void main_(string[] args)
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
    getopt(args,
           config.bundling,
           "v|verbose+", &verbosity,
           "s|snapshots|snapshot", &snapshots,
           "f|final", &finalSnapshot,
           "a|all", &allSnapshots,
           "b|binary-format", &binaryFormat);

    if (verbosity > 0) writefln("lmr %s: Begin program.", cmdName);

    initConfiguration(); // To read in GlobalConfig
    nBlocks = GlobalConfig.nFluidBlocks;
    fileFmt = GlobalConfig.flow_format;
    variables = readVariablesFromMetadata(lmrCfg.limiterMetadataFile); 

    auto availSnapshots = determineAvailableSnapshots();
    auto snaps2process = determineSnapshotsToProcess(availSnapshots, snapshots, allSnapshots, finalSnapshot);

    /*
     * Now write vtk files for each snapshot
     */
    if (verbosity > 0) writefln("lmr %s: Writing VTK files to disk.", cmdName);

    ensure_directory_is_present(lmrCfg.vtkDir);
    File pvdFile = begin_PVD_file(lmrCfg.vtkDir~"/"~lmrCfg.limiterPrefix~".pvd");
    foreach (snap; snaps2process) {
        // We can't process snapshot 0000, no limiter values computed
        if (snap == format(lmrCfg.snapshotIdxFmt, 0))
            continue;
        if (verbosity > 1) writefln("lmr %s: Writing snapshot %s to disk.", cmdName, snap);
        // We need to load a flow solution to get access to the grid and number of cells
        auto soln = new FlowSolution(to!int(snap), GlobalConfig.nFluidBlocks);
        string pvtuFileName = lmrCfg.limiterPrefix ~ "-" ~ snap ~ ".pvtu";
        File pvtuFile = begin_PVTU_file(lmrCfg.vtkDir ~ "/" ~ pvtuFileName, variables);
        foreach (jb; 0 .. nBlocks) {
            readValuesFromFile(data, limiterFilename(to!int(snap), jb), variables, soln.flowBlocks[jb].ncells, fileFmt);
            if (verbosity > 2) writefln("lmr %s: Writing block %d for snapshot %s to disk.", cmdName, jb, snap);
            string vtuFileName = lmrCfg.limiterPrefix ~ "-" ~ format(lmrCfg.blkIdxFmt, jb) ~ "-" ~ snap ~ ".vtu";
            add_dataset_to_PVD_file(pvdFile, to!double(snap), vtuFileName);
            add_piece_to_PVTU_file(pvtuFile, vtuFileName);
            writeVTUfile(data, soln.gridBlocks[jb], variables, lmrCfg.vtkDir~"/"~vtuFileName, binaryFormat);
        }
        finish_PVTU_file(pvtuFile);
    }
    finish_PVD_file(pvdFile);

    if (verbosity > 0) writefln("lmr %s: Done.", cmdName);

    return;
}

void writeVTUfile(double[][] data, Grid grid, string[] variables, string fileName, bool binaryFormat)
{
    auto fp = File(fileName, "wb"); // We may be writing some binary data.
    ubyte[] binary_data_string;
    ubyte[] binary_data;
    int binary_data_offset = 0;
    bool two_D = (grid.dimensions == 2);
    size_t NumberOfPoints = grid.nvertices;
    if (data.length != grid.ncells) {
        string msg = text("Mismatch between grid and data grid.ncells=",
                          grid.ncells, " data.length=", data.length);
        throw new FlowSolverException(msg);
    }
    size_t NumberOfCells = data.length;
    fp.write("<VTKFile type=\"UnstructuredGrid\" byte_order=\"BigEndian\">\n");
    fp.write("<UnstructuredGrid>");
    fp.writef("<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
              NumberOfPoints, NumberOfCells);
    //
    fp.write("<Points>\n");
    fp.write(" <DataArray type=\"Float32\" NumberOfComponents=\"3\"");
    if (binaryFormat) {
        fp.writef(" format=\"appended\" offset=\"%d\">", binary_data_offset);
        binary_data.length=0;
    } else {
        fp.write(" format=\"ascii\">\n");
    }
    foreach (i; 0 .. grid.nvertices) {
        float x = uflowz(grid[i].x.re);
        float y = uflowz(grid[i].y.re);
        float z = uflowz(grid[i].z.re);
        if (binaryFormat) {
            binary_data ~= nativeToBigEndian(x);
            binary_data ~= nativeToBigEndian(y);
            binary_data ~= nativeToBigEndian(z);
        } else {
            fp.writef(" %.18e %.18e %.18e\n", x,y,z);
        }
    }
    fp.write(" </DataArray>\n");
    if (binaryFormat) {
        uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
        binary_data_string ~= nativeToBigEndian(binary_data_count);
        binary_data_string ~= binary_data;
        binary_data_offset += 4 + binary_data.length;
    }
    fp.write("</Points>\n");
    //
    fp.write("<Cells>\n");
    fp.write(" <DataArray type=\"Int32\" Name=\"connectivity\"");
    if (binaryFormat) {
        fp.writef(" format=\"appended\" offset=\"%d\">", binary_data_offset);
        binary_data.length = 0;
    } else {
        fp.write(" format=\"ascii\">\n");
    }
    foreach (i; 0 .. grid.ncells) {
        auto ids = grid.get_vtx_id_list_for_cell(i);
        if (binaryFormat) {
            foreach (id; ids) { binary_data ~= nativeToBigEndian(to!int32_t(id)); }
        } else {
            foreach (id; ids) { fp.writef(" %d", id); }
            fp.write("\n");
        }
    } // end foreach i
    fp.write(" </DataArray>\n");
    if (binaryFormat) {
        uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
        binary_data_string ~= nativeToBigEndian(binary_data_count);
        binary_data_string ~= binary_data;
        binary_data_offset += 4 + binary_data.length;
    }
    //
    fp.write(" <DataArray type=\"Int32\" Name=\"offsets\"");
    if (binaryFormat) {
        fp.writef(" format=\"appended\" offset=\"%d\">", binary_data_offset);
        binary_data.length = 0;
    } else {
        fp.write(" format=\"ascii\">\n");
    }
    // Since all of the point-lists are concatenated, these offsets into the connectivity
    // array specify the end of each cell.
    size_t conn_offset = 0;
    foreach (i; 0 .. grid.ncells) {
        conn_offset += grid.number_of_vertices_for_cell(i);
        if (binaryFormat) {
            binary_data ~= nativeToBigEndian(to!int32_t(conn_offset));
        } else {
            fp.writef(" %d\n", conn_offset);
        }
    } // end foreach i
    fp.write(" </DataArray>\n");
    if (binaryFormat) {
        uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
        binary_data_string ~= nativeToBigEndian(binary_data_count);
        binary_data_string ~= binary_data;
        binary_data_offset += 4 + binary_data.length;
    }
    //
    fp.write(" <DataArray type=\"UInt8\" Name=\"types\"");
    if (binaryFormat) {
        fp.writef(" format=\"appended\" offset=\"%d\">", binary_data_offset);
        binary_data.length = 0;
    } else {
        fp.write(" format=\"ascii\">\n");
    }
    foreach (i; 0 .. grid.ncells) {
        int type_value = grid.vtk_element_type_for_cell(i);
        if (binaryFormat) {
            binary_data ~= nativeToBigEndian(to!uint8_t(type_value));
        } else {
            fp.writef(" %d\n", type_value);
        }
    } // end foreach i
    fp.write(" </DataArray>\n");
    if (binaryFormat) {
        uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
        binary_data_string ~= nativeToBigEndian(binary_data_count);
        binary_data_string ~= binary_data;
        binary_data_offset += 4 + binary_data.length;
    }
    fp.write("</Cells>\n");
    //
    fp.write("<CellData>\n");
    // Write variables from the dictionary.
    foreach (ivar, var; variables) {
        fp.writef(" <DataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\"", var);
        if (binaryFormat) {
            fp.writef(" format=\"appended\" offset=\"%d\">", binary_data_offset);
            binary_data.length = 0;
        } else {
            fp.write(" format=\"ascii\">\n");
        }
        foreach (i; 0 .. grid.ncells) {
            if (binaryFormat) {
                binary_data ~= nativeToBigEndian(to!float(uflowz(data[i][ivar])));
            } else {
                fp.writef(" %.18e\n", uflowz(data[i][ivar]));
            }
        } // end foreach i
        fp.write(" </DataArray>\n");
        if (binaryFormat) {
            uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
            binary_data_string ~= nativeToBigEndian(binary_data_count);
            binary_data_string ~= binary_data;
            binary_data_offset += 4 + binary_data.length;
        }
    } // end foreach var
    fp.write("</CellData>\n");
    fp.write("</Piece>\n");
    fp.write("</UnstructuredGrid>\n");
    if (binaryFormat) {
        fp.write("<AppendedData encoding=\"raw\">\n");
        fp.write('_');
        fp.rawWrite(binary_data_string);
        fp.write("</AppendedData>\n");
    }
    fp.write("</VTKFile>\n");
    fp.close();
    return;
}
