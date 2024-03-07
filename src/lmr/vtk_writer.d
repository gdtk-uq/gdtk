// vtk_writer.d
// Functions to write the flow solution data in VTK format.
//
// Author: Peter J. and Rowan G.
// 2021-01-05 Extracted from postprocess.d

module vtk_writer;

import std.math;
import std.stdio;
import std.conv;
import std.format;
import std.string;
import std.regex;
import std.algorithm;
import std.bitmanip;
import std.stdint;
import std.range;
import ntypes.complex;
import nm.number;
import gzip;
import fileutil;
import geom;
import gas;
import globalconfig;
import flowsolution : FluidBlockLite;
import solidsolution : SolidBlockLite;
import lmrexceptions : LmrPostProcessingException;


File begin_PVD_file(string fileName)
{
    // Start a Paraview collection file.
    // For each time index, this justs lists the name of the top-level .pvtu file.
    File f = File(fileName, "w");
    f.write("<?xml version=\"1.0\"?>\n");
    f.write("<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    f.write("<Collection>\n");
    return f;
}

void add_dataset_to_PVD_file(File f, double timeStamp, string vtuFileName)
{
    f.writef("<DataSet timestep=\"%.18e\" group=\"\" part=\"0\" file=\"%s\"/>\n",
             timeStamp, vtuFileName);
}

void finish_PVD_file(File f)
{
    f.write("</Collection>\n");
    f.write("</VTKFile>\n");
    f.close();
}

File begin_PVTU_file(string fileName, string[] variableNames)
{
    File f = File(fileName, "w");
    f.write("<VTKFile type=\"PUnstructuredGrid\">\n");
    f.write("<PUnstructuredGrid GhostLevel=\"0\">");
    f.write("<PPoints>\n");
    f.write(" <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n");
    f.write("</PPoints>\n");
    f.write("<PCellData>\n");
    foreach (var; variableNames) {
        f.writef(" <PDataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\"/>\n", var);
    }
    if (canFind(variableNames,"vel.x")) {
        f.write(" <PDataArray Name=\"vel.vector\" type=\"Float32\" NumberOfComponents=\"3\"/>\n");
    }
    if (canFind(variableNames,"c.x")) {
        f.write(" <PDataArray Name=\"c.vector\" type=\"Float32\" NumberOfComponents=\"3\"/>\n");
    }
    if (canFind(variableNames,"B.x")) {
        f.write(" <PDataArray Name=\"B.vector\" type=\"Float32\" NumberOfComponents=\"3\"/>\n");
    }
    f.write("</PCellData>\n");
    return f;
} // end begin_PVTU_file()

void add_piece_to_PVTU_file(File f, string fileName)
{
    f.writef("<Piece Source=\"%s\"/>\n", fileName);
}

void finish_PVTU_file(File f)
{
    f.write("</PUnstructuredGrid>\n");
    f.write("</VTKFile>\n");
    f.close();
}

void write_VTU_file(FluidBlockLite flow, Grid grid, string fileName, bool binary_format)
// Write the cell-centred flow data from a single block (index jb)
// as an unstructured grid of finite-volume cells.
{
    auto fp = File(fileName, "wb"); // We may be writing some binary data.
    ubyte[] binary_data_string;
    ubyte[] binary_data;
    int binary_data_offset = 0;
    bool two_D = (grid.dimensions == 2);
    size_t NumberOfPoints = grid.nvertices;
    if (flow.ncells != grid.ncells) {
        string msg = text("Mismatch between grid and flow grid.ncells=",
                          grid.ncells, " flow.ncells=", flow.ncells);
        throw new LmrPostProcessingException(msg);
    }
    size_t NumberOfCells = flow.ncells;
    fp.write("<VTKFile type=\"UnstructuredGrid\" byte_order=\"BigEndian\">\n");
    fp.write("<UnstructuredGrid>");
    fp.writef("<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
              NumberOfPoints, NumberOfCells);
    //
    fp.write("<Points>\n");
    fp.write(" <DataArray type=\"Float32\" NumberOfComponents=\"3\"");
    if (binary_format) {
        fp.writef(" format=\"appended\" offset=\"%d\">", binary_data_offset);
        binary_data.length=0;
    } else {
        fp.write(" format=\"ascii\">\n");
    }
    foreach (i; 0 .. grid.nvertices) {
        float x = uflowz(grid[i].x.re);
        float y = uflowz(grid[i].y.re);
        float z = uflowz(grid[i].z.re);
        if (binary_format) {
            binary_data ~= nativeToBigEndian(x);
            binary_data ~= nativeToBigEndian(y);
            binary_data ~= nativeToBigEndian(z);
        } else {
            fp.writef(" %.18e %.18e %.18e\n", x,y,z);
        }
    }
    fp.write(" </DataArray>\n");
    if (binary_format) {
        uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
        binary_data_string ~= nativeToBigEndian(binary_data_count);
        binary_data_string ~= binary_data;
        binary_data_offset += 4 + binary_data.length;
    }
    fp.write("</Points>\n");
    //
    fp.write("<Cells>\n");
    fp.write(" <DataArray type=\"Int32\" Name=\"connectivity\"");
    if (binary_format) {
        fp.writef(" format=\"appended\" offset=\"%d\">", binary_data_offset);
        binary_data.length = 0;
    } else {
        fp.write(" format=\"ascii\">\n");
    }
    foreach (i; 0 .. grid.ncells) {
        auto ids = grid.get_vtx_id_list_for_cell(i);
        if (binary_format) {
            foreach (id; ids) { binary_data ~= nativeToBigEndian(to!int32_t(id)); }
        } else {
            foreach (id; ids) { fp.writef(" %d", id); }
            fp.write("\n");
        }
    } // end foreach i
    fp.write(" </DataArray>\n");
    if (binary_format) {
        uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
        binary_data_string ~= nativeToBigEndian(binary_data_count);
        binary_data_string ~= binary_data;
        binary_data_offset += 4 + binary_data.length;
    }
    //
    fp.write(" <DataArray type=\"Int32\" Name=\"offsets\"");
    if (binary_format) {
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
        if (binary_format) {
            binary_data ~= nativeToBigEndian(to!int32_t(conn_offset));
        } else {
            fp.writef(" %d\n", conn_offset);
        }
    } // end foreach i
    fp.write(" </DataArray>\n");
    if (binary_format) {
        uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
        binary_data_string ~= nativeToBigEndian(binary_data_count);
        binary_data_string ~= binary_data;
        binary_data_offset += 4 + binary_data.length;
    }
    //
    fp.write(" <DataArray type=\"UInt8\" Name=\"types\"");
    if (binary_format) {
        fp.writef(" format=\"appended\" offset=\"%d\">", binary_data_offset);
        binary_data.length = 0;
    } else {
        fp.write(" format=\"ascii\">\n");
    }
    foreach (i; 0 .. grid.ncells) {
        int type_value = grid.vtk_element_type_for_cell(i);
        if (binary_format) {
            binary_data ~= nativeToBigEndian(to!uint8_t(type_value));
        } else {
            fp.writef(" %d\n", type_value);
        }
    } // end foreach i
    fp.write(" </DataArray>\n");
    if (binary_format) {
        uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
        binary_data_string ~= nativeToBigEndian(binary_data_count);
        binary_data_string ~= binary_data;
        binary_data_offset += 4 + binary_data.length;
    }
    fp.write("</Cells>\n");
    //
    fp.write("<CellData>\n");
    // Write variables from the dictionary.
    foreach (var; flow.variableNames) {
        fp.writef(" <DataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\"", var);
        if (binary_format) {
            fp.writef(" format=\"appended\" offset=\"%d\">", binary_data_offset);
            binary_data.length = 0;
        } else {
            fp.write(" format=\"ascii\">\n");
        }
        foreach (i; 0 .. flow.ncells) {
            if (binary_format) {
                binary_data ~= nativeToBigEndian(to!float(uflowz(flow[var,i])));
            } else {
                fp.writef(" %.18e\n", uflowz(flow[var,i]));
            }
        } // end foreach i
        fp.write(" </DataArray>\n");
        if (binary_format) {
            uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
            binary_data_string ~= nativeToBigEndian(binary_data_count);
            binary_data_string ~= binary_data;
            binary_data_offset += 4 + binary_data.length;
        }
    } // end foreach var
    //
    // Write the special variables:
    // i.e. variables constructed from those in the dictionary.
    if (canFind(flow.variableNames, "vel.x") && canFind(flow.variableNames, "vel.y")) {
        fp.write(" <DataArray Name=\"vel.vector\" type=\"Float32\" NumberOfComponents=\"3\"");
        if (binary_format) {
            fp.writef(" format=\"appended\" offset=\"%d\">", binary_data_offset);
            binary_data.length = 0;
        } else {
            fp.write(" format=\"ascii\">\n");
        }

        bool isThreeDimensional = canFind(flow.variableNames, "vel.z");
        foreach (i; 0 .. flow.ncells) {
            float x = uflowz(flow["vel.x",i]);
            float y = uflowz(flow["vel.y",i]);
            float z = (isThreeDimensional) ? uflowz(flow["vel.z",i]) : 0.0;
            if (binary_format) {
                binary_data ~= nativeToBigEndian(x);
                binary_data ~= nativeToBigEndian(y);
                binary_data ~= nativeToBigEndian(z);
            } else {
                fp.writef(" %.18e %.18e %.18e\n", x, y, z);
            }
        } // end foreach i
        fp.write(" </DataArray>\n");
        if (binary_format) {
            uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
            binary_data_string ~= nativeToBigEndian(binary_data_count);
            binary_data_string ~= binary_data;
            binary_data_offset += 4 + binary_data.length;
        }
    }
    //
    if (canFind(flow.variableNames, "c.x")) {
        fp.write(" <DataArray Name=\"c.vector\" type=\"Float32\" NumberOfComponents=\"3\"");
        if (binary_format) {
            fp.writef(" format=\"appended\" offset=\"%d\">", binary_data_offset);
            binary_data.length = 0;
        } else {
            fp.write(" format=\"ascii\">\n");
        }
        foreach (i; 0 .. flow.ncells) {
            float x = uflowz(flow["c.x",i]);
            float y = uflowz(flow["c.y",i]);
            float z = uflowz(flow["c.z",i]);
            if (binary_format) {
                binary_data ~= nativeToBigEndian(x);
                binary_data ~= nativeToBigEndian(y);
                binary_data ~= nativeToBigEndian(z);
            } else {
                fp.writef(" %.18e %.18e %.18e\n", x, y, z);
            }
        } // end foreach i
        fp.write(" </DataArray>\n");
        if (binary_format) {
            uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
            binary_data_string ~= nativeToBigEndian(binary_data_count);
            binary_data_string ~= binary_data;
            binary_data_offset += 4 + binary_data.length;
        }
    } // if canFind c.x
    //
    if (canFind(flow.variableNames, "B.x")) {
        fp.write(" <DataArray Name=\"B.vector\" type=\"Float32\" NumberOfComponents=\"3\"");
        if (binary_format) {
            fp.writef(" format=\"appended\" offset=\"%d\">", binary_data_offset);
            binary_data.length = 0;
        } else {
            fp.write(" format=\"ascii\">\n");
        }
        foreach (i; 0 .. flow.ncells) {
            float x = uflowz(flow["B.x",i]);
            float y = uflowz(flow["B.y",i]);
            float z = uflowz(flow["B.z",i]);
            if (binary_format) {
                binary_data ~= nativeToBigEndian(x);
                binary_data ~= nativeToBigEndian(y);
                binary_data ~= nativeToBigEndian(z);
            } else {
                fp.writef(" %.18e %.18e %.18e\n", x, y, z);
            }
        } // end foreach i
        fp.write(" </DataArray>\n");
        if (binary_format) {
            uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
            binary_data_string ~= nativeToBigEndian(binary_data_count);
            binary_data_string ~= binary_data;
            binary_data_offset += 4 + binary_data.length;
        }
    } // if canFind B.x
    //
    fp.write("</CellData>\n");
    fp.write("</Piece>\n");
    fp.write("</UnstructuredGrid>\n");
    if (binary_format) {
        fp.write("<AppendedData encoding=\"raw\">\n");
        fp.write('_');
        fp.rawWrite(binary_data_string);
        fp.write("</AppendedData>\n");
    }
    fp.write("</VTKFile>\n");
    fp.close();
    return;
} // end write_VTU_file()


// This version is for the solid domain.
void write_VTU_file(SolidBlockLite solid, StructuredGrid grid, string fileName, bool binary_format)
// Write the cell-centred flow data from a single block (index jb)
// as an unstructured grid of finite-volume cells.
{
    auto fp = File(fileName, "wb"); // We may be writing some binary data.
    //auto solid = soln.solidBlocks[jb];
    //auto grid = soln.gridBlocks[jb];
    ubyte[] binary_data_string;
    ubyte[] binary_data;
    int binary_data_offset = 0;
    size_t niv = grid.niv; size_t njv = grid.njv; size_t nkv = grid.nkv;
    size_t nic = solid.nic; size_t njc = solid.njc; size_t nkc = solid.nkc;
    bool two_D = (nkv == 1);
    size_t NumberOfPoints = niv * njv * nkv;
    size_t NumberOfCells = nic * njc * nkc;
    fp.write("<VTKFile type=\"UnstructuredGrid\" byte_order=\"BigEndian\">\n");
    fp.write("<UnstructuredGrid>");
    fp.writef("<Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
              NumberOfPoints, NumberOfCells);
    //
    fp.write("<Points>\n");
    fp.write(" <DataArray type=\"Float32\" NumberOfComponents=\"3\"");
    if (binary_format) {
        fp.writef(" format=\"appended\" offset=\"%d\">", binary_data_offset);
        binary_data.length=0;
    } else {
        fp.write(" format=\"ascii\">\n");
    }
    size_t vtx_number = 0;
    size_t[][][] vtx_id;
    vtx_id.length = niv;
    foreach (i; 0 .. niv) {
        vtx_id[i].length = njv;
        foreach (j; 0 .. njv) {
            vtx_id[i][j].length = nkv;
        }
    }
    foreach (k; 0 .. nkv) {
        foreach (j; 0 .. njv) {
            foreach (i; 0 .. niv) {
                vtx_id[i][j][k] = vtx_number;
                float x = uflowz(grid[i,j,k].x.re);
                float y = uflowz(grid[i,j,k].y.re);
                float z = uflowz(grid[i,j,k].z.re);
                if (binary_format) {
                    binary_data ~= nativeToBigEndian(x);
                    binary_data ~= nativeToBigEndian(y);
                    binary_data ~= nativeToBigEndian(z);
                } else {
                    fp.writef(" %.18e %.18e %.18e\n", x,y,z);
                }
                vtx_number += 1;
            }
        }
    }
    fp.write(" </DataArray>\n");
    if (binary_format) {
        uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
        binary_data_string ~= nativeToBigEndian(binary_data_count);
        binary_data_string ~= binary_data;
        binary_data_offset += 4 + binary_data.length;
    }
    fp.write("</Points>\n");
    //
    fp.write("<Cells>\n");
    fp.write(" <DataArray type=\"Int32\" Name=\"connectivity\"");
    if (binary_format) {
        fp.writef(" format=\"appended\" offset=\"%d\">", binary_data_offset);
        binary_data.length = 0;
    } else {
        fp.write(" format=\"ascii\">\n");
    }
    foreach (k; 0 .. nkc) {
        foreach (j; 0 .. njc) {
            foreach (i; 0 .. nic) {
                if (two_D) {
                    auto ids = [vtx_id[i][j][k], vtx_id[i+1][j][k],
                                vtx_id[i+1][j+1][k], vtx_id[i][j+1][k]];
                    if (binary_format) {
                        foreach (id; ids) { binary_data ~= nativeToBigEndian(to!int32_t(id)); }
                    } else {
                        fp.writef(" %d %d %d %d\n", ids[0], ids[1], ids[2], ids[3]);
                    }
                } else {
                    auto ids = [vtx_id[i][j][k], vtx_id[i+1][j][k],
                                vtx_id[i+1][j+1][k], vtx_id[i][j+1][k],
                                vtx_id[i][j][k+1], vtx_id[i+1][j][k+1],
                                vtx_id[i+1][j+1][k+1], vtx_id[i][j+1][k+1]];
                    if (binary_format) {
                        foreach (id; ids) { binary_data ~= nativeToBigEndian(to!int32_t(id)); }
                    } else {
                        fp.writef(" %d %d %d %d %d %d %d %d\n", ids[0], ids[1], ids[2],
                                  ids[3], ids[4], ids[5], ids[6], ids[7]);
                    }
                }
            } // end foreach i
        } // end foreach j
    } // end foreach k
    fp.write(" </DataArray>\n");
    if (binary_format) {
        uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
        binary_data_string ~= nativeToBigEndian(binary_data_count);
        binary_data_string ~= binary_data;
        binary_data_offset += 4 + binary_data.length;
    }
    //
    fp.write(" <DataArray type=\"Int32\" Name=\"offsets\"");
    if (binary_format) {
        fp.writef(" format=\"appended\" offset=\"%d\">", binary_data_offset);
        binary_data.length = 0;
    } else {
        fp.write(" format=\"ascii\">\n");
    }
    // Since all of the point-lists are concatenated, these offsets into the connectivity
    // array specify the end of each cell.
    foreach (k; 0 .. nkc) {
        foreach (j; 0 .. njc) {
            foreach (i; 0 .. nic) {
                size_t conn_offset;
                if (two_D) {
                    conn_offset = 4*(1+i+j*nic);
                } else {
                    conn_offset = 8*(1+i+j*nic+k*(nic*njc));
                }
                if (binary_format) {
                    binary_data ~= nativeToBigEndian(to!int32_t(conn_offset));
                } else {
                    fp.writef(" %d\n", conn_offset);
                }
            } // end foreach i
        } // end foreach j
    } // end foreach k
    fp.write(" </DataArray>\n");
    if (binary_format) {
        uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
        binary_data_string ~= nativeToBigEndian(binary_data_count);
        binary_data_string ~= binary_data;
        binary_data_offset += 4 + binary_data.length;
    }
    //
    fp.write(" <DataArray type=\"UInt8\" Name=\"types\"");
    if (binary_format) {
        fp.writef(" format=\"appended\" offset=\"%d\">", binary_data_offset);
        binary_data.length = 0;
    } else {
        fp.write(" format=\"ascii\">\n");
    }
    int type_value;
    if (two_D) {
        type_value = 9; // VTK_QUAD
    } else {
        type_value = 12; // VTK_HEXAHEDRON
    }
    foreach (k; 0 .. nkc) {
        foreach (j; 0 .. njc) {
            foreach (i; 0 .. nic) {
                if (binary_format) {
                    binary_data ~= nativeToBigEndian(to!uint8_t(type_value));
                } else {
                    fp.writef(" %d\n", type_value);
                }
            } // end foreach i
        } // end foreach j
    } // end foreach k
    fp.write(" </DataArray>\n");
    if (binary_format) {
        uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
        binary_data_string ~= nativeToBigEndian(binary_data_count);
        binary_data_string ~= binary_data;
        binary_data_offset += 4 + binary_data.length;
    }
    fp.write("</Cells>\n");
    //
    fp.write("<CellData>\n");
    // Write variables from the dictionary.
    foreach (var; solid.variableNames) {
        fp.writef(" <DataArray Name=\"%s\" type=\"Float32\" NumberOfComponents=\"1\"", var);
        if (binary_format) {
            fp.writef(" format=\"appended\" offset=\"%d\">", binary_data_offset);
            binary_data.length = 0;
        } else {
            fp.write(" format=\"ascii\">\n");
        }
        foreach (k; 0 .. nkc) {
            foreach (j; 0 .. njc) {
                foreach (i; 0 .. nic) {
                    if (binary_format) {
                        binary_data ~= nativeToBigEndian(to!float(uflowz(solid[var,i,j,k])));
                    } else {
                        fp.writef(" %.18e\n", uflowz(solid[var,i,j,k]));
                    }
                } // end foreach i
            } // end foreach j
        } // end foreach k
        fp.write(" </DataArray>\n");
        if (binary_format) {
            uint32_t binary_data_count = to!uint32_t(binary_data.length); // 4-byte count of bytes
            binary_data_string ~= nativeToBigEndian(binary_data_count);
            binary_data_string ~= binary_data;
            binary_data_offset += 4 + binary_data.length;
        }
    } // end foreach var
    //
    fp.write("</CellData>\n");
    fp.write("</Piece>\n");
    fp.write("</UnstructuredGrid>\n");
    if (binary_format) {
        fp.write("<AppendedData encoding=\"raw\">\n");
        fp.write('_');
        fp.rawWrite(binary_data_string);
        fp.write("</AppendedData>\n");
    }
    fp.write("</VTKFile>\n");
    fp.close();
    return;
} // end write_VTU_file()

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
