// fluidblockio.d
// Old Input/Output functions for the flow solution data.
//
// 2021-Mar-12: PJ
// Gather existing I/O functions into this module and
// do the minimum to keep them working with the current code.
//
// These functions had been built in a piece-meal fashion
// while the Eilmer4 code evolved over the few years and
// thinking about the data file format was mainly at the cell-level.
// In Daryl's new IO design, the flow-data input and output
// is principally at the block-level.  Quite a change.
//
// This module gathers together the old input-output functions so that
// it is easier to preserve them while making room for the new IO code.
//

module fluidblockio_old;

import fvcell;

import std.conv;
import std.string;
import std.regex;
import std.array;
import std.format;
import std.stdio;
import std.math;
import std.algorithm;
import gzip;
import util.lua;
import util.lua_service;
import nm.complex;
import nm.number;
import nm.bbla;
import geom;
import geom.luawrap;
import gas;
import globalconfig;
import fvcore;
import flowstate;
import luaflowstate;
import turbulence;
import fluidblock;
import sfluidblock;
import ufluidblock;
import flowsolution;

// ----------------------------------------------------------------------
// Flow variables, as they might appear in the data files.
// We define their names here with a single, definitive spelling.
// Although having symbolic names for the variables might seem excessive,
// we hope to gain some benefit from the compiler being able to check them.
enum FlowVar {
    pos_x, pos_y, pos_z, volume,
    rho, vel_x, vel_y, vel_z,
    B_x, B_y, B_z, divB,
    psi, quality, p, a, mu, k,
    mu_t, k_t, S,
    Q_rad_org, f_rad_org, Q_rE_rad,
    dt_local, dt_chem, u, T, dt_therm
};

@nogc
string flowVarName(FlowVar var)
{
    final switch(var) {
    case FlowVar.pos_x: return "pos.x";
    case FlowVar.pos_y: return "pos.y";
    case FlowVar.pos_z: return "pos.z";
    case FlowVar.volume: return "volume";
    case FlowVar.rho: return "rho";
    case FlowVar.vel_x: return "vel.x";
    case FlowVar.vel_y: return "vel.y";
    case FlowVar.vel_z: return "vel.z";
    case FlowVar.B_x: return "B.x";
    case FlowVar.B_y: return "B.y";
    case FlowVar.B_z: return "B.z";
    case FlowVar.divB: return "divB";
    case FlowVar.psi: return "psi";
    case FlowVar.quality: return "quality";
    case FlowVar.p: return "p";
    case FlowVar.a: return "a";
    case FlowVar.mu: return "mu";
    case FlowVar.k: return "k";
    case FlowVar.mu_t: return "mu_t";
    case FlowVar.k_t: return "k_t";
    case FlowVar.S: return "S";
    case FlowVar.Q_rad_org: return "Q_rad_org";
    case FlowVar.f_rad_org: return "f_rad_org";
    case FlowVar.Q_rE_rad: return "Q_rE_rad";
    case FlowVar.dt_local: return "dt_local";
    case FlowVar.dt_chem: return "dt_chem";
    case FlowVar.u: return "u";
    case FlowVar.T: return "T";
    case FlowVar.dt_therm: return "dt_therm";
    } // end switch(var)
} // end FlowVarName
string massfName(GasModel gmodel, int i) {
    auto name = cast(char[]) gmodel.species_name(i);
    name = tr(name, " \t", "--", "s"); // Replace internal whitespace with dashes.
    return "massf[" ~ to!string(i) ~ "]-" ~ to!string(name);
}
string k_modesName(int i) { return "k_modes[" ~ to!string(i) ~ "]"; }
string u_modesName(int i) { return "u_modes[" ~ to!string(i) ~ "]"; }
string T_modesName(int i) { return "T_modes[" ~ to!string(i) ~ "]"; }

string[] build_flow_variable_list()
{
    // Returns a list of variable names in the order of the fixed-layout data files.
    // This function needs to be consistent with the cell-data reading and writing
    // functions toward the end of fvcell.d.
    string[] list;
    list ~= [flowVarName(FlowVar.pos_x),
             flowVarName(FlowVar.pos_y),
             flowVarName(FlowVar.pos_z),
             flowVarName(FlowVar.volume),
             flowVarName(FlowVar.rho),
             flowVarName(FlowVar.vel_x),
             flowVarName(FlowVar.vel_y),
             flowVarName(FlowVar.vel_z)];
    if (GlobalConfig.MHD) {
        list ~= [flowVarName(FlowVar.B_x),
                 flowVarName(FlowVar.B_y),
                 flowVarName(FlowVar.B_z),
                 flowVarName(FlowVar.divB)];
    }
    if (GlobalConfig.MHD && GlobalConfig.divergence_cleaning) { list ~= flowVarName(FlowVar.psi); }
    if (GlobalConfig.include_quality) { list ~= flowVarName(FlowVar.quality); }
    list ~= [flowVarName(FlowVar.p),
             flowVarName(FlowVar.a),
             flowVarName(FlowVar.mu),
             flowVarName(FlowVar.k)];
    foreach(i; 0 .. GlobalConfig.gmodel_master.n_modes) { list ~= k_modesName(i); }
    list ~= [flowVarName(FlowVar.mu_t),
             flowVarName(FlowVar.k_t),
             flowVarName(FlowVar.S)];
    if (GlobalConfig.radiation) {
        list ~= [flowVarName(FlowVar.Q_rad_org),
                 flowVarName(FlowVar.f_rad_org),
                 flowVarName(FlowVar.Q_rE_rad)];
    }
    foreach(i; 0 .. GlobalConfig.turb_model.nturb) list ~= GlobalConfig.turb_model.primitive_variable_name(i);
    foreach(i; 0 .. GlobalConfig.gmodel_master.n_species) { list ~= [massfName(GlobalConfig.gmodel_master, i)]; }
    if (GlobalConfig.gmodel_master.n_species > 1) { list ~= flowVarName(FlowVar.dt_chem); }
    list ~= [flowVarName(FlowVar.u), flowVarName(FlowVar.T)];
    foreach(i; 0 .. GlobalConfig.gmodel_master.n_modes) {
        list ~= [u_modesName(i), T_modesName(i)];
    }
    if (GlobalConfig.gmodel_master.n_modes > 0) { list ~= flowVarName(FlowVar.dt_therm); }
    if (GlobalConfig.with_local_time_stepping) { list ~= flowVarName(FlowVar.dt_local); }
    return list;
} // end build_flow_variable_list()


//--------------------------------- Block IO ----------------------------------------
//
double read_solution(FluidBlock blk, string filename, bool overwrite_geometry_data)
{
    auto sblk = cast(SFluidBlock) blk;
    if (sblk) { return read_solution(sblk, filename, overwrite_geometry_data); }
    auto ublk = cast(UFluidBlock) blk;
    if (ublk) { return read_solution(ublk, filename, overwrite_geometry_data); }
    throw new Error("Oops, unknown type of FluidBlock");
}

double read_solution(SFluidBlock blk, string filename, bool overwrite_geometry_data)
// Note that the position data is read into grid-time-level 0
// by scan_values_from_string().
// Returns sim_time from file.
// Keep in sync with write_initial_flow_file() in flowstate.d
// and write_solution below.
{
    if (blk.myConfig.verbosity_level > 1) { writeln("read_solution(): Start block ", blk.id); }
    double sim_time; // to be read from file
    string myLabel;
    int nvariables;
    string[] variableNames;
    switch (blk.myConfig.flow_format) {
    case "gziptext": goto default;
    case "rawbinary":
        File fin = File(filename, "rb");
        string expected_header = "structured_grid_flow 1.0";
        char[] found_header = new char[expected_header.length];
        fin.rawRead(found_header);
        if (found_header != expected_header) {
            throw new FlowSolverException("SFluidBlock.read_solution from raw_binary_file: " ~
                                          "unexpected header: " ~ to!string(found_header));
        }
        int[1] int1; fin.rawRead(int1);
        int label_length = int1[0];
        if (label_length > 0) {
            char[] found_label = new char[label_length];
            fin.rawRead(found_label);
            myLabel = to!string(found_label);
        }
        double[1] dbl1; fin.rawRead(dbl1); sim_time = dbl1[0];
        fin.rawRead(int1); nvariables = int1[0];
        foreach(i; 0 .. nvariables) {
            char[] varname; fin.rawRead(int1); varname.length = int1[0];
            fin.rawRead(varname);
        }
        int[4] int4; fin.rawRead(int4);
        int dimensions_read = int4[0];
        size_t nic_read = int4[1]; size_t njc_read = int4[2]; size_t nkc_read = int4[3];
        if (dimensions_read != blk.myConfig.dimensions) {
            string msg = text("dimensions found: " ~ to!string(dimensions_read));
            throw new FlowSolverException(msg);
        }
        if (blk.nic != nic_read || blk.njc != njc_read ||
            nkc_read != ((blk.myConfig.dimensions == 3) ? blk.nkc : 1)) {
            string msg = text("For block[", blk.id, "] we have a mismatch in solution size.",
                              " Have read nic=", nic_read, " njc=", njc_read, " nkc=", nkc_read);
            throw new FlowSolverException(msg);
        }
        foreach (k; 0 .. blk.nkc) {
            foreach (j; 0 .. blk.njc) {
                foreach (i; 0 .. blk.nic) {
                    blk.get_cell(i,j,k).read_values_from_raw_binary(fin, overwrite_geometry_data);
                }
            }
        }
        break;
    default:
        auto byLine = new GzipByLine(filename);
        auto line = byLine.front; byLine.popFront();
        string format_version;
        formattedRead(line, "structured_grid_flow %s", &format_version);
        if (format_version != "1.0") {
            string msg = text("File format version found: " ~ format_version);
            throw new FlowSolverException(msg);
        }
        line = byLine.front; byLine.popFront();
        formattedRead(line, "label: %s", &myLabel);
        line = byLine.front; byLine.popFront();
        formattedRead(line, "sim_time: %g", &sim_time);
        line = byLine.front; byLine.popFront();
        formattedRead(line, "variables: %d", &nvariables);
        line = byLine.front; byLine.popFront();
        variableNames = line.strip().split();
        foreach (i; 0 .. variableNames.length) { variableNames[i] = variableNames[i].strip("\""); }
        // If all of the variables line up,
        // it is probably faster to use fixed-order, so default to that.
        bool useFixedOrder = true;
        foreach (i, varName; blk.myConfig.flow_variable_list) {
            if (!std.algorithm.canFind(variableNames, varName)) {
                throw new Exception("Could not find variable: " ~ varName);
            }
            if (!std.algorithm.equal(variableNames[i], varName)) {
                useFixedOrder = false;
            }
        }
        line = byLine.front; byLine.popFront();
        size_t dimensions_read, nic_read, njc_read, nkc_read;
        formattedRead(line, "dimensions: %d", &dimensions_read);
        if (dimensions_read != blk.myConfig.dimensions) {
            string msg = text("dimensions found: " ~ to!string(dimensions_read));
            throw new FlowSolverException(msg);
        }
        line = byLine.front; byLine.popFront();
        formattedRead(line, "nicell: %d", &nic_read);
        line = byLine.front; byLine.popFront();
        formattedRead(line, "njcell: %d", &njc_read);
        line = byLine.front; byLine.popFront();
        formattedRead(line, "nkcell: %d", &nkc_read);
        if (nic_read != blk.nic || njc_read != blk.njc ||
            nkc_read != ((blk.myConfig.dimensions == 3) ? blk.nkc : 1)) {
            string msg = text("For block[", blk.id, "] we have a mismatch in solution size.",
                              " Have read nic=", nic_read, " njc=", njc_read, " nkc=", nkc_read);
            throw new FlowSolverException(msg);
        }
        foreach (k; 0 .. blk.nkc) {
            foreach (j; 0 .. blk.njc) {
                foreach (i; 0 .. blk.nic) {
                    line = byLine.front; byLine.popFront();
                    blk.get_cell(i,j,k).scan_values_from_string(line, variableNames, useFixedOrder,
                                                                blk.myConfig.gmodel, overwrite_geometry_data);
                }
            }
        }
    } // end switch flow_format
    if (blk.myConfig.verbosity_level > 1) { writeln("read_solution(): Done block ", blk.id); }
    return sim_time;
} // end read_solution() for structured-grid block

double read_solution(UFluidBlock blk, string filename, bool overwrite_geometry_data)
// Note that this function needs to be kept in sync with the BlockFlow class
// over in flowsolution.d and with write_solution() below and with
// write_initial_usg_flow_file_from_lua() in luaflowstate.d.
// Returns sim_time from file.
{
    if (blk.myConfig.verbosity_level > 1) { writeln("read_solution(): Start block ", blk.id); }
    double sim_time; // to be read from file
    string myLabel;
    size_t nvariables;
    string[] variableNames;
    int my_dimensions;
    size_t nc;
    switch (blk.myConfig.flow_format) {
    case "gziptext": goto default;
    case "rawbinary":
        File fin = File(filename, "rb");
        string expected_header = "unstructured_grid_flow 1.0";
        char[] found_header = new char[expected_header.length];
        fin.rawRead(found_header);
        if (found_header != expected_header) {
            throw new FlowSolverException("UFluidBlock.read_solution from raw_binary_file: " ~
                                          "unexpected header: " ~ to!string(found_header));
        }
        int[1] int1; fin.rawRead(int1);
        int label_length = int1[0];
        if (label_length > 0) {
            char[] found_label = new char[label_length];
            fin.rawRead(found_label);
            myLabel = to!string(found_label);
        }
        double[1] dbl1; fin.rawRead(dbl1); sim_time = dbl1[0];
        fin.rawRead(int1); nvariables = int1[0];
        foreach(i; 0 .. nvariables) {
            char[] varname; fin.rawRead(int1); varname.length = int1[0];
            fin.rawRead(varname);
        }
        int[2] int2; fin.rawRead(int2); my_dimensions = int2[0]; nc = int2[1];
        if (my_dimensions != blk.myConfig.dimensions) {
            string msg = text("dimensions found: " ~ to!string(my_dimensions));
            throw new FlowSolverException(msg);
        }
        if (nc != blk.ncells) {
            string msg = text("For block[", blk.id, "] we have a mismatch in solution size.",
                              " Have read nc=", nc, " ncells=", blk.ncells);
            throw new FlowSolverException(msg);
        }
        foreach (i; 0 .. blk.ncells) {
            blk.cells[i].read_values_from_raw_binary(fin, overwrite_geometry_data);
        }
        break;
    default:
        auto byLine = new GzipByLine(filename);
        auto line = byLine.front; byLine.popFront();
        string format_version;
        formattedRead(line, "unstructured_grid_flow %s", &format_version);
        if (format_version != "1.0") {
            string msg = text("UFluidBlock.read_solution(): " ~ "format version found: "
                              ~ format_version);
            throw new FlowSolverException(msg);
        }
        line = byLine.front; byLine.popFront();
        formattedRead(line, "label: %s", &myLabel);
        line = byLine.front; byLine.popFront();
        formattedRead(line, "sim_time: %g", &sim_time);
        line = byLine.front; byLine.popFront();
        formattedRead(line, "variables: %d", &nvariables);
        line = byLine.front; byLine.popFront();
        variableNames = line.strip().split();
        foreach (i; 0 .. variableNames.length) { variableNames[i] = variableNames[i].strip("\""); }
        // If all of the variables line up,
        // it is probably faster to use fixed-order, so default to that.
        bool useFixedOrder = true;
        foreach (i, varName; blk.myConfig.flow_variable_list) {
            if (!std.algorithm.canFind(variableNames, varName)) {
                throw new Exception("Could not find variable: " ~ varName);
            }
            if (!std.algorithm.equal(variableNames[i], varName)) {
                useFixedOrder = false;
            }
        }
        line = byLine.front; byLine.popFront();
        formattedRead(line, "dimensions: %d", &my_dimensions);
        line = byLine.front; byLine.popFront();
        formattedRead(line, "ncells: %d", &nc);
        if (nc != blk.ncells) {
            string msg = text("For block[", blk.id, "] we have a mismatch in solution size.",
                              " Have read nc=", nc, " ncells=", blk.ncells);
            throw new FlowSolverException(msg);
        }
        foreach (i; 0 .. blk.ncells) {
            line = byLine.front; byLine.popFront();
            blk.cells[i].scan_values_from_string(line, variableNames, useFixedOrder,
                                                 blk.myConfig.gmodel, overwrite_geometry_data);
        }
    } // end switch flow_format
    return sim_time;
} // end read_solution() for unstructured-grid block


void write_solution(FluidBlock blk, string filename, double sim_time)
{
    auto sblk = cast(SFluidBlock) blk;
    if (sblk) { return write_solution(sblk, filename, sim_time); }
    auto ublk = cast(UFluidBlock) blk;
    if (ublk) { return write_solution(ublk, filename, sim_time); }
    throw new Error("Oops, unknown type of FluidBlock");
}

void write_solution(SFluidBlock blk, string filename, double sim_time)
// Write the flow solution (i.e. the primary variables at the cell centers)
// for a single block.
// The text format is almost Tecplot POINT format; the raw binary format has the same layout.
// Keep in sync with write_initial_flow_file() in flowstate.d and read_solution above.
{
    if (blk.myConfig.verbosity_level > 1) { writeln("write_solution(): Start block ", blk.id); }
    switch (blk.myConfig.flow_format) {
    case "gziptext": goto default;
    case "rawbinary":
        File outfile = File(filename, "wb");
        int[1] int1; int[4] int4; double[1] dbl1; // buffer arrays
        string header = "structured_grid_flow 1.0";
        outfile.rawWrite(to!(char[])(header));
        int1[0] = to!int(blk.label.length); outfile.rawWrite(int1);
        if (blk.label.length > 0) { outfile.rawWrite(to!(char[])(blk.label)); }
        dbl1[0] = sim_time; outfile.rawWrite(dbl1);
        int1[0] = to!int(blk.myConfig.flow_variable_list.length); outfile.rawWrite(int1);
        foreach(varname; blk.myConfig.flow_variable_list) {
            int1[0] = to!int(varname.length); outfile.rawWrite(int1);
            outfile.rawWrite(to!(char[])(varname));
        }
        int4[0] = blk.myConfig.dimensions;
        int4[1] = to!int(blk.nic); int4[2] = to!int(blk.njc); int4[3] = to!int(blk.nkc);
        outfile.rawWrite(int4);
        foreach (k; 0 .. blk.nkc) {
            foreach (j; 0 .. blk.njc) {
                foreach (i; 0 .. blk.nic) {
                    blk.get_cell(i,j,k).write_values_to_raw_binary(outfile);
                }
            }
        }
        outfile.close();
        break;
    default:
        auto outfile = new GzipOut(filename);
        auto writer = appender!string();
        formattedWrite(writer, "structured_grid_flow 1.0\n");
        formattedWrite(writer, "label: %s\n", blk.label);
        formattedWrite(writer, "sim_time: %.18e\n", sim_time);
        formattedWrite(writer, "variables: %d\n", blk.myConfig.flow_variable_list.length);
        foreach(varname; blk.myConfig.flow_variable_list) { formattedWrite(writer, " \"%s\"", varname); }
        formattedWrite(writer, "\n");
        formattedWrite(writer, "dimensions: %d\n", blk.myConfig.dimensions);
        formattedWrite(writer, "nicell: %d\n", blk.nic);
        formattedWrite(writer, "njcell: %d\n", blk.njc);
        formattedWrite(writer, "nkcell: %d\n", blk.nkc);
        outfile.compress(writer.data);
        foreach (k; 0 .. blk.nkc) {
            foreach (j; 0 .. blk.njc) {
                foreach (i; 0 .. blk.nic) {
                    outfile.compress(" " ~ blk.get_cell(i,j,k).write_values_to_string() ~ "\n");
                }
            }
        }
        outfile.finish();
    } // end switch flow_format
    if (blk.myConfig.verbosity_level > 1) { writeln("write_solution(): Done block ", blk.id); }
} // end write_solution() for structured-grid block


void write_solution(UFluidBlock blk, string filename, double sim_time)
// Write the flow solution (i.e. the primary variables at the cell centers)
// for a single block.
// Keep this function in sync with
// write_initial_usg_flow_file_from_lua() from luaflowstate.d and
// write_initial_flow_file() from flowstate.d.
{
    if (blk.myConfig.verbosity_level > 1) { writeln("write_solution(): Start block ", blk.id); }
    switch (blk.myConfig.flow_format) {
    case "gziptext": goto default;
    case "rawbinary":
        File outfile = File(filename, "wb");
        int[1] int1; int[2] int2; double[1] dbl1; // buffer arrays
        string header = "unstructured_grid_flow 1.0";
        outfile.rawWrite(to!(char[])(header));
        int1[0] = to!int(blk.label.length); outfile.rawWrite(int1);
        if (blk.label.length > 0) { outfile.rawWrite(to!(char[])(blk.label)); }
        dbl1[0] = sim_time; outfile.rawWrite(dbl1);
        int1[0] = to!int(blk.myConfig.flow_variable_list.length); outfile.rawWrite(int1);
        foreach(varname; blk.myConfig.flow_variable_list) {
            int1[0] = to!int(varname.length); outfile.rawWrite(int1);
            outfile.rawWrite(to!(char[])(varname));
        }
        int2[0] = blk.myConfig.dimensions; int2[1] = to!int(blk.ncells); outfile.rawWrite(int2);
        foreach(cell; blk.cells) { cell.write_values_to_raw_binary(outfile); }
        outfile.close();
        break;
    default:
        auto outfile = new GzipOut(filename);
        auto writer = appender!string();
        formattedWrite(writer, "unstructured_grid_flow 1.0\n");
        formattedWrite(writer, "label: %s\n", blk.label);
        formattedWrite(writer, "sim_time: %.18e\n", sim_time);
        formattedWrite(writer, "variables: %d\n", blk.myConfig.flow_variable_list.length);
        // Variable list for cell on one line.
        foreach(varname; blk.myConfig.flow_variable_list) {
            formattedWrite(writer, " \"%s\"", varname);
        }
        formattedWrite(writer, "\n");
        // Numbers of cells
        formattedWrite(writer, "dimensions: %d\n", blk.myConfig.dimensions);
        formattedWrite(writer, "ncells: %d\n", blk.ncells);
        outfile.compress(writer.data);
        // The actual cell data.
        foreach(cell; blk.cells) {
            outfile.compress(" " ~ cell.write_values_to_string() ~ "\n");
        }
        outfile.finish();
    } // end switch flow_format
} // end write_solution() for unstructured-grid block


void write_residuals(FluidBlock blk, string filename)
{
    auto sblk = cast(SFluidBlock) blk;
    if (sblk) { return write_residuals(sblk, filename); }
    auto ublk = cast(UFluidBlock) blk;
    if (ublk) { return write_residuals(ublk, filename); }
    throw new Error("Oops, unknown type of FluidBlock");
}

void write_residuals(SFluidBlock blk, string filename)
{
    auto outfile = new GzipOut(filename);
    auto writer = appender!string();
    formattedWrite(writer, "structured_grid_residuals 1.0\n");
    formattedWrite(writer, "label: %s\n", blk.label);
    // RJG: Fix me, we'll need to make this variable based on conservation
    //      equations being solved.
    //      For now, assume the simplest: single-species ideal gas in two dimensions
    formattedWrite(writer, "variables: %d\n", 4);
    string[] varnames = ["density", "x-momentum", "y-momentum", "energy"];
    foreach (var; varnames) {
        formattedWrite(writer, " \"%s\"", var);
    }
    formattedWrite(writer, "\n");
    formattedWrite(writer, "dimensions: %d\n", blk.myConfig.dimensions);
    formattedWrite(writer, "nicell: %d\n", blk.nic);
    formattedWrite(writer, "njcell: %d\n", blk.njc);
    formattedWrite(writer, "nkcell: %d\n", blk.nkc);
    outfile.compress(writer.data);
    foreach (cell; blk.cells) {
        outfile.compress(" " ~ cell.write_residuals_to_string() ~ "\n");
    }
    outfile.finish();
}

void write_residuals(UFluidBlock blk, string filename)
{
    auto outfile = new GzipOut(filename);
    auto writer = appender!string();
    formattedWrite(writer, "unstructured_grid_residuals 1.0\n");
    formattedWrite(writer, "label: %s\n", blk.label);
    // RJG: Fix me, we'll need to make this variable based on conservation
    //      equations being solved.
    //      For now, assume the simplest: single-species ideal gas in two dimensions
    formattedWrite(writer, "variables: %d\n", 4);
    string[] varnames = ["density", "x-momentum", "y-momentum", "energy"];
    foreach (var; varnames) {
        formattedWrite(writer, " \"%s\"", var);
    }
    formattedWrite(writer, "\n");
    formattedWrite(writer, "dimensions: %d\n", blk.myConfig.dimensions);
    formattedWrite(writer, "ncells: %d\n", blk.ncells);
    outfile.compress(writer.data);
    foreach (cell; blk.cells) {
        outfile.compress(" " ~ cell.write_residuals_to_string() ~ "\n");
    }
    outfile.finish();
}

void write_DFT(FluidBlock blk, string filename)
{
    auto outfile = new GzipOut(filename);
    auto writer = appender!string();
    formattedWrite(writer, "# Discrete Fourier Transform, Cell-by-Cell\n");
    outfile.compress(writer.data);
    foreach (cell; blk.cells) {
        outfile.compress(" " ~ cell.write_DFT_to_string() ~ "\n");
    }
    outfile.finish();
}

version(shape_sensitivity) {
    void write_adjoint_variables(FluidBlock blk, string filename)
    {
        auto sblk = cast(SFluidBlock) blk;
        if (sblk) { return write_adjoint_variables(sblk, filename); }
        auto ublk = cast(UFluidBlock) blk;
        if (ublk) { return write_adjoint_variables(ublk, filename); }
        throw new Error("Oops, unknown type of FluidBlock");
    }

    void write_adjoint_variables(SFluidBlock blk, string filename)
    {
        throw new FlowSolverException("write_adjoint_variables not implemented for Structured blocks.");
    }

    void write_adjoint_variables(UFluidBlock blk, string filename)
    {
        auto outfile = new GzipOut(filename);
        auto writer = appender!string();
        formattedWrite(writer, "unstructured_grid_adjoint_variables 1.0\n");
        formattedWrite(writer, "label: %s\n", blk.label);
        // RJG: Fix me, we'll need to make this variable based on conservation
        //      equations being solved.
        //      For now, assume the simplest: single-species ideal gas in two dimensions
        formattedWrite(writer, "variables: %d\n", 4);
        string[] varnames = ["density", "x-momentum", "y-momentum", "total-energy"];
        foreach (var; varnames) {
            formattedWrite(writer, " \"%s\"", var);
        }
        formattedWrite(writer, "\n");
        formattedWrite(writer, "dimensions: %d\n", blk.myConfig.dimensions);
        formattedWrite(writer, "ncells: %d\n", blk.ncells);
        outfile.compress(writer.data);
        // RJG: FIX ME
        //      As per hard-coded assumption made above.
        int np = 4;
        foreach (i; 0 .. blk.ncells) {
            auto s = format!"%.18e %.18e %.18e %.18e\n"(blk.psi[np*i+blk.MASS].re,
                                                        blk.psi[np*i+blk.X_MOM].re,
                                                        blk.psi[np*i+blk.Y_MOM].re,
                                                        blk.psi[np*i+blk.TOT_ENERGY].re);
            outfile.compress(s);
        }
        outfile.finish();
    }
}


void write_initial_flow_file(string fileName, ref StructuredGrid grid,
                             in FlowState fs, double t0, double omegaz, GasModel gmodel)
// Keep in sync with SFluidBlock.write_solution.
// Note that we assume the FlowState fs to be in a nonrotating frame.
{
    // Numbers of cells derived from numbers of vertices in grid.
    auto nicell = grid.niv - 1;
    auto njcell = grid.njv - 1;
    auto nkcell = grid.nkv - 1;
    if (GlobalConfig.dimensions == 2) nkcell = 1;
    //
    auto myfs = new FlowState(fs);
    //
    // Write the data for the whole structured block.
    switch (GlobalConfig.flow_format) {
    case "gziptext": goto default;
    case "rawbinary":
        File outfile = File(fileName, "wb");
        int[1] int1; int[4] int4; double[1] dbl1; // buffer arrays
        string header = "structured_grid_flow 1.0";
        outfile.rawWrite(to!(char[])(header));
        int1[0] = to!int(grid.label.length); outfile.rawWrite(int1);
        if (grid.label.length > 0) { outfile.rawWrite(to!(char[])(grid.label)); }
        dbl1[0] = t0; outfile.rawWrite(dbl1); // sim_time
        int1[0] = to!int(GlobalConfig.flow_variable_list.length); outfile.rawWrite(int1);
        foreach(varname; GlobalConfig.flow_variable_list) {
            int1[0] = to!int(varname.length); outfile.rawWrite(int1);
            outfile.rawWrite(to!(char[])(varname));
        }
        int4[0] = to!int(GlobalConfig.dimensions);
        int4[1] = to!int(nicell); int4[2] = to!int(njcell); int4[3] = to!int(nkcell);
        outfile.rawWrite(int4);
        foreach (k; 0 .. nkcell) {
            foreach (j; 0 .. njcell) {
                foreach (i; 0 .. nicell) {
                    Vector3 p000 = *grid[i,j,k];
                    Vector3 p100 = *grid[i+1,j,k];
                    Vector3 p110 = *grid[i+1,j+1,k];
                    Vector3 p010 = *grid[i,j+1,k];
                    Vector3 pos;
                    number volume, iLen, jLen, kLen;
                    if (GlobalConfig.dimensions == 2) {
                        number xyplane_area;
                        xyplane_quad_cell_properties(p000, p100, p110, p010, pos, xyplane_area, iLen, jLen, kLen);
                        volume = xyplane_area * ((GlobalConfig.axisymmetric) ? pos.y : to!number(1.0) );
                    } else if (GlobalConfig.dimensions == 3) {
                        Vector3 p001 = *grid[i,j,k+1];
                        Vector3 p101 = *grid[i+1,j,k+1];
                        Vector3 p111 = *grid[i+1,j+1,k+1];
                        Vector3 p011 = *grid[i,j+1,k+1];
                        hex_cell_properties(p000, p100, p110, p010, p001, p101, p111, p011, pos, volume, iLen, jLen, kLen);
                        if (omegaz != 0.0) { myfs.vel.set(fs.vel); into_rotating_frame(myfs.vel, pos, omegaz); }
                    } else {
                        throw new Exception("GlobalConfig.dimensions not 2 or 3.");
                    }
                    cell_data_to_raw_binary(outfile, pos, volume, myfs,
                                            to!number(0.0), to!number(0.0), to!number(0.0),
                                            GlobalConfig.with_local_time_stepping,
                                            -1.0, -1.0, -1.0,
                                            GlobalConfig.include_quality,
                                            GlobalConfig.MHD,
                                            GlobalConfig.divergence_cleaning,
                                            GlobalConfig.radiation,
                                            GlobalConfig.turb_model.nturb);
                }
            }
        }
        outfile.close();
        break;
    default:
        auto outfile = new GzipOut(fileName);
        auto writer = appender!string();
        formattedWrite(writer, "structured_grid_flow 1.0\n");
        formattedWrite(writer, "label: %s\n", grid.label);
        formattedWrite(writer, "sim_time: %.18e\n", t0);
        formattedWrite(writer, "variables: %d\n", GlobalConfig.flow_variable_list.length);
        // Variable list for cell on one line.
        foreach(varname; GlobalConfig.flow_variable_list) {
            formattedWrite(writer, " \"%s\"", varname);
        }
        formattedWrite(writer, "\n");
        // Numbers of cells
        formattedWrite(writer, "dimensions: %d\n", GlobalConfig.dimensions);
        formattedWrite(writer, "nicell: %d\n", nicell);
        formattedWrite(writer, "njcell: %d\n", njcell);
        formattedWrite(writer, "nkcell: %d\n", nkcell);
        outfile.compress(writer.data);
        // The actual cell data.
        foreach (k; 0 .. nkcell) {
            foreach (j; 0 .. njcell) {
                foreach (i; 0 .. nicell) {
                    Vector3 p000 = *grid[i,j,k];
                    Vector3 p100 = *grid[i+1,j,k];
                    Vector3 p110 = *grid[i+1,j+1,k];
                    Vector3 p010 = *grid[i,j+1,k];
                    Vector3 pos;
                    number volume, iLen, jLen, kLen;
                    if (GlobalConfig.dimensions == 2) {
                        number xyplane_area;
                        xyplane_quad_cell_properties(p000, p100, p110, p010, pos, xyplane_area, iLen, jLen, kLen);
                        volume = xyplane_area * ((GlobalConfig.axisymmetric) ? pos.y : to!number(1.0) );
                    } else if (GlobalConfig.dimensions == 3) {
                        Vector3 p001 = *grid[i,j,k+1];
                        Vector3 p101 = *grid[i+1,j,k+1];
                        Vector3 p111 = *grid[i+1,j+1,k+1];
                        Vector3 p011 = *grid[i,j+1,k+1];
                        hex_cell_properties(p000, p100, p110, p010, p001, p101, p111, p011, pos, volume, iLen, jLen, kLen);
                        if (omegaz != 0.0) { myfs.vel.set(fs.vel); into_rotating_frame(myfs.vel, pos, omegaz); }
                    } else {
                        throw new Exception("GlobalConfig.dimensions not 2 or 3.");
                    }
                    outfile.compress(" " ~ cell_data_as_string(pos, volume, myfs,
                                                               to!number(0.0), to!number(0.0), to!number(0.0),
                                                               GlobalConfig.with_local_time_stepping, -1.0, -1.0, -1.0,
                                                               GlobalConfig.include_quality,
                                                               GlobalConfig.MHD,
                                                               GlobalConfig.divergence_cleaning,
                                                               GlobalConfig.radiation,
                                                               GlobalConfig.turb_model.nturb) ~ "\n");
                }
            }
        }
        outfile.finish();
    } // end switch flow_format
    return;
} // end write_initial_flow_file() StructuredGrid version

void write_initial_flow_file(string fileName, ref UnstructuredGrid grid,
                             in FlowState fs, double t0, double omegaz, GasModel gmodel)
// Keep in sync with UFluidBlock.write_solution.
// Note that we assume the FlowState fs to be in a nonrotating frame.
{
    // Numbers of cells derived from numbers of vertices in grid.
    auto ncells = grid.ncells;
    //
    // Write the data for the whole unstructured block.
    switch (GlobalConfig.flow_format) {
    case "gziptext": goto default;
    case "rawbinary":
        File outfile = File(fileName, "wb");
        int[1] int1; int[2] int2; double[1] dbl1; // buffer arrays
        string header = "unstructured_grid_flow 1.0";
        outfile.rawWrite(to!(char[])(header));
        int1[0] = to!int(grid.label.length); outfile.rawWrite(int1);
        if (grid.label.length > 0) { outfile.rawWrite(to!(char[])(grid.label)); }
        dbl1[0] = t0; outfile.rawWrite(dbl1); // sim_time
        int1[0] = to!int(GlobalConfig.flow_variable_list.length); outfile.rawWrite(int1);
        foreach(varname; GlobalConfig.flow_variable_list) {
            int1[0] = to!int(varname.length); outfile.rawWrite(int1);
            outfile.rawWrite(to!(char[])(varname));
        }
        int2[0] = to!int(GlobalConfig.dimensions);
        int2[1] = to!int(ncells);
        outfile.rawWrite(int2);
        foreach (i; 0 .. ncells) {
            Vector3 pos = Vector3(0.0, 0.0, 0.0);
            foreach (id; grid.cells[i].vtx_id_list) { pos += grid.vertices[id]; }
            pos /= to!number(grid.cells[i].vtx_id_list.length);
            number volume = 0.0;
            if (omegaz != 0.0) {
                throw new Error("Oops, we have not yet implemented rotating-frame code here.");
            }
            cell_data_to_raw_binary(outfile, pos, volume, fs,
                                    to!number(0.0), to!number(0.0), to!number(0.0),
                                    GlobalConfig.with_local_time_stepping, -1.0, -1.0, -1.0,
                                    GlobalConfig.include_quality,
                                    GlobalConfig.MHD,
                                    GlobalConfig.divergence_cleaning,
                                    GlobalConfig.radiation,
                                    GlobalConfig.turb_model.nturb);
        }
        outfile.close();
        break;
    default:
        auto outfile = new GzipOut(fileName);
        auto writer = appender!string();
        formattedWrite(writer, "unstructured_grid_flow 1.0\n");
        formattedWrite(writer, "label: %s\n", grid.label);
        formattedWrite(writer, "sim_time: %.18e\n", t0);
        formattedWrite(writer, "variables: %d\n", GlobalConfig.flow_variable_list.length);
        // Variable list for cell on one line.
        foreach(varname; GlobalConfig.flow_variable_list) {
            formattedWrite(writer, " \"%s\"", varname);
        }
        formattedWrite(writer, "\n");
        // Numbers of cells
        formattedWrite(writer, "dimensions: %d\n", GlobalConfig.dimensions);
        formattedWrite(writer, "ncells: %d\n", ncells);
        outfile.compress(writer.data);
        // The actual cell data.
        foreach (i; 0 .. ncells) {
            Vector3 pos = Vector3(0.0, 0.0, 0.0);
            foreach (id; grid.cells[i].vtx_id_list) { pos += grid.vertices[id]; }
            pos /= to!number(grid.cells[i].vtx_id_list.length);
            number volume = 0.0;
            if (omegaz != 0.0) {
                throw new Error("Oops, we have not yet implemented rotating-frame code here.");
            }
            outfile.compress(" " ~ cell_data_as_string(pos, volume, fs,
                                                       to!number(0.0), to!number(0.0), to!number(0.0),
                                                       GlobalConfig.with_local_time_stepping, -1.0, -1.0, -1.0,
                                                       GlobalConfig.include_quality,
                                                       GlobalConfig.MHD,
                                                       GlobalConfig.divergence_cleaning,
                                                       GlobalConfig.radiation,
                                                       GlobalConfig.turb_model.nturb) ~ "\n");
        }
        outfile.finish();
    } // end switch flow_format
    return;
} // end write_initial_flow_file() UnstructuredGrid version


extern(C) int write_initial_sg_flow_file_from_lua(lua_State* L)
{
    auto fname = to!string(luaL_checkstring(L, 1));
    auto grid = checkStructuredGrid(L, 2);
    double t0 = luaL_checknumber(L, 4);
    double omegaz = luaL_checknumber(L, 5);
    if (GlobalConfig.flow_variable_list.length == 0) {
        foreach(varname; build_flow_variable_list()) { GlobalConfig.flow_variable_list ~= varname; }
    }
    FlowState fs;
    // Test if we have a simple flow state or something more exotic
    if ( isObjType(L, 3, "_FlowState") ) {
        fs = checkFlowState(L, 3);
        write_initial_flow_file(fname, grid, fs, t0, omegaz, GlobalConfig.gmodel_master);
        lua_settop(L, 0); // clear stack
        return 0;
    }
    // Else, we might have a callable lua function
    if (lua_isfunction(L, 3)) {
        // Assume we can use the function and also assume that it provides the
        // flowstate with the velocity already in the rotating frame, if omegaz
        // for the block is nonzero.
        //
        // A lot of code borrowed from flowstate.d
        // Keep in sync with write_initial_flow_file() function in that file.
        //
        // Numbers of cells derived from numbers of vertices in grid.
        auto nicell = grid.niv - 1;
        auto njcell = grid.njv - 1;
        auto nkcell = grid.nkv - 1;
        if (GlobalConfig.dimensions == 2) nkcell = 1;
        //
        // Write the data for the whole structured block.
        auto gmodel = GlobalConfig.gmodel_master;
        switch (GlobalConfig.flow_format) {
        case "gziptext": goto default;
        case "rawbinary":
            File outfile = File(fname, "wb");
            int[1] int1; int[4] int4; double[1] dbl1; // buffer arrays
            string header = "structured_grid_flow 1.0";
            outfile.rawWrite(to!(char[])(header));
            int1[0] = to!int(grid.label.length); outfile.rawWrite(int1);
            if (grid.label.length > 0) { outfile.rawWrite(to!(char[])(grid.label)); }
            dbl1[0] = t0; outfile.rawWrite(dbl1); // sim_time
            int1[0] = to!int(GlobalConfig.flow_variable_list.length); outfile.rawWrite(int1);
            foreach(varname; GlobalConfig.flow_variable_list) {
                int1[0] = to!int(varname.length); outfile.rawWrite(int1);
                outfile.rawWrite(to!(char[])(varname));
            }
            int4[0] = to!int(GlobalConfig.dimensions);
            int4[1] = to!int(nicell); int4[2] = to!int(njcell); int4[3] = to!int(nkcell);
            outfile.rawWrite(int4);
            foreach (k; 0 .. nkcell) {
                foreach (j; 0 .. njcell) {
                    foreach (i; 0 .. nicell) {
                        Vector3 p000 = *grid[i,j,k];
                        Vector3 p100 = *grid[i+1,j,k];
                        Vector3 p110 = *grid[i+1,j+1,k];
                        Vector3 p010 = *grid[i,j+1,k];
                        Vector3 pos;
                        number volume, iLen, jLen, kLen;
                        if (GlobalConfig.dimensions == 2) {
                            number xyplane_area;
                            xyplane_quad_cell_properties(p000, p100, p110, p010, pos, xyplane_area, iLen, jLen, kLen);
                            volume = xyplane_area * ((GlobalConfig.axisymmetric) ? pos.y : to!number(1.0) );
                        } else if (GlobalConfig.dimensions == 3) {
                            Vector3 p001 = *grid[i,j,k+1];
                            Vector3 p101 = *grid[i+1,j,k+1];
                            Vector3 p111 = *grid[i+1,j+1,k+1];
                            Vector3 p011 = *grid[i,j+1,k+1];
                            hex_cell_properties(p000, p100, p110, p010, p001, p101, p111, p011, pos, volume, iLen, jLen, kLen);
                        } else {
                            throw new Exception("GlobalConfig.dimensions not 2 or 3.");
                        }
                        // Now grab flow state via Lua function call.
                        // If the block is in a rotating frame with omegaz != 0.0,
                        // we presume that the Lua function will provide the velocity
                        // components relative to the rotating frame.
                        lua_pushvalue(L, 3);
                        lua_pushnumber(L, pos.x);
                        lua_pushnumber(L, pos.y);
                        lua_pushnumber(L, pos.z);
                        if (lua_pcall(L, 3, 1, 0) != 0) {
                            string errMsg = "Error in Lua function call for setting FlowState\n";
                            errMsg ~= "as a function of position (x, y, z).\n";
                            luaL_error(L, errMsg.toStringz);
                        }
                        if (lua_istable(L, -1)) {
                            fs = makeFlowStateFromTable(L, lua_gettop(L));
                        } else {
                            fs = checkFlowState(L, -1);
                        }
                        if (!fs) {
                            string errMsg = "Error in from Lua function call for setting FlowState\n";
                            errMsg ~= "as a function of position (x, y, z).\n";
                            errMsg ~= "The returned object is not a proper _FlowState object or table.";
                            luaL_error(L, errMsg.toStringz);
                        }
                        cell_data_to_raw_binary(outfile, pos, volume, fs,
                                                to!number(0.0), to!number(0.0), to!number(0.0),
                                                GlobalConfig.with_local_time_stepping, -1.0, -1.0, -1.0,
                                                GlobalConfig.include_quality,
                                                GlobalConfig.MHD,
                                                GlobalConfig.divergence_cleaning,
                                                GlobalConfig.radiation,
                                                GlobalConfig.turb_model.nturb);
                    }
                }
            }
            outfile.close();
            break;
        default:
            auto outfile = new GzipOut(fname);
            auto writer = appender!string();
            formattedWrite(writer, "structured_grid_flow 1.0\n");
            formattedWrite(writer, "label: %s\n", grid.label);
            formattedWrite(writer, "sim_time: %.18e\n", t0);
            formattedWrite(writer, "variables: %d\n", GlobalConfig.flow_variable_list.length);
            // Variable list for cell on one line.
            foreach(varname; GlobalConfig.flow_variable_list) {
                formattedWrite(writer, " \"%s\"", varname);
            }
            formattedWrite(writer, "\n");
            // Numbers of cells
            formattedWrite(writer, "dimensions: %d\n", GlobalConfig.dimensions);
            formattedWrite(writer, "nicell: %d\n", nicell);
            formattedWrite(writer, "njcell: %d\n", njcell);
            formattedWrite(writer, "nkcell: %d\n", nkcell);
            outfile.compress(writer.data);
            // The actual cell data.
            foreach (k; 0 .. nkcell) {
                foreach (j; 0 .. njcell) {
                    foreach (i; 0 .. nicell) {
                        Vector3 p000 = *grid[i,j,k];
                        Vector3 p100 = *grid[i+1,j,k];
                        Vector3 p110 = *grid[i+1,j+1,k];
                        Vector3 p010 = *grid[i,j+1,k];
                        Vector3 pos;
                        number volume, iLen, jLen, kLen;
                        if (GlobalConfig.dimensions == 2) {
                            number xyplane_area;
                            xyplane_quad_cell_properties(p000, p100, p110, p010, pos, xyplane_area, iLen, jLen, kLen);
                            volume = xyplane_area * ((GlobalConfig.axisymmetric) ? pos.y : to!number(1.0) );
                        } else if (GlobalConfig.dimensions == 3) {
                            Vector3 p001 = *grid[i,j,k+1];
                            Vector3 p101 = *grid[i+1,j,k+1];
                            Vector3 p111 = *grid[i+1,j+1,k+1];
                            Vector3 p011 = *grid[i,j+1,k+1];
                            hex_cell_properties(p000, p100, p110, p010, p001, p101, p111, p011, pos, volume, iLen, jLen, kLen);
                        } else {
                            throw new Exception("GlobalConfig.dimensions not 2 or 3.");
                        }
                        // Now grab flow state via Lua function call.
                        // If the block is in a rotating frame with omegaz != 0.0,
                        // we presume that the Lua function will provide the velocity
                        // components relative to the rotating frame.
                        lua_pushvalue(L, 3);
                        lua_pushnumber(L, pos.x);
                        lua_pushnumber(L, pos.y);
                        lua_pushnumber(L, pos.z);
                        if (lua_pcall(L, 3, 1, 0) != 0) {
                            string errMsg = "Error in Lua function call for setting FlowState\n";
                            errMsg ~= "as a function of position (x, y, z).\n";
                            luaL_error(L, errMsg.toStringz);
                        }
                        if (lua_istable(L, -1)) {
                            fs = makeFlowStateFromTable(L, lua_gettop(L));
                        } else {
                            fs = checkFlowState(L, -1);
                        }
                        if (!fs) {
                            string errMsg = "Error in from Lua function call for setting FlowState\n";
                            errMsg ~= "as a function of position (x, y, z).\n";
                            errMsg ~= "The returned object is not a proper _FlowState object or suitable table.";
                            luaL_error(L, errMsg.toStringz);
                        }
                        outfile.compress(" " ~ cell_data_as_string(pos, volume, fs,
                                                                   to!number(0.0), to!number(0.0), to!number(0.0),
                                                                   GlobalConfig.with_local_time_stepping, -1.0, -1.0, -1.0,
                                                                   GlobalConfig.include_quality,
                                                                   GlobalConfig.MHD,
                                                                   GlobalConfig.divergence_cleaning,
                                                                   GlobalConfig.radiation,
                                                                   GlobalConfig.turb_model.nturb) ~ "\n");
                    } // end foreach i
                } // end foreach j
            } // end foreach k
            outfile.finish();
        } // end switch flow_format
        lua_settop(L, 0); // clear stack
        return 0;
    } // end if lua_isfunction
    lua_settop(L, 0); // clear stack
    return -1;
} // end write_initial_sg_flow_file_from_lua()

extern(C) int write_initial_usg_flow_file_from_lua(lua_State* L)
{
    auto fname = to!string(luaL_checkstring(L, 1));
    auto grid = checkUnstructuredGrid(L, 2);
    double t0 = luaL_checknumber(L, 4);
    double omegaz = luaL_checknumber(L, 5);
    if (GlobalConfig.flow_variable_list.length == 0) {
        foreach(varname; build_flow_variable_list()) { GlobalConfig.flow_variable_list ~= varname; }
    }
    FlowState fs;
    // Test if we have a simple flow state or something more exotic
    if ( isObjType(L, 3, "_FlowState") ) {
        fs = checkFlowState(L, 3);
        write_initial_flow_file(fname, grid, fs, t0, omegaz, GlobalConfig.gmodel_master);
        lua_settop(L, 0); // clear stack
        return 0;
    }
    // Else, we might have a callable lua function
    if ( lua_isfunction(L, 3) ) {
        // Assume we can use the function then.
        // A lot of code borrowed from flowstate.d
        // Keep in sync with write_initial_flow_file() function in that file.
        //
        // Numbers of cells derived from numbers of vertices in grid.
        auto ncells = grid.ncells;
        //
        // Write the data for the whole unstructured block.
        auto gmodel = GlobalConfig.gmodel_master;
        switch (GlobalConfig.flow_format) {
        case "gziptext": goto default;
        case "rawbinary":
            File outfile = File(fname, "wb");
            int[1] int1; int[2] int2; double[1] dbl1; // buffer arrays
            string header = "unstructured_grid_flow 1.0";
            outfile.rawWrite(to!(char[])(header));
            int1[0] = to!int(grid.label.length); outfile.rawWrite(int1);
            if (grid.label.length > 0) { outfile.rawWrite(to!(char[])(grid.label)); }
            dbl1[0] = t0; outfile.rawWrite(dbl1); // sim_time
            int1[0] = to!int(GlobalConfig.flow_variable_list.length); outfile.rawWrite(int1);
            foreach(varname; GlobalConfig.flow_variable_list) {
                int1[0] = to!int(varname.length); outfile.rawWrite(int1);
                outfile.rawWrite(to!(char[])(varname));
            }
            int2[0] = to!int(GlobalConfig.dimensions);
            int2[1] = to!int(ncells);
            outfile.rawWrite(int2);
            foreach (i; 0 .. ncells) {
                Vector3 pos = Vector3(0.0, 0.0, 0.0);
                foreach (id; grid.cells[i].vtx_id_list) { pos += grid.vertices[id]; }
                pos /= grid.cells[i].vtx_id_list.length;
                number volume = 0.0;
                // Now grab flow state via Lua function call
                lua_pushvalue(L, 3);
                lua_pushnumber(L, pos.x);
                lua_pushnumber(L, pos.y);
                lua_pushnumber(L, pos.z);
                if (lua_pcall(L, 3, 1, 0) != 0) {
                    string errMsg = "Error in Lua function call for setting FlowState\n";
                    errMsg ~= "as a function of position (x, y, z).\n";
                    errMsg ~= format("LUA ERROR: %s\n", lua_tostring(L, -1));
                    luaL_error(L, errMsg.toStringz);
                }
                if (lua_istable(L, -1)) {
                    fs = makeFlowStateFromTable(L, lua_gettop(L));
                } else {
                    fs = checkFlowState(L, -1);
                }
                if (!fs) {
                    string errMsg = "Error in from Lua function call for setting FlowState\n";
                    errMsg ~= "as a function of position (x, y, z).\n";
                    errMsg ~= "The returned object is not a proper _FlowState object or a suitable table.";
                    luaL_error(L, errMsg.toStringz);
                }
                cell_data_to_raw_binary(outfile, pos, volume, fs,
                                        to!number(0.0), to!number(0.0), to!number(0.0),
                                        GlobalConfig.with_local_time_stepping, -1.0, -1.0, -1.0,
                                        GlobalConfig.include_quality,
                                        GlobalConfig.MHD,
                                        GlobalConfig.divergence_cleaning,
                                        GlobalConfig.radiation,
                                        GlobalConfig.turb_model.nturb);
            }
            outfile.close();
        break;
        default:
            auto outfile = new GzipOut(fname);
            auto writer = appender!string();
            formattedWrite(writer, "unstructured_grid_flow 1.0\n");
            formattedWrite(writer, "label: %s\n", grid.label);
            formattedWrite(writer, "sim_time: %.18e\n", t0);
            formattedWrite(writer, "variables: %d\n", GlobalConfig.flow_variable_list.length);
            // Variable list for cell on one line.
            foreach(varname; GlobalConfig.flow_variable_list) {
                formattedWrite(writer, " \"%s\"", varname);
            }
            formattedWrite(writer, "\n");
            // Numbers of cells
            formattedWrite(writer, "dimensions: %d\n", GlobalConfig.dimensions);
            formattedWrite(writer, "ncells: %d\n", ncells);
            outfile.compress(writer.data);
            // The actual cell data.
            foreach (i; 0 .. ncells) {
                Vector3 pos = Vector3(0.0, 0.0, 0.0);
                foreach (id; grid.cells[i].vtx_id_list) { pos += grid.vertices[id]; }
                pos /= grid.cells[i].vtx_id_list.length;
                number volume = 0.0;
                // Now grab flow state via Lua function call
                lua_pushvalue(L, 3);
                lua_pushnumber(L, pos.x);
                lua_pushnumber(L, pos.y);
                lua_pushnumber(L, pos.z);
                if (lua_pcall(L, 3, 1, 0) != 0) {
                    string errMsg = "Error in Lua function call for setting FlowState\n";
                    errMsg ~= "as a function of position (x, y, z).\n";
                    errMsg ~= format("LUA ERROR: %s\n", lua_tostring(L, -1));
                    luaL_error(L, errMsg.toStringz);
                }
                if (lua_istable(L, -1)) {
                    fs = makeFlowStateFromTable(L, lua_gettop(L));
                } else {
                    fs = checkFlowState(L, -1);
                }
                if (!fs) {
                    string errMsg = "Error in from Lua function call for setting FlowState\n";
                    errMsg ~= "as a function of position (x, y, z).\n";
                    errMsg ~= "The returned object is not a proper _FlowState object or suitable table.";
                    luaL_error(L, errMsg.toStringz);
                }
                outfile.compress(" " ~ cell_data_as_string(pos, volume, fs,
                                                           to!number(0.0), to!number(0.0), to!number(0.0),
                                                           GlobalConfig.with_local_time_stepping, -1.0, -1.0, -1.0,
                                                           GlobalConfig.include_quality,
                                                           GlobalConfig.MHD,
                                                           GlobalConfig.divergence_cleaning,
                                                           GlobalConfig.radiation,
                                                           GlobalConfig.turb_model.nturb) ~ "\n");
            } // end foreach i
            outfile.finish();
        }
        lua_settop(L, 0); // clear stack
        return 0;
    } // end if lua_isfunction
    lua_settop(L, 0); // clear stack
    return -1;
} // end write_initial_usg_flow_file_from_lua()


void read_block_data_from_file(BlockFlow blk, string filename, Grid_t gridType, string flow_format)
{
    string myLabel;
    size_t nvariables;
    // Read in the flow data for a single block.
    //
    // Keep in sync with:
    // 1. SFluidBlock.write_solution(),
    // 2. UFluidBlock.write_solution()
    // 3. write_initial_sg_flow_file_from_lua() in luaflowstate.d
    // 4. write_initial_usg_flow_file_from_lua() in luaflowstate.d.
    //
    switch (flow_format) {
    case "gziptext": goto default;
    case "rawbinary":
        File fin = File(filename, "rb");
        string expected_header;
        final switch (gridType) {
        case Grid_t.structured_grid: expected_header = "structured_grid_flow 1.0"; break;
        case Grid_t.unstructured_grid: expected_header = "unstructured_grid_flow 1.0";
        }
        char[] found_header = new char[expected_header.length];
        fin.rawRead(found_header);
        if (found_header != expected_header) {
            throw new FlowSolverException("BlockFlow constructor, read_solution from raw_binary_file: "
                                          ~ "unexpected header: " ~ to!string(found_header));
        }
        int[1] int1; fin.rawRead(int1);
        int label_length = int1[0];
        if (label_length > 0) {
            char[] found_label = new char[label_length];
            fin.rawRead(found_label);
            myLabel = to!string(found_label);
        }
        double[1] dbl1; fin.rawRead(dbl1); blk.sim_time = dbl1[0];
        fin.rawRead(int1); nvariables = int1[0];
        foreach(i; 0 .. nvariables) {
            char[] varname; fin.rawRead(int1); varname.length = int1[0];
            fin.rawRead(varname); blk.variableNames ~= to!string(varname);
        }
        fin.rawRead(int1); int my_dimensions = int1[0];
        final switch (gridType) {
        case Grid_t.structured_grid:
            int[3] int3; fin.rawRead(int3);
            blk.nic = int3[0]; blk.njc = int3[1]; blk.nkc = int3[2];
            blk.ncells = blk.nic*blk.njc*blk.nkc;
            break;
        case Grid_t.unstructured_grid:
            fin.rawRead(int1); blk.ncells = int1[0];
            blk.nic = blk.ncells; blk.njc = 1; blk.nkc = 1;
        }
        // Assume the remaining data is all double type and is in standard cell order.
        blk._data.length = blk.ncells;
        foreach (i; 0 .. blk.ncells) {
            blk._data[i].length = nvariables;
            fin.rawRead(blk._data[i]);
        }
        fin.close();
        break;
    default:
        string[] tokens;
        auto byLine = new GzipByLine(filename);
        auto line = byLine.front; byLine.popFront();
        string format_version;
        final switch (gridType) {
        case Grid_t.structured_grid:
            formattedRead(line, "structured_grid_flow %s", &format_version);
            break;
        case Grid_t.unstructured_grid:
            formattedRead(line, "unstructured_grid_flow %s", &format_version);
        }
        if (format_version != "1.0") {
            string msg = text("BlockFlow.read_solution(): " ~
                              "format version found: " ~ format_version);
            throw new FlowSolverException(msg);
        }
        line = byLine.front; byLine.popFront();
        formattedRead(line, "label: %s", &myLabel);
        line = byLine.front; byLine.popFront();
        formattedRead(line, "sim_time: %g", &blk.sim_time);
        line = byLine.front; byLine.popFront();
        formattedRead(line, "variables: %d", &nvariables);
        line = byLine.front; byLine.popFront();
        blk.variableNames = line.strip().split();
        foreach (ref var; blk.variableNames) { var = var.replaceAll(regex("\""), ""); }
        line = byLine.front; byLine.popFront();
        formattedRead(line, "dimensions: %d", &blk.dimensions);
        final switch (gridType) {
        case Grid_t.structured_grid:
            line = byLine.front; byLine.popFront();
            formattedRead(line, "nicell: %d", &blk.nic);
            line = byLine.front; byLine.popFront();
            formattedRead(line, "njcell: %d", &blk.njc);
            line = byLine.front; byLine.popFront();
            formattedRead(line, "nkcell: %d", &blk.nkc);
            blk.ncells = blk.nic*blk.njc*blk.nkc;
            break;
        case Grid_t.unstructured_grid:
            line = byLine.front; byLine.popFront();
            formattedRead(line, "ncells: %d", &blk.ncells);
            blk.nic = blk.ncells; blk.njc = 1; blk.nkc = 1;
        }
        // Scan the remainder of the file, extracting our data.
        // Assume it is in standard cell order.
        blk._data.length = blk.ncells;
        foreach (i; 0 .. blk.ncells) {
            line = byLine.front; byLine.popFront();
            tokens = line.strip().split();
            assert(tokens.length == blk.variableNames.length,
                   "wrong number of items for variable data");
            blk._data[i].length = blk.variableNames.length;
            foreach (ivar; 0 .. blk.variableNames.length) {
                blk._data[i][ivar] = to!double(tokens[ivar]);
            }
        } // foreach i
        byLine.range.f.close();
    } // end switch flow_format
    foreach (i; 0 .. blk.variableNames.length) { blk.variableIndex[blk.variableNames[i]] = i; }
} // end read_block_data_from_file()

void read_extra_vars_from_file(BlockFlow blk, string filename, string[] extraVars)
{
    string[] tokens;
    auto byLine = new GzipByLine(filename);
    auto line = byLine.front; byLine.popFront();
    string format_version;
    string dummyName;
    final switch (blk.gridType) {
    case Grid_t.structured_grid:
        formattedRead(line, "%s  %s", &dummyName, &format_version);
        break;
    case Grid_t.unstructured_grid:
        formattedRead(line, "%s %s", &dummyName, &format_version);
    }
    if (format_version != "1.0") {
        string msg = text("BlockFlow.read_extra_vars_from_file(): " ~
                          "format version found: " ~ format_version);
        throw new FlowSolverException(msg);
    }
    // Throw away label line.
    byLine.popFront();
    // Throw away variables line.
    byLine.popFront();
    // Read actual list of variables.
    line = byLine.front; byLine.popFront();
    string[] varNames = line.strip().split();
    if (varNames.length != extraVars.length) {
        string msg = text("BlockFlow.read_extra_vars_from_file(): " ~
                          "number of variables in file does not match user-supplied variables list.");
        throw new FlowSolverException(msg);
    }
    foreach (var; extraVars) {
        blk.variableIndex[var] = blk.variableNames.length;
        blk.variableNames ~= var;
    }
    line = byLine.front; byLine.popFront();
    formattedRead(line, "dimensions: %d", &blk.dimensions);
    final switch (blk.gridType) {
    case Grid_t.structured_grid:
        line = byLine.front; byLine.popFront();
        formattedRead(line, "nicell: %d", &blk.nic);
        line = byLine.front; byLine.popFront();
        formattedRead(line, "njcell: %d", &blk.njc);
        line = byLine.front; byLine.popFront();
        formattedRead(line, "nkcell: %d", &blk.nkc);
        blk.ncells = blk.nic*blk.njc*blk.nkc;
        break;
    case Grid_t.unstructured_grid:
        line = byLine.front; byLine.popFront();
        formattedRead(line, "ncells: %d", &blk.ncells);
        blk.nic = blk.ncells; blk.njc = 1; blk.nkc = 1;
    }
    if (blk._data.length != blk.ncells) {
        string msg = text("BlockFlow.read_extra_vars_from_file(): " ~
                          "number of cells in file does not match number of existing cells in block.");
        throw new FlowSolverException(msg);
    }
    // Scan the remainder of the file, extracting our data.
    // Assume it is in standard cell order.
    foreach (i; 0 .. blk.ncells) {
        line = byLine.front; byLine.popFront();
        tokens = line.strip().split();
        assert(tokens.length == extraVars.length,
               "wrong number of items for variable data");
        foreach (ivar; 0 .. extraVars.length) {
            blk._data[i] ~= to!double(tokens[ivar]);
        }
    } // foreach i
    byLine.range.f.close();
} // end read_extra_vars_from_file()


//------------------------------------- Cell IO -------------------------------------------
//
void scan_values_from_string(FVCell c, string buffer, ref string[] varNameList, bool fixedOrder,
                             GasModel gmodel, bool overwrite_geometry_data)
// Note that the position data is read into grid_time_level 0.
{
    Vector3 new_pos;
    number new_volume;
    if (fixedOrder) {
        scan_cell_data_from_fixed_order_string
            (buffer,
             new_pos, new_volume, c.fs,
             c.Q_rad_org, c.f_rad_org, c.Q_rE_rad,
             c.myConfig.with_local_time_stepping, c.dt_local, c.dt_chem, c.dt_therm,
             c.myConfig.include_quality, c.myConfig.MHD,
             c.myConfig.divergence_cleaning, c.myConfig.radiation,
             c.myConfig.turb_model.nturb);
    } else {
        scan_cell_data_from_variable_order_string
            (buffer, varNameList, gmodel, c.myConfig.turb_model,
             new_pos, new_volume, c.fs,
             c.Q_rad_org, c.f_rad_org, c.Q_rE_rad,
             c.myConfig.with_local_time_stepping, c.dt_local, c.dt_chem, c.dt_therm,
             c.myConfig.include_quality, c.myConfig.MHD,
             c.myConfig.divergence_cleaning, c.myConfig.radiation);
    }
    if (overwrite_geometry_data) {
        c.pos[0].set(new_pos);
        c.volume[0] = new_volume;
    }
} // end scan_values_from_string()

void read_values_from_raw_binary(FVCell c, ref File fin, bool overwrite_geometry_data)
// Note that the position data is read into grid_time_level 0.
{
    Vector3 new_pos;
    number new_volume;
    raw_binary_to_cell_data(fin, new_pos, new_volume, c.fs,
                            c.Q_rad_org, c.f_rad_org, c.Q_rE_rad,
                            c.myConfig.with_local_time_stepping, c.dt_local, c.dt_chem, c.dt_therm,
                            c.myConfig.include_quality, c.myConfig.MHD,
                            c.myConfig.divergence_cleaning, c.myConfig.radiation,
                            c.myConfig.turb_model.nturb);
    if (overwrite_geometry_data) {
        c.pos[0].set(new_pos);
        c.volume[0] = new_volume;
    }
} // end read_values_from_raw_binary()

string write_values_to_string(const(FVCell) c)
{
    return cell_data_as_string(c.pos[0], c.volume[0], c.fs,
                               c.Q_rad_org, c.f_rad_org, c.Q_rE_rad,
                               c.myConfig.with_local_time_stepping, c.dt_local, c.dt_chem, c.dt_therm,
                               c.myConfig.include_quality, c.myConfig.MHD,
                               c.myConfig.divergence_cleaning, c.myConfig.radiation,
                               c.myConfig.turb_model.nturb);
} // end write_values_to_string()

void write_values_to_raw_binary(const(FVCell) c, ref File fout)
{
    cell_data_to_raw_binary(fout, c.pos[0], c.volume[0], c.fs,
                            c.Q_rad_org, c.f_rad_org, c.Q_rE_rad,
                            c.myConfig.with_local_time_stepping, c.dt_local, c.dt_chem, c.dt_therm,
                            c.myConfig.include_quality, c.myConfig.MHD,
                            c.myConfig.divergence_cleaning, c.myConfig.radiation,
                            c.myConfig.turb_model.nturb);
} // end write_values_to_raw_binary()

string write_residuals_to_string(const(FVCell) c)
{
    auto writer = appender!string();
    version(complex_numbers) {
        formattedWrite(writer, "%.18e %.18e %.18e %.18e",
                       -c.dUdt[0].mass.re, -c.dUdt[0].momentum.x.re,
                       -c.dUdt[0].momentum.y.re, -c.dUdt[0].total_energy.re);
    } else {
        formattedWrite(writer, "%.18e %.18e %.18e %.18e",
                       -c.dUdt[0].mass, -c.dUdt[0].momentum.x,
                       -c.dUdt[0].momentum.y, -c.dUdt[0].total_energy);
    }
    return writer.data;
}

// begin write_DFT_to_string()
string write_DFT_to_string(FVCell c)
{
    auto writer = appender!string();
    foreach(i; 0 .. c.myConfig.DFT_n_modes) {
        formattedWrite(writer, "%.18e %1.8e ", c.DFT_local_real[i], c.DFT_local_imag[i]);
    }
    return writer.data;
} // end write_DFT_to_string


string cell_data_as_string(ref const(Vector3) pos, number volume, ref const(FlowState) fs,
                           number Q_rad_org, number f_rad_org, number Q_rE_rad,
                           bool with_local_time_stepping, double dt_local, double dt_chem, double dt_therm,
                           bool include_quality,
                           bool MHDflag, bool divergence_cleaning,
                           bool radiation, size_t nturb)
{
    // We'll treat this function as the master definition of the data format.
    auto writer = appender!string();
    version(complex_numbers) {
        // For complex_numbers, we presently write out only the real parts.
        // [TODO] Maybe we should write full complex numbers.
        formattedWrite(writer, "%.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e",
                       pos.x.re, pos.y.re, pos.z.re, volume.re, fs.gas.rho.re,
                       fs.vel.x.re, fs.vel.y.re, fs.vel.z.re);
        version(MHD) {
            if (MHDflag) { formattedWrite(writer, " %.18e %.18e %.18e %.18e", fs.B.x.re, fs.B.y.re, fs.B.z.re, fs.divB.re); }
            if (MHDflag && divergence_cleaning) { formattedWrite(writer, " %.18e", fs.psi.re); }
        } else {
            assert(!MHDflag, "inappropriate MHDflag");
        }
        if (include_quality) { formattedWrite(writer, " %.18e", fs.gas.quality.re); }
        formattedWrite(writer, " %.18e %.18e %.18e", fs.gas.p.re, fs.gas.a.re, fs.gas.mu.re);
        formattedWrite(writer, " %.18e", fs.gas.k.re);
        version(multi_T_gas) {
            foreach (kvalue; fs.gas.k_modes) { formattedWrite(writer, " %.18e", kvalue.re); }
        }
        formattedWrite(writer, " %.18e %.18e %.18e", fs.mu_t.re, fs.k_t.re, fs.S.re);
        if (radiation) { formattedWrite(writer, " %.18e %.18e %.18e", Q_rad_org.re, f_rad_org.re, Q_rE_rad.re); }
        version(turbulence) {
            foreach(it; 0 .. nturb){
                formattedWrite(writer, " %.18e", fs.turb[it].re);
            }
        }
        version(multi_species_gas) {
            foreach (massfvalue; fs.gas.massf) { formattedWrite(writer, " %.18e", massfvalue.re); }
            if (fs.gas.massf.length > 1) { formattedWrite(writer, " %.18e", dt_chem); }
        } else {
            formattedWrite(writer, " %.18e", 1.0); // single-species mass fraction
        }
        formattedWrite(writer, " %.18e %.18e", fs.gas.u.re, fs.gas.T.re);
        version(multi_T_gas) {
            foreach (imode; 0 .. fs.gas.u_modes.length) {
                formattedWrite(writer, " %.18e %.18e", fs.gas.u_modes[imode].re, fs.gas.T_modes[imode].re);
            }
            if (fs.gas.u_modes.length > 0) { formattedWrite(writer, " %.18e", dt_therm); }
        }
        if (with_local_time_stepping) formattedWrite(writer, " %.18e", dt_local);
    } else {
        // version double_numbers
        formattedWrite(writer, "%.18e %.18e %.18e %.18e %.18e %.18e %.18e %.18e",
                       pos.x, pos.y, pos.z, volume, fs.gas.rho,
                       fs.vel.x, fs.vel.y, fs.vel.z);
        version(MHD) {
            if (MHDflag) { formattedWrite(writer, " %.18e %.18e %.18e %.18e", fs.B.x, fs.B.y, fs.B.z, fs.divB); }
            if (MHDflag && divergence_cleaning) { formattedWrite(writer, " %.18e", fs.psi); }
        } else {
            assert(!MHDflag, "inappropriate MHDflag");
        }
        if (include_quality) { formattedWrite(writer, " %.18e", fs.gas.quality); }
        formattedWrite(writer, " %.18e %.18e %.18e", fs.gas.p, fs.gas.a, fs.gas.mu);
        formattedWrite(writer, " %.18e", fs.gas.k);
        version(multi_T_gas) {
            foreach (kvalue; fs.gas.k_modes) { formattedWrite(writer, " %.18e", kvalue); }
        }
        formattedWrite(writer, " %.18e %.18e %.18e", fs.mu_t, fs.k_t, fs.S);
        if (radiation) { formattedWrite(writer, " %.18e %.18e %.18e", Q_rad_org, f_rad_org, Q_rE_rad); }
        version(turbulence) {
            foreach(it; 0 .. nturb){
                formattedWrite(writer, " %.18e", fs.turb[it]);
            }
        }
        version(multi_species_gas) {
            foreach (massfvalue; fs.gas.massf) { formattedWrite(writer, " %.18e", massfvalue); }
            if (fs.gas.massf.length > 1) { formattedWrite(writer, " %.18e", dt_chem); }
        } else {
            formattedWrite(writer, " %.18e", 1.0); // single-species mass fraction
        }
        formattedWrite(writer, " %.18e %.18e", fs.gas.u, fs.gas.T);
        version(multi_T_gas) {
            foreach (imode; 0 .. fs.gas.u_modes.length) {
                formattedWrite(writer, " %.18e %.18e", fs.gas.u_modes[imode], fs.gas.T_modes[imode]);
            }
            if (fs.gas.u_modes.length > 0) { formattedWrite(writer, " %.18e", dt_therm); }
        }
        if (with_local_time_stepping) { formattedWrite(writer, " %.18e", dt_local); }
    } // end version double_numbers
    return writer.data;
} // end cell_data_as_string()

void cell_data_to_raw_binary(ref File fout,
                             ref const(Vector3) pos, number volume, ref const(FlowState) fs,
                             number Q_rad_org, number f_rad_org, number Q_rE_rad,
                             bool with_local_time_stepping, double dt_local, double dt_chem, double dt_therm,
                             bool include_quality, bool MHDflag, bool divergence_cleaning,
                             bool radiation, size_t nturb)
{
    // This function should match function cell_data_as_string()
    // which is considered the master definition of the data format.
    // We have tried to keep the same code layout.  There is some history.
    // 2017-09-02:
    // Change to using all double values so that the FlowSolution reader becomes simpler.
    //
    version(complex_numbers) {
        // For complex_numbers, we presently write out only the real parts.
        // [TODO] Maybe we should write full complex numbers.
        // Fixed-length buffers to hold data for sending to binary file.
        double[1] dbl1; double[2] dbl2; double[3] dbl3; double[4] dbl4;
        //
        dbl4[0] = pos.x.re; dbl4[1] = pos.y.re; dbl4[2] = pos.z.re; dbl4[3] = volume.re;
        fout.rawWrite(dbl4);
        dbl4[0] = fs.gas.rho.re; dbl4[1] = fs.vel.x.re; dbl4[2] = fs.vel.y.re; dbl4[3] = fs.vel.z.re;
        fout.rawWrite(dbl4);
        version(MHD) {
            if (MHDflag) {
                dbl4[0] = fs.B.x.re; dbl4[1] = fs.B.y.re; dbl4[2] = fs.B.z.re; dbl4[3] = fs.divB.re;
                fout.rawWrite(dbl4);
            }
            if (MHDflag && divergence_cleaning) { dbl1[0] = fs.psi.re; fout.rawWrite(dbl1); }
        } else {
            assert(!MHDflag, "inappropriate MHDflag");
        }
        if (include_quality) { dbl1[0] = fs.gas.quality.re; fout.rawWrite(dbl1); }
        dbl4[0] = fs.gas.p.re; dbl4[1] = fs.gas.a.re; dbl4[2] = fs.gas.mu.re; dbl4[3] = fs.gas.k.re;
        fout.rawWrite(dbl4);
        foreach (kvalue; fs.gas.k_modes) { dbl1[0] = kvalue.re; fout.rawWrite(dbl1); }
        dbl2[0] = fs.mu_t.re; dbl2[1] = fs.k_t.re; fout.rawWrite(dbl2);
        dbl1[0] = fs.S.re; fout.rawWrite(dbl1);
        if (radiation) {
            dbl3[0] = Q_rad_org.re; dbl3[1] = f_rad_org.re; dbl3[2] = Q_rE_rad.re;
            fout.rawWrite(dbl3);
        }
        version(turbulence) {
            foreach(it; 0 .. nturb){
                dbl1[0] = fs.turb[it].re; fout.rawWrite(dbl1);
            }
        }
        version(multi_species_gas) {
            foreach (mf; fs.gas.massf) { dbl1[0] = mf.re; fout.rawWrite(dbl1); }
            if (fs.gas.massf.length > 1) { dbl1[0] = dt_chem; fout.rawWrite(dbl1); }
        } else {
            dbl1[0] = 1.0; fout.rawWrite(dbl1); // single-species mass fraction
        }
        dbl2[0] = fs.gas.u.re; dbl2[1] = fs.gas.T.re; fout.rawWrite(dbl2);
        version(multi_T_gas) {
            foreach (imode; 0 .. fs.gas.u_modes.length) {
                dbl2[0] = fs.gas.u_modes[imode].re; dbl2[1] = fs.gas.T_modes[imode].re;
                fout.rawWrite(dbl2);
            }
            if (fs.gas.u_modes.length > 0) { dbl1[0] = dt_therm; fout.rawWrite(dbl1); }
        }
        if (with_local_time_stepping) { dbl1[0] = dt_local; fout.rawWrite(dbl1); }
    } else {
        // version double_numbers
        // Fixed-length buffers to hold data for sending to binary file.
        double[1] dbl1; double[2] dbl2; double[3] dbl3; double[4] dbl4;
        //
        dbl4[0] = pos.x; dbl4[1] = pos.y; dbl4[2] = pos.z; dbl4[3] = volume;
        fout.rawWrite(dbl4);
        dbl4[0] = fs.gas.rho; dbl4[1] = fs.vel.x; dbl4[2] = fs.vel.y; dbl4[3] = fs.vel.z;
        fout.rawWrite(dbl4);
        version(MHD) {
            if (MHDflag) {
                dbl4[0] = fs.B.x; dbl4[1] = fs.B.y; dbl4[2] = fs.B.z; dbl4[3] = fs.divB;
                fout.rawWrite(dbl4);
            }
            if (MHDflag && divergence_cleaning) { dbl1[0] = fs.psi; fout.rawWrite(dbl1); }
        } else {
            assert(!MHDflag, "inappropriate MHDflag");
        }
        if (include_quality) { dbl1[0] = fs.gas.quality; fout.rawWrite(dbl1); }
        dbl4[0] = fs.gas.p; dbl4[1] = fs.gas.a; dbl4[2] = fs.gas.mu; dbl4[3] = fs.gas.k;
        fout.rawWrite(dbl4);
        version(multi_species_gas) {
            foreach (kvalue; fs.gas.k_modes) { dbl1[0] = kvalue; fout.rawWrite(dbl1); }
        }
        dbl2[0] = fs.mu_t; dbl2[1] = fs.k_t; fout.rawWrite(dbl2);
        dbl1[0] = fs.S; fout.rawWrite(dbl1);
        if (radiation) {
            dbl3[0] = Q_rad_org; dbl3[1] = f_rad_org; dbl3[2] = Q_rE_rad;
            fout.rawWrite(dbl3);
        }
        version(turbulence) {
            foreach(it; 0 .. nturb){
                dbl1[0] = fs.turb[it]; fout.rawWrite(dbl1);
            }
        }
        version(multi_species_gas) {
            fout.rawWrite(fs.gas.massf);
            if (fs.gas.massf.length > 1) { dbl1[0] = dt_chem; fout.rawWrite(dbl1); }
        } else {
            dbl1[0] = 1.0; fout.rawWrite(dbl1); // single-species mass fraction
        }
        dbl2[0] = fs.gas.u; dbl2[1] = fs.gas.T; fout.rawWrite(dbl2);
        version(multi_T_gas) {
            foreach (imode; 0 .. fs.gas.u_modes.length) {
                dbl2[0] = fs.gas.u_modes[imode]; dbl2[1] = fs.gas.T_modes[imode];
                fout.rawWrite(dbl2);
            }
            if (fs.gas.u_modes.length > 0) { dbl1[0] = dt_therm; fout.rawWrite(dbl1); }
        }
        if (with_local_time_stepping) { dbl1[0] = dt_local; fout.rawWrite(dbl1); }
    } // end version double_numbers
    return;
} // end cell_data_to_raw_binary()

void scan_cell_data_from_fixed_order_string
(string buffer,
 ref Vector3 pos, ref number volume, ref FlowState fs,
 ref number Q_rad_org, ref number f_rad_org, ref number Q_rE_rad,
 bool with_local_time_stepping, ref double dt_local, ref double dt_chem, ref double dt_therm,
 bool include_quality, bool MHDflag, bool divergence_cleaning, bool radiation, size_t nturb)
{
    // This function needs to be kept consistent with cell_data_as_string() above.
    auto items = split(buffer);
    version(complex_numbers) {
        // For complex_numbers, we presently set only the real parts.
        // [TODO] Maybe we should read full complex numbers.
        pos.refx = Complex!double(items.front); items.popFront();
        pos.refy = Complex!double(items.front); items.popFront();
        pos.refz = Complex!double(items.front); items.popFront();
        volume = Complex!double(items.front); items.popFront();
        fs.gas.rho = Complex!double(items.front); items.popFront();
        fs.vel.refx = Complex!double(items.front); items.popFront();
        fs.vel.refy = Complex!double(items.front); items.popFront();
        fs.vel.refz = Complex!double(items.front); items.popFront();
        version(MHD) {
            if (MHDflag) {
                fs.B.refx = Complex!double(items.front); items.popFront();
                fs.B.refy = Complex!double(items.front); items.popFront();
                fs.B.refz = Complex!double(items.front); items.popFront();
                fs.divB = Complex!double(items.front); items.popFront();
                if (divergence_cleaning) {
                    fs.psi = Complex!double(items.front); items.popFront();
                } else {
                    fs.psi = 0.0;
                }
            } else {
                fs.B.clear(); fs.psi = 0.0; fs.divB = 0.0;
            }
        }
        if (include_quality) {
            fs.gas.quality = Complex!double(items.front); items.popFront();
        } else {
            fs.gas.quality = 1.0;
        }
        fs.gas.p = Complex!double(items.front); items.popFront();
        fs.gas.a = Complex!double(items.front); items.popFront();
        fs.gas.mu = Complex!double(items.front); items.popFront();
        fs.gas.k = Complex!double(items.front); items.popFront();
        version(multi_T_gas) {
            foreach(i; 0 .. fs.gas.k_modes.length) {
                fs.gas.k_modes[i] = Complex!double(items.front); items.popFront();
            }
        }
        fs.mu_t = Complex!double(items.front); items.popFront();
        fs.k_t = Complex!double(items.front); items.popFront();
        fs.S = Complex!double(items.front); items.popFront();
        if (radiation) {
            Q_rad_org = Complex!double(items.front); items.popFront();
            f_rad_org = Complex!double(items.front); items.popFront();
            Q_rE_rad = Complex!double(items.front); items.popFront();
        } else {
            Q_rad_org = 0.0; f_rad_org = 0.0; Q_rE_rad = 0.0;
        }
        version(turbulence) {
            foreach(it; 0 .. nturb) {
                fs.turb[it] = Complex!double(items.front); items.popFront();
            }
        }
        version(multi_species_gas) {
            foreach(i; 0 .. fs.gas.massf.length) {
                fs.gas.massf[i] = Complex!double(items.front); items.popFront();
            }
            if (fs.gas.massf.length > 1) {
                dt_chem = to!double(items.front); items.popFront();
            }
        } else {
            items.popFront(); // discard the single-species mass fraction, assumed 1.0
        }
        fs.gas.u = Complex!double(items.front); items.popFront();
        fs.gas.T = Complex!double(items.front); items.popFront();
        version(multi_T_gas) {
            foreach(i; 0 .. fs.gas.u_modes.length) {
                fs.gas.u_modes[i] = Complex!double(items.front); items.popFront();
                fs.gas.T_modes[i] = Complex!double(items.front); items.popFront();
            }
            if (fs.gas.u_modes.length > 0) {
                dt_therm = to!double(items.front); items.popFront();
            }
        }
        if (with_local_time_stepping) { dt_local = to!double(items.front); items.popFront(); }
    } else {
        // version double_numbers
        pos.refx = to!double(items.front); items.popFront();
        pos.refy = to!double(items.front); items.popFront();
        pos.refz = to!double(items.front); items.popFront();
        volume = to!double(items.front); items.popFront();
        fs.gas.rho = to!double(items.front); items.popFront();
        fs.vel.refx = to!double(items.front); items.popFront();
        fs.vel.refy = to!double(items.front); items.popFront();
        fs.vel.refz = to!double(items.front); items.popFront();
        version(MHD) {
            if (MHDflag) {
                fs.B.refx = to!double(items.front); items.popFront();
                fs.B.refy = to!double(items.front); items.popFront();
                fs.B.refz = to!double(items.front); items.popFront();
                fs.divB = to!double(items.front); items.popFront();
                if (divergence_cleaning) {
                    fs.psi = to!double(items.front); items.popFront();
                } else {
                    fs.psi = 0.0;
                }
            } else {
                fs.B.clear(); fs.psi = 0.0; fs.divB = 0.0;
            }
        }
        if (include_quality) {
            fs.gas.quality = to!double(items.front); items.popFront();
        } else {
            fs.gas.quality = 1.0;
        }
        fs.gas.p = to!double(items.front); items.popFront();
        fs.gas.a = to!double(items.front); items.popFront();
        fs.gas.mu = to!double(items.front); items.popFront();
        fs.gas.k = to!double(items.front); items.popFront();
        version(multi_T_gas) {
            foreach(i; 0 .. fs.gas.k_modes.length) {
                fs.gas.k_modes[i] = to!double(items.front); items.popFront();
            }
        }
        fs.mu_t = to!double(items.front); items.popFront();
        fs.k_t = to!double(items.front); items.popFront();
        fs.S = to!double(items.front); items.popFront();
        if (radiation) {
            Q_rad_org = to!double(items.front); items.popFront();
            f_rad_org = to!double(items.front); items.popFront();
            Q_rE_rad = to!double(items.front); items.popFront();
        } else {
            Q_rad_org = 0.0; f_rad_org = 0.0; Q_rE_rad = 0.0;
        }
        version(turbulence) {
            foreach(it; 0 .. nturb) {
                fs.turb[it] = to!double(items.front); items.popFront();
            }
        }
        version(multi_species_gas) {
            foreach(i; 0 .. fs.gas.massf.length) {
                fs.gas.massf[i] = to!double(items.front); items.popFront();
            }
            if (fs.gas.massf.length > 1) {
                dt_chem = to!double(items.front); items.popFront();
            }
        } else {
            items.popFront(); // discard single-species mass fraction
        }
        fs.gas.u = to!double(items.front); items.popFront();
        fs.gas.T = to!double(items.front); items.popFront();
        version(multi_T_gas) {
            foreach(i; 0 .. fs.gas.u_modes.length) {
                fs.gas.u_modes[i] = to!double(items.front); items.popFront();
                fs.gas.T_modes[i] = to!double(items.front); items.popFront();
            }
            if (fs.gas.u_modes.length > 0) {
                dt_therm = to!double(items.front); items.popFront();
            }
        }
        if (with_local_time_stepping) { dt_local = to!double(items.front); items.popFront(); }
    } // end version double_numbers
} // end scan_values_from_fixed_order_string()

void scan_cell_data_from_variable_order_string
(string buffer, string[] varNameList, GasModel gmodel, const ref TurbulenceModel tm,
 ref Vector3 pos, ref number volume, ref FlowState fs,
 ref number Q_rad_org, ref number f_rad_org, ref number Q_rE_rad,
 bool with_local_time_stepping, ref double dt_local, ref double dt_chem, ref double dt_therm,
 bool include_quality, bool MHDflag, bool divergence_cleaning, bool radiation)
{
    // This function uses the list of variable names read from the file
    // to work out which data item to assign to each variable.
    version(complex_numbers) {
        Complex!double[] values;
        // Note that we expect only the real part of each item to be in the string.
        foreach (item; buffer.strip().split()) { values ~= to!(Complex!double)(item); }
    } else {
        double[] values;
        foreach (item; buffer.strip().split()) { values ~= to!double(item); }
    }
    pos.refx = values[countUntil(varNameList, flowVarName(FlowVar.pos_x))];
    pos.refy = values[countUntil(varNameList, flowVarName(FlowVar.pos_y))];
    pos.refz = values[countUntil(varNameList, flowVarName(FlowVar.pos_z))];
    volume = values[countUntil(varNameList, flowVarName(FlowVar.volume))];
    fs.gas.rho = values[countUntil(varNameList, flowVarName(FlowVar.rho))];
    fs.vel.refx = values[countUntil(varNameList, flowVarName(FlowVar.vel_x))];
    fs.vel.refy = values[countUntil(varNameList, flowVarName(FlowVar.vel_y))];
    fs.vel.refz = values[countUntil(varNameList, flowVarName(FlowVar.vel_z))];
    version(MHD) {
        if (MHDflag) {
            fs.B.refx = values[countUntil(varNameList, flowVarName(FlowVar.B_x))];
            fs.B.refy = values[countUntil(varNameList, flowVarName(FlowVar.B_y))];
            fs.B.refz = values[countUntil(varNameList, flowVarName(FlowVar.B_z))];
            fs.divB = values[countUntil(varNameList, flowVarName(FlowVar.divB))];
            if (divergence_cleaning) {
                fs.psi = values[countUntil(varNameList, flowVarName(FlowVar.psi))];
            } else {
                fs.psi = 0.0;
            }
        } else {
            fs.B.clear(); fs.psi = 0.0; fs.divB = 0.0;
        }
    }
    if (include_quality) {
        fs.gas.quality = values[countUntil(varNameList, flowVarName(FlowVar.quality))];
    } else {
        fs.gas.quality = 1.0;
    }
    fs.gas.p = values[countUntil(varNameList, flowVarName(FlowVar.p))];
    fs.gas.a = values[countUntil(varNameList, flowVarName(FlowVar.a))];
    fs.gas.mu = values[countUntil(varNameList, flowVarName(FlowVar.mu))];
    fs.gas.k = values[countUntil(varNameList, flowVarName(FlowVar.k))];
    version(multi_T_gas) {
        foreach(i; 0 .. fs.gas.k_modes.length) {
            fs.gas.k_modes[i] = values[countUntil(varNameList, k_modesName(to!int(i)))];
        }
    }
    fs.mu_t = values[countUntil(varNameList, flowVarName(FlowVar.mu_t))];
    fs.k_t = values[countUntil(varNameList, flowVarName(FlowVar.k_t))];
    fs.S = values[countUntil(varNameList, flowVarName(FlowVar.S))];
    if (radiation) {
        Q_rad_org = values[countUntil(varNameList, flowVarName(FlowVar.Q_rad_org))];
        f_rad_org = values[countUntil(varNameList, flowVarName(FlowVar.f_rad_org))];
        Q_rE_rad = values[countUntil(varNameList, flowVarName(FlowVar.Q_rE_rad))];
    } else {
        Q_rad_org = 0.0; f_rad_org = 0.0; Q_rE_rad = 0.0;
    }
    version(turbulence) {
        foreach(i; 0 .. tm.nturb) {
            fs.turb[i] = values[countUntil(varNameList, tm.primitive_variable_name(i))];
        }
    }
    version(multi_species_gas) {
        foreach(i; 0 .. fs.gas.massf.length) {
            fs.gas.massf[i] = values[countUntil(varNameList, massfName(gmodel, to!int(i)))];
        }
        if (fs.gas.massf.length > 1) {
            dt_chem = values[countUntil(varNameList, flowVarName(FlowVar.dt_chem))].re;
        }
    }
    fs.gas.u = values[countUntil(varNameList, flowVarName(FlowVar.u))];
    fs.gas.T = values[countUntil(varNameList, flowVarName(FlowVar.T))];
    version(multi_T_gas) {
        foreach(i; 0 .. fs.gas.u_modes.length) {
            fs.gas.u_modes[i] = values[countUntil(varNameList, u_modesName(to!int(i)))];
            fs.gas.T_modes[i] = values[countUntil(varNameList, T_modesName(to!int(i)))];
        }
        if (fs.gas.u_modes.length > 0) {
            dt_therm = values[countUntil(varNameList, flowVarName(FlowVar.dt_therm))].re;
        }
    }
    if (with_local_time_stepping) { dt_local = values[countUntil(varNameList, flowVarName(FlowVar.dt_local))].re; }
} // end scan_values_from_variable_order_string()

void raw_binary_to_cell_data(ref File fin,
                             ref Vector3 pos, ref number volume, ref FlowState fs,
                             ref number Q_rad_org, ref number f_rad_org, ref number Q_rE_rad,
                             bool with_local_time_stepping, ref double dt_local, ref double dt_chem, ref double dt_therm,
                             bool include_quality, bool MHDflag, bool divergence_cleaning,
                             bool radiation, size_t nturb)
{
    // This function needs to be kept consistent with cell_data_to_raw_binary() above.
    //
    version(complex_numbers) {
        // For complex_numbers, we presently set only the real parts.
        // [TODO] Maybe we should read full complex numbers.
        // Fixed-length buffers to hold data while reading binary file.
        double[1] dbl1; double[2] dbl2; double[3] dbl3; double[4] dbl4;
        fin.rawRead(dbl4);
        pos.set(dbl4[0], dbl4[1], dbl4[2]);
        volume = dbl4[3];
        fin.rawRead(dbl4);
        fs.gas.rho = dbl4[0];
        fs.vel.set(dbl4[1], dbl4[2], dbl4[3]);
        version(MHD) {
            if (MHDflag) {
                fin.rawRead(dbl4);
                fs.B.set(dbl4[0], dbl4[1], dbl4[2]);
                fs.divB = dbl4[3];
                if (divergence_cleaning) {
                    fin.rawRead(dbl1); fs.psi = dbl1[0];
                } else {
                    fs.psi = 0.0;
                }
            } else {
                fs.B.clear(); fs.psi = 0.0; fs.divB = 0.0;
            }
        }
        if (include_quality) {
            fin.rawRead(dbl1); fs.gas.quality = dbl1[0];
        } else {
            fs.gas.quality = 1.0;
        }
        fin.rawRead(dbl4);
        fs.gas.p = dbl4[0]; fs.gas.a = dbl4[1]; fs.gas.mu = dbl4[2]; fs.gas.k = dbl4[3];
        version(multi_T_gas) {
            foreach(i; 0 .. fs.gas.k_modes.length) {
                fin.rawRead(dbl1); fs.gas.k_modes[i] = dbl1[0];
            }
        }
        fin.rawRead(dbl2); fs.mu_t = dbl2[0]; fs.k_t = dbl2[1];
        fin.rawRead(dbl1); fs.S = dbl1[0];
        if (radiation) {
            fin.rawRead(dbl3);
            Q_rad_org = dbl3[0]; f_rad_org = dbl3[1]; Q_rE_rad = dbl3[2];
        } else {
            Q_rad_org = 0.0; f_rad_org = 0.0; Q_rE_rad = 0.0;
        }
        fin.rawRead(dbl2); // tke, omega
        version(turbulence) {
            foreach(i; 0 .. nturb){
                fin.rawRead(dbl1); fs.turb[i] = dbl1[0];
            }
        }
        version(multi_species_gas) {
            foreach (i; 0 .. fs.gas.massf.length) {
                fin.rawRead(dbl1); fs.gas.massf[i] = dbl1[0];
            }
            if (fs.gas.massf.length > 1) { fin.rawRead(dbl1); dt_chem = dbl1[0]; }
        } else {
            fin.rawRead(dbl1); // single-species mass fraction discarded
        }
        fin.rawRead(dbl2);
        fs.gas.u = dbl2[0]; fs.gas.T = dbl2[1];
        version(multi_T_gas) {
            foreach(i; 0 .. fs.gas.u_modes.length) {
                fin.rawRead(dbl2); fs.gas.u_modes[i] = dbl2[0]; fs.gas.T_modes[i] = dbl2[1];
            }
            if (fs.gas.u_modes.length > 0) { fin.rawRead(dbl1); dt_therm = dbl1[0]; }
        }
        if (with_local_time_stepping) { fin.rawRead(dbl1); dt_local = dbl1[0]; }
    } else {
        // double_numbers
        // Fixed-length buffers to hold data while reading binary file.
        double[1] dbl1; double[2] dbl2; double[3] dbl3; double[4] dbl4;
        fin.rawRead(dbl4);
        pos.set(dbl4[0], dbl4[1], dbl4[2]);
        volume = dbl4[3];
        fin.rawRead(dbl4);
        fs.gas.rho = dbl4[0];
        fs.vel.set(dbl4[1], dbl4[2], dbl4[3]);
        version(MHD) {
            if (MHDflag) {
                fin.rawRead(dbl4);
                fs.B.set(dbl4[0], dbl4[1], dbl4[2]);
                fs.divB = dbl4[3];
                if (divergence_cleaning) {
                    fin.rawRead(dbl1); fs.psi = dbl1[0];
                } else {
                    fs.psi = 0.0;
                }
            } else {
                fs.B.clear(); fs.psi = 0.0; fs.divB = 0.0;
            }
        }
        if (include_quality) {
            fin.rawRead(dbl1); fs.gas.quality = dbl1[0];
        } else {
            fs.gas.quality = 1.0;
        }
        fin.rawRead(dbl4);
        fs.gas.p = dbl4[0]; fs.gas.a = dbl4[1]; fs.gas.mu = dbl4[2]; fs.gas.k = dbl4[3];
        version(multi_T_gas) {
            foreach(i; 0 .. fs.gas.k_modes.length) {
                fin.rawRead(dbl1); fs.gas.k_modes[i] = dbl1[0];
            }
        }
        fin.rawRead(dbl2); fs.mu_t = dbl2[0]; fs.k_t = dbl2[1];
        fin.rawRead(dbl1); fs.S = dbl1[0];
        if (radiation) {
            fin.rawRead(dbl3);
            Q_rad_org = dbl3[0]; f_rad_org = dbl3[1]; Q_rE_rad = dbl3[2];
        } else {
            Q_rad_org = 0.0; f_rad_org = 0.0; Q_rE_rad = 0.0;
        }
        version(turbulence) {
            foreach(i; 0 .. nturb){
                fin.rawRead(dbl1); fs.turb[i] = dbl1[0];
            }
        }
        version(multi_species_gas) {
            fin.rawRead(fs.gas.massf);
            if (fs.gas.massf.length > 1) { fin.rawRead(dbl1); dt_chem = dbl1[0]; }
        } else {
            fin.rawRead(dbl1); // single-species mass fraction discarded
        }
        fin.rawRead(dbl2);
        fs.gas.u = dbl2[0]; fs.gas.T = dbl2[1];
        version(multi_T_gas) {
            foreach(i; 0 .. fs.gas.u_modes.length) {
                fin.rawRead(dbl2); fs.gas.u_modes[i] = dbl2[0]; fs.gas.T_modes[i] = dbl2[1];
            }
            if (fs.gas.u_modes.length > 0) { fin.rawRead(dbl1); dt_therm = dbl1[0]; }
        }
        if (with_local_time_stepping) { fin.rawRead(dbl1); dt_local = dbl1[0]; }
    } // end version double_numbers
} // end raw_binary_to_cell_data()
