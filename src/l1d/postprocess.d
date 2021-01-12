// postprocess.d for the Lagrangian 1D Gas Dynamics, also known as L1d4.
// PA Jacobs
// 2020-04-08
//
module postprocess;

import std.conv;
import std.stdio;
import std.string;
import std.format;
import std.json;
import std.file;
import std.algorithm;
import std.math;

import json_helper;
import gas;
import config;
import tube;
import gasslug;
import lcell;
import piston;
import simcore;
import misc;


void generate_xt_dataset_gnuplot(string varName, int tindxStart, int tindxEnd, bool takeLog)
{
    writeln("Postprocessing to produce an xt-dataset for flow variable=", varName);
    double[int] times = readTimesFile();
    int[] tindices = times.keys();
    tindices.sort();
    if (!tindices.canFind(tindxStart)) { tindxStart = tindices[0]; }
    tindxEnd = min(tindxEnd, tindices[$-1]);
    size_t ntimes = 0;
    foreach (tindx; tindxStart .. tindxEnd+1) {
        if (tindices.canFind(tindx)) { ntimes += 1; }
    }
    // We need just enough of the configuration to set up the gas slug array.
    string dirName = L1dConfig.job_name;
    JSONValue jsonData = readJSONfile(dirName~"/config.json");
    auto configData = jsonData["config"];
    L1dConfig.gas_model_files = getJSONstringarray(configData, "gas_model_files", []);
    foreach (i, fileName; L1dConfig.gas_model_files) {
        auto gm = init_gas_model(fileName);
        gmodels ~= gm;
    }
    L1dConfig.nslugs = getJSONint(configData, "nslugs", 0);
    foreach (i; 0 .. L1dConfig.nslugs) {
        auto myData = jsonData[format("slug_%d", i)];
        size_t indx = gasslugs.length;
        gasslugs ~= new GasSlug(indx, myData);
    }
    // Build a GNUPlot-compatible file for each gas slug.
    foreach (i, s; gasslugs) {
        writeln("  Read state data for slug ", i);
        File fp = File(L1dConfig.job_name~format("/slug-%04d-cells.data", i), "r");
        File fpv = File(format("slug-%04d-xtdata-%s.data", i, varName), "w");
        string header = format("# x  t  ");
        header ~= (takeLog) ? format("log10(%s)", varName) : varName;
        fpv.writeln(header);
        // Add some metadata as GNUPlot comments.
        // Although not needed for GNUPlot, these may aid other programs.
        fpv.writefln("# ntimes %d", ntimes);
        fpv.writefln("# ncells %d", s.cells.length);
        // Write the data proper.
        foreach (tindx; tindxStart .. tindxEnd+1) {
            if (tindices.canFind(tindx)) {
                writeln("  Get data for tindx=", tindx);
                s.read_cell_data(fp, tindx);
                foreach (c; s.cells) {
                    double v;
                    switch (varName) {
                    case "p": v = c.gas.p; break;
                    case "T": v = c.gas.T; break;
                    case "rho": v = c.gas.rho; break;
                    case "vel": v = c.vel; break;
                    default: v = 0.0;
                    }
                    if (takeLog) { v = log10(v); }
                    fpv.writefln("%e %e %e", c.xmid, times[tindx], v);
                } // foreach c
                fpv.writeln(""); // blank line after block of data
            }
        } // foreach tindx
        fpv.close();
        fp.close();
    } // foreach gasslug
    return;
} // end generate_xt_dataset_gnuplot()


void generate_xt_dataset_vtk(int tindxStart, int tindxEnd, bool milliSec)
{
    writeln("Produce an xt-dataset in legacy VTK format for all flow variables");
    double[int] times = readTimesFile();
    int[] tindices = times.keys();
    tindices.sort();
    if (!tindices.canFind(tindxStart)) { tindxStart = tindices[0]; }
    tindxEnd = min(tindxEnd, tindices[$-1]);
    size_t ntimes = 0;
    foreach (tindx; tindxStart .. tindxEnd+1) {
        if (tindices.canFind(tindx)) { ntimes += 1; }
    }
    // We need just enough of the configuration to set up the gas slug array.
    string dirName = L1dConfig.job_name;
    JSONValue jsonData = readJSONfile(dirName~"/config.json");
    auto configData = jsonData["config"];
    L1dConfig.gas_model_files = getJSONstringarray(configData, "gas_model_files", []);
    foreach (i, fileName; L1dConfig.gas_model_files) {
        auto gm = init_gas_model(fileName);
        gmodels ~= gm;
    }
    L1dConfig.nslugs = getJSONint(configData, "nslugs", 0);
    foreach (i; 0 .. L1dConfig.nslugs) {
        auto myData = jsonData[format("slug_%d", i)];
        size_t indx = gasslugs.length;
        gasslugs ~= new GasSlug(indx, myData);
    }
    // Build a VTK legacy data file for each gas slug.
    // We are going to assume that the data can be accumulated in a structured-grid format
    // such that there are constant number of cells in the slug at each time.
    foreach (i, s; gasslugs) {
        writeln("  Read state data for slug ", i);
        size_t ncells = 0;
        File fp = File(L1dConfig.job_name~format("/slug-%04d-cells.data", i), "r");
        string[] varNames = ["x", "t", "p", "T", "rho", "vel"];
        double[][][string] data;
        foreach (v; varNames) { data[v].length = ntimes; }
        size_t jt = 0;
        foreach (tindx; tindxStart .. tindxEnd+1) {
            if (tindices.canFind(tindx)) {
                writeln("  Get data for tindx=", tindx);
                s.read_cell_data(fp, tindx);
                if (ncells == 0) {
                    ncells = s.cells.length;
                } else {
                    if (ncells != s.cells.length) {
                        throw new Error(format("Number of cells has changed: was %d now %d.",
                                               ncells, s.cells.length));
                    }
                }
                foreach (v; varNames) { data[v][jt].length = ncells; }
                foreach (ix, c; s.cells) {
                    data["x"][jt][ix] = c.xmid;
                    data["t"][jt][ix] = (milliSec) ? times[tindx]*1000.0 : times[tindx];
                    data["p"][jt][ix] = c.gas.p;
                    data["T"][jt][ix] = c.gas.T;
                    data["rho"][jt][ix] = c.gas.rho;
                    data["vel"][jt][ix] = c.vel;
                } // foreach c
                jt += 1;
            }
        } // foreach tindx
        fp.close();
        File fpv = File(format("slug-%04d-xtdata.vtk", i), "w");
        fpv.writeln("# vtk DataFile Version 2.0");
        fpv.writefln("# job: %s slug: %d time-units: %s",
                     L1dConfig.job_name, i, (milliSec)?"ms":"s");
        fpv.writefln("ASCII");
        fpv.writefln("DATASET STRUCTURED_GRID");
        fpv.writefln("DIMENSIONS %d %d 1", ncells, ntimes);
        fpv.writefln("POINTS %d double", ncells*ntimes);
        foreach (nt; 0 .. ntimes) {
            foreach (ix; 0 .. ncells) {
                fpv.writefln("%g %g 0.0", data["x"][nt][ix], data["t"][nt][ix]);
            }
        }
        fpv.writefln("POINT_DATA %d", ncells*ntimes);
        foreach (v; varNames) {
            fpv.writefln("SCALARS %s double 1", v);
            fpv.writeln("LOOKUP_TABLE default");
            foreach (nt; 0 .. ntimes) {
                foreach (ix; 0 .. ncells) {
                    fpv.writefln("%g", data[v][nt][ix]);
                }
            }
        } // end foreach v
        fpv.close();
    } // foreach gasslug
    return;
} // end generate_xt_dataset_vtk()


void extract_time_slice(int tindx)
{
    double[int] times = readTimesFile();
    int[] tindices = times.keys();
    tindices.sort();
    tindx = min(max(tindx, tindices[0]), tindices[$-1]);
    writeln("Postprocessing to extract slug data at a tindx=", tindx);
    // We need just enough of the configuration to set up the gas slug array.
    string dirName = L1dConfig.job_name;
    JSONValue jsonData = readJSONfile(dirName~"/config.json");
    auto configData = jsonData["config"];
    L1dConfig.gas_model_files = getJSONstringarray(configData, "gas_model_files", []);
    foreach (i, fileName; L1dConfig.gas_model_files) {
        auto gm = init_gas_model(fileName);
        gmodels ~= gm;
    }
    L1dConfig.nslugs = getJSONint(configData, "nslugs", 0);
    foreach (i; 0 .. L1dConfig.nslugs) {
        auto myData = jsonData[format("slug_%d", i)];
        size_t indx = gasslugs.length;
        gasslugs ~= new GasSlug(indx, myData);
    }
    foreach (i, s; gasslugs) {
        writeln("  Read state data for slug ", i);
        string fileName = L1dConfig.job_name ~ format("/slug-%04d-faces.data", i);
        File fp = File(fileName, "r");
        s.read_face_data(fp, tindx);
        fp.close();
        fileName = L1dConfig.job_name ~ format("/slug-%04d-cells.data", i);
        fp = File(fileName, "r");
        s.read_cell_data(fp, tindx);
        fp.close();
        writeln("  Writing state data for slug ", i);
        fileName = format("slug-%04d-tindx-%04d-faces.data", i, tindx);
        fp = File(fileName, "w");
        s.write_face_data(fp, tindx, true);
        fp.close();
        fileName = format("slug-%04d-tindx-%04d-cells.data", i, tindx);
        fp = File(fileName, "w");
        s.write_cell_data(fp, tindx, true);
        fp.close();
    }
    return;
} // end extract_time_slice()


void assemble_piston_history(int pindx)
{
    // We are going to insert the values for times into the record.
    writeln("Postprocessing to assemble history for piston=", pindx);
    double[int] times = readTimesFile();
    string fileName = L1dConfig.job_name ~ format("/piston-%04d.data", pindx);
    File fp = File(fileName, "r");
    fileName = format("piston-%04d-history.data", pindx);
    File fph = File(fileName, "w");
    fph.writeln("# tindx  t  x  vel  is_restrain  brakes_on  hit_buffer");
    string txt = fp.readln().chomp(); // Discard header line
    while (!fp.eof()) {
        txt = fp.readln().chomp();
        if (txt.length > 0) {
            int tindx; double x, vel; int is_restrain, brakes_on, hit_buffer;
            txt.formattedRead!"%d %e %e %d %d %d"(tindx, x, vel, is_restrain,
                                                  brakes_on, hit_buffer);
            fph.writefln("%d %e %e %e %d %d %d", tindx, times[tindx], x, vel,
                         is_restrain, brakes_on, hit_buffer);
        }
    }
    fp.close();
    fph.close();
    return;
} // end assemble_piston_history()

void trim_solution_files(int tindxEnd)
{
    writeln("Postprocessing to remove solution data after tindx= ", tindxEnd);
    double[int] times = readTimesFile();
    int[] tindices = times.keys();
    tindices.sort();
    if (tindxEnd >= tindices[$-1]) {
        writeln("Nothing to do since final tindx= ", tindices[$-1]);
        return;
    }
    double timeEnd = times[tindxEnd];
    writeln("Corresponding end time= ", timeEnd);
    //
    // We use just enough of the configuration to know which slug,
    // diaphragm, piston and history files to work on.
    string dirName = L1dConfig.job_name;
    JSONValue jsonData = readJSONfile(dirName~"/config.json");
    auto configData = jsonData["config"];
    //
    string fileName;
    string backupFileName;
    File fp_dest;
    File fp_src;
    string txt;
    fileName = L1dConfig.job_name ~ "/times.data.backup";
    if (std.file.exists(fileName)) {
        writeln("Can see a times.data.backup file already.");
        writeln("  Quitting without trimming files.");
        return;
    }
    //
    int nslugs = getJSONint(configData, "nslugs", 0);
    foreach (i; 0 .. nslugs) {
        writeln("  Trim faces file for slug ", i);
        fileName = L1dConfig.job_name ~ format("/slug-%04d-faces.data", i);
        backupFileName = fileName ~ ".backup";
        std.file.rename(fileName, backupFileName);
        fp_src = File(backupFileName, "r");
        fp_dest = File(fileName, "w");
        txt = fp_src.readln().chomp(); // header line
        fp_dest.writeln(txt);
        while (!fp_src.eof()) {
            txt = fp_src.readln().chomp();
            if (txt.length > 0 && txt.canFind("tindx")) {
                int tindx = to!int(txt.split()[2]);
                if (tindx > tindxEnd) { break; }
            }
            fp_dest.writeln(txt);
        }
        fp_src.close();
        fp_dest.close();
        //
        writeln("  Trim cells file for slug ", i);
        fileName = L1dConfig.job_name ~ format("/slug-%04d-cells.data", i);
        backupFileName = fileName ~ ".backup";
        std.file.rename(fileName, backupFileName);
        fp_src = File(backupFileName, "r");
        fp_dest = File(fileName, "w");
        txt = fp_src.readln().chomp(); // header line
        fp_dest.writeln(txt);
        while (!fp_src.eof()) {
            txt = fp_src.readln().chomp();
            if (txt.length > 0 && txt.canFind("tindx")) {
                int tindx = to!int(txt.split()[2]);
                if (tindx > tindxEnd) { break; }
            }
            fp_dest.writeln(txt);
        }
        fp_src.close();
        fp_dest.close();
    }
    int npistons = getJSONint(configData, "npistons", 0);
    foreach (i; 0 .. npistons) {
        writeln("  Trim file for piston ", i);
        fileName = L1dConfig.job_name ~ format("/piston-%04d.data", i);
        backupFileName = fileName ~ ".backup";
        std.file.rename(fileName, backupFileName);
        fp_src = File(backupFileName, "r");
        fp_dest = File(fileName, "w");
        txt = fp_src.readln().chomp(); // header line
        fp_dest.writeln(txt);
        while (!fp_src.eof()) {
            txt = fp_src.readln().chomp();
            if (txt.length > 0) {
                int tindx = to!int(txt.split()[0]);
                if (tindx > tindxEnd) { break; }
            }
            fp_dest.writeln(txt);
        }
        fp_src.close();
        fp_dest.close();
    }
    int ndiaphragms = getJSONint(configData, "ndiaphragms", 0);
    foreach (i; 0 .. ndiaphragms) {
        writeln("  Trim file for diaphragm ", i);
        fileName = L1dConfig.job_name ~ format("/diaphragm-%04d.data", i);
        backupFileName = fileName ~ ".backup";
        std.file.rename(fileName, backupFileName);
        fp_src = File(backupFileName, "r");
        fp_dest = File(fileName, "w");
        txt = fp_src.readln().chomp(); // header line
        fp_dest.writeln(txt);
        while (!fp_src.eof()) {
            txt = fp_src.readln().chomp();
            if (txt.length > 0) {
                int tindx = to!int(txt.split()[0]);
                if (tindx > tindxEnd) { break; }
            }
            fp_dest.writeln(txt);
        }
        fp_src.close();
        fp_dest.close();
    }
    int hloc_n = getJSONint(configData, "hloc_n", 0);
    foreach (i; 0 .. hloc_n) {
        writeln("  Trim history file for location ", i);
        fileName = L1dConfig.job_name ~ format("/history-loc-%04d.data", i);
        backupFileName = fileName ~ ".backup";
        std.file.rename(fileName, backupFileName);
        fp_src = File(backupFileName, "r");
        fp_dest = File(fileName, "w");
        txt = fp_src.readln().chomp(); // header line
        fp_dest.writeln(txt);
        while (!fp_src.eof()) {
            txt = fp_src.readln().chomp();
            if (txt.length > 0) {
                double tme = to!double(txt.split()[0]);
                if (tme > timeEnd) { break; }
            }
            fp_dest.writeln(txt);
        }
        fp_src.close();
        fp_dest.close();
    }
    writeln("  Trim times file.");
    fileName = L1dConfig.job_name ~ "/times.data";
    backupFileName = fileName ~ ".backup";
    std.file.rename(fileName, backupFileName);
    fp_src = File(backupFileName, "r");
    fp_dest = File(fileName, "w");
    txt = fp_src.readln().chomp(); // header line
    fp_dest.writeln(txt);
    while (!fp_src.eof()) {
        txt = fp_src.readln().chomp();
        if (txt.length > 0) {
            int tindx = to!int(txt.split()[0]);
            if (tindx > tindxEnd) { break; }
        }
        fp_dest.writeln(txt);
    }
    fp_src.close();
    fp_dest.close();
    //
    writeln("  Trim energies file.");
    fileName = L1dConfig.job_name ~ "/energies.data";
    backupFileName = fileName ~ ".backup";
    std.file.rename(fileName, backupFileName);
    fp_src = File(backupFileName, "r");
    fp_dest = File(fileName, "w");
    txt = fp_src.readln().chomp(); // header line
    fp_dest.writeln(txt);
    while (!fp_src.eof()) {
        txt = fp_src.readln().chomp();
        if (txt.length > 0) {
            double tme = to!double(txt.split()[0]);
            if (tme > timeEnd) { break; }
        }
        fp_dest.writeln(txt);
    }
    fp_src.close();
    fp_dest.close();
    return;
} // end trim_solution_files()
