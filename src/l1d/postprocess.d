// postprocess.d for the Lagrangian 1D Gas Dynamics, also known as L1d4.
// PA Jacobs
// 2020-04-08
// 2022-12-08 Add cell history for David Mee.
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
import core.stdc.math: HUGE_VAL;

import util.json_helper;
import gas;
import config;
import tube;
import gasslug;
import lcell;
import piston;
import simcore;
import misc;


void generate_cell_history(double x0, int tindxStart, int tindxEnd, bool milliSec)
{
    writeln("Postprocessing to produce the history for a cell, starting at x=", x0);
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
    assert(gmodels.length == 0, "Gas models array is not empty.");
    foreach (fileName; L1dConfig.gas_model_files) { gmodels ~= init_gas_model(fileName); }
    L1dConfig.nslugs = getJSONint(configData, "nslugs", 0);
    assert(gasslugs.length == 0, "Gas slugs array is not empty.");
    foreach (i; 0 .. L1dConfig.nslugs) {
        auto myData = jsonData[format("slug_%d", i)];
        gasslugs ~= new GasSlug(i, myData);
    }
    // Look through all of the gas slugs and pick the cell closest to the pin.
    size_t selected_slug;
    size_t selected_cell;
    double near_dist = HUGE_VAL;
    foreach (islug, s; gasslugs) {
        writeln("Read state data for slug ", islug);
        File fp = File(L1dConfig.job_name~format("/slug-%04d-cells.data", islug), "r");
        int tindx = tindxStart;
        if (tindices.canFind(tindx)) {
            writeln("Get data for tindx=", tindx);
            s.read_cell_data(fp, tindx);
        } else {
            throw new Error(format("Specified tindx=%d is not available.", tindx));
        }
        fp.close();
        foreach (icell, c; s.cells) {
            double dist = fabs(c.xmid - x0);
            if (dist < near_dist) {
                selected_slug = islug;
                selected_cell = icell;
                near_dist = dist;
            }
        }
        writefln("Selected slug=%d cell=%d at initial-position=%g",
                 selected_slug, selected_cell,
                 gasslugs[selected_slug].cells[selected_cell].xmid);
    } // foreach gasslug
    //
    // Re-open the individual slug and write the data for just the selected cell.
    string fileName = format("slug-%04d-cell-%04d.data", selected_slug, selected_cell);
    File fp = File(L1dConfig.job_name~format("/slug-%04d-cells.data", selected_slug), "r");
    File fpv = File(fileName, "w");
    string header = "# t" ~ ((milliSec) ? "(ms)" : "(s)") ~ " xmid p T rho vel";
    fpv.writeln(header);
    foreach (tindx; tindxStart .. tindxEnd+1) {
        if (tindices.canFind(tindx)) {
            writeln("  Get data for tindx=", tindx);
            auto s = gasslugs[selected_slug];
            s.read_cell_data(fp, tindx);
            auto c = s.cells[selected_cell];
            double t = times[tindx]; if (milliSec) { t *= 1000.0; }
            fpv.writefln("%e %e %e %e %e %e", t, c.xmid, c.gas.p, c.gas.T, c.gas.rho, c.vel);
        }
    } // foreach tindx
    fpv.close();
    fp.close();
    writeln("Have written cell history to file: ", fileName);
    return;
}


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
    assert(gmodels.length == 0, "Gas models array is not empty.");
    foreach (fileName; L1dConfig.gas_model_files) { gmodels ~= init_gas_model(fileName); }
    L1dConfig.nslugs = getJSONint(configData, "nslugs", 0);
    assert(gasslugs.length == 0, "Gas slugs array is not empty.");
    foreach (i; 0 .. L1dConfig.nslugs) {
        auto myData = jsonData[format("slug_%d", i)];
        gasslugs ~= new GasSlug(i, myData);
    }
    // Build a GNUPlot-compatible file for each gas slug.
    foreach (islug, s; gasslugs) {
        writeln("  Read state data for slug ", islug);
        File fp = File(L1dConfig.job_name~format("/slug-%04d-cells.data", islug), "r");
        File fpv = File(format("slug-%04d-xtdata-%s.data", islug, varName), "w");
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


void generate_xt_dataset(int tindxStart, int tindxEnd, bool milliSec, string fmt)
{
    writefln("Produce an xt-dataset in %s format for all flow variables", fmt);
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
    assert(gmodels.length == 0, "Gas models array is not empty.");
    foreach (fileName; L1dConfig.gas_model_files) { gmodels ~= init_gas_model(fileName); }
    L1dConfig.nslugs = getJSONint(configData, "nslugs", 0);
    assert(gasslugs.length == 0, "Gas slugs array is not empty.");
    foreach (i; 0 .. L1dConfig.nslugs) {
        auto myData = jsonData[format("slug_%d", i)];
        gasslugs ~= new GasSlug(i, myData);
    }
    // Build an xt-data file for each gas slug.
    // We are going to assume that the data can be accumulated in a structured-grid format
    // such that there are constant number of cells in the slug at each time.
    foreach (islug, s; gasslugs) {
        writeln("  Read state data for slug ", islug);
        size_t ncells = 0;
        File fpcells = File(L1dConfig.job_name~format("/slug-%04d-cells.data", islug), "r");
        File fpfaces = File(L1dConfig.job_name~format("/slug-%04d-faces.data", islug), "r");
        string[] varNames = ["x", "t", "p", "T", "rho", "vel"];
        string[] varUnits = ["m", (milliSec)?"ms":"s", "Pa", "K", "kg/m^3", "m/s"];
        assert(varNames.length == varUnits.length, "Mismatch in variable names and units.");
        double[][][string] data;
        foreach (v; varNames) { data[v].length = ntimes; }
        double[] simTimes; simTimes.length = ntimes;
        double[] xL; xL.length = ntimes;
        double[] xR; xR.length = ntimes;
        //
        size_t jt = 0;
        foreach (tindx; tindxStart .. tindxEnd+1) {
            if (tindices.canFind(tindx)) {
                writeln("  Get data for tindx=", tindx);
                double simTime = (milliSec) ? times[tindx]*1000.0 : times[tindx];
                //
                s.read_cell_data(fpcells, tindx);
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
                    data["t"][jt][ix] = simTime;
                    data["p"][jt][ix] = c.gas.p;
                    data["T"][jt][ix] = c.gas.T;
                    data["rho"][jt][ix] = c.gas.rho;
                    data["vel"][jt][ix] = c.vel;
                } // foreach c
                //
                s.read_face_data(fpfaces, tindx);
                xL[jt] = s.faces[0].x;
                xR[jt] = s.faces[$-1].x;
                simTimes[jt] = simTime;
                //
                jt += 1;
            }
        } // foreach tindx
        fpcells.close();
        fpfaces.close();
        //
        switch (fmt) {
        case "VTK":
            File fpv = File(format("slug-%04d-xtdata.vtk", islug), "w");
            fpv.writeln("# vtk DataFile Version 2.0");
            fpv.writefln("# job: %s slug: %d time-units: %s",
                         L1dConfig.job_name, islug, (milliSec)?"ms":"s");
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
            foreach (i, v; varNames) {
                fpv.writefln("SCALARS %s-%s double 1", v, varUnits[i]);
                fpv.writeln("LOOKUP_TABLE default");
                foreach (nt; 0 .. ntimes) {
                    foreach (ix; 0 .. ncells) {
                        fpv.writefln("%g", data[v][nt][ix]);
                    }
                }
            } // end foreach v
            fpv.close();
            break;
        case "JSON":
            File fpj = File(format("slug-%04d-xtdata.json", islug), "w");
            fpj.writeln("{");
            fpj.writefln("\"jobName\": \"%s\",", L1dConfig.job_name);
            fpj.writefln("\"slug\": %d,", islug);
            fpj.writefln("\"ntimes\": %d,", ntimes);
            fpj.writefln("\"ncells\": %d,", ncells);
            fpj.write("\"varNames\": [");
            foreach (i, v; varNames) { fpj.writef("\"%s\"%s", v, (i<(varNames.length-1))?", ":""); }
            fpj.writeln("],");
            fpj.write("\"varUnits\": [");
            foreach (i, v; varUnits) { fpj.writef("\"%s\"%s", v, (i<(varUnits.length-1))?", ":""); }
            fpj.writeln("],");
            foreach (i, v; varNames) {
                fpj.writefln("\"%s\": [", v);
                foreach (nt; 0 .. ntimes) {
                    fpj.write("  [");
                    foreach (ix; 0 .. ncells) { fpj.writef("%g%s", data[v][nt][ix], (ix<(ncells-1))?", ":""); }
                    fpj.writefln("]%s", (nt<(ntimes-1))?", ":"");
                }
                fpj.writefln("],");
            }
            fpj.write("\"xL\": [");
            foreach (i, v; xL) { fpj.writef("%g%s", v, (i<(xL.length-1))?", ":""); }
            fpj.writeln("],");
            fpj.write("\"xR\": [");
            foreach (i, v; xR) { fpj.writef("%g%s", v, (i<(xR.length-1))?", ":""); }
            fpj.writeln("],");
            fpj.write("\"simTimes\": [");
            foreach (i, v; simTimes) { fpj.writef("%g%s", v, (i<(simTimes.length-1))?", ":""); }
            fpj.writeln("],");
            fpj.writeln("\"dummy\": 0");
            fpj.writeln("}");
            fpj.close();
            break;
        default:
            throw new Error(format("Invalid format specified: %s", fmt));
        }
    } // foreach gasslug
    return;
} // end generate_xt_dataset()


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
    fph.writeln("# tindx  t  x  vel  is_restrain  brakes_on  on_buffer");
    string txt = fp.readln().chomp(); // Discard header line
    while (!fp.eof()) {
        txt = fp.readln().chomp();
        if (txt.length > 0) {
            int tindx; double x, vel; int is_restrain, brakes_on, on_buffer;
            txt.formattedRead!"%d %e %e %d %d %d"(tindx, x, vel, is_restrain,
                                                  brakes_on, on_buffer);
            fph.writefln("%d %e %e %e %d %d %d", tindx, times[tindx], x, vel,
                         is_restrain, brakes_on, on_buffer);
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
