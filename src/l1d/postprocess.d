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


void generate_xt_dataset(string varName, int tindxStart, int tindxEnd, bool takeLog)
{
    writeln("Postprocessing to produce an xt-dataset for flow variable=", varName);
    double[int] times = readTimesFile();
    int[] tindices = times.keys();
    tindices.sort();
    if (!tindices.canFind(tindxStart)) { tindxStart = tindices[0]; }
    tindxEnd = min(tindxEnd, tindices[$-1]);
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
        if (takeLog) {
            header ~= format("log10(%s)", varName);
        } else {
            header ~= varName;
        }
        fpv.writeln(header);
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
        fp.close();
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
