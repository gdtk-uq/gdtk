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

import json_helper;
import gas;
import config;
import tube;
import gasslug;
import lcell;
import piston;
import simcore;


void generate_xt_dataset(string varName, int tindxStart, int tindxEnd)
{
    writeln("Postprocessing to produce an xt-dataset.");
    return;
}


void extract_time_slice(int tindx)
{
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
        writeln("Set up and read state data for slug ", i);
        string fileName = L1dConfig.job_name ~ format("/slug-%04d-faces.data", i);
        File fp = File(fileName, "r");
        s.read_face_data(fp, tindx);
        fp.close();
        fileName = L1dConfig.job_name ~ format("/slug-%04d-cells.data", i);
        fp = File(fileName, "r");
        s.read_cell_data(fp, tindx);
        fp.close();
        writeln("Writing state data for slug ", i);
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


double[int] readTimesFile()
{
    // Returns the associative array of time values.
    double[int] times;
    string fileName = L1dConfig.job_name ~ "/times.data";
    File fp = File(fileName, "r");
    string txt = fp.readln().chomp(); // Discard header line
    double previous_time = 0.0;
    while (!fp.eof()) {
        txt = fp.readln().chomp();
        if (txt.length > 0) {
            int tindx; double tme;
            txt.formattedRead!"%d %e"(tindx, tme);
            if (tme < previous_time) {
                writeln("Warning: at tindx=%d, time=%e but previous=%e",
                        tindx, tme, previous_time);
            }
            times[tindx] = tme;
            previous_time = tme;
        }
    }
    fp.close();
    return times;
} // end readTimesFile()
