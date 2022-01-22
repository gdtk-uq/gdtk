// calc.d
//
// PA Jacobs
// 2022-01-22
//
module calc;

import std.conv;
import std.stdio;
import std.string;
import std.json;
import std.file;
import std.datetime;
import std.format;
import std.range;
import std.math;
import std.algorithm;

import json_helper;
import geom;
import config;
import stream;

__gshared static StreamTube[] streams;

struct ProgressData {
    int step = 0;
    double x = 0.0;
    int steps_since_last_plot_write = 0;
    SysTime wall_clock_start;
}

__gshared static ProgressData progress;

int init_calculation()
{
    string dirName = Config.job_name;
    JSONValue configData = readJSONfile(dirName~"/config.json");
    Config.title = getJSONstring(configData, "title", "");
    Config.gas_model_file = getJSONstring(configData, "gas_model_file", "");
    Config.reaction_file_1 = getJSONstring(configData, "reaction_files_1", "");
    Config.reaction_file_2 = getJSONstring(configData, "reaction_file_2", "");
    Config.reacting = getJSONbool(configData, "reacting", false);
    Config.T_frozen = getJSONdouble(configData, "T_frozen", 300.0);
    Config.axisymmetric = getJSONbool(configData, "axisymmetric", false);
    Config.max_x = getJSONdouble(configData, "max_x", 0.0);
    Config.max_step = getJSONint(configData, "max_step", 0);
    Config.dx = getJSONdouble(configData, "dx", 0.0);
    Config.cfl = getJSONdouble(configData, "cfl", 0.5);
    Config.cfl_count = getJSONint(configData, "cfl_count", 10);
    Config.print_count = getJSONint(configData, "print_count", 50);
    Config.plot_count = getJSONint(configData, "plot_count", 10);
    Config.x_order = getJSONint(configData, "x_order", 2);
    Config.n_streams = getJSONint(configData, "n_streams", 1);
    if (Config.verbosity_level >= 1) {
        writeln("Config:");
        writefln("  title= \"%s\"", Config.title);
        writeln("  gas_model_files= ", Config.gas_model_file);
        writeln("  reaction_files_1= ", Config.reaction_file_1);
        writeln("  reaction_files_2= ", Config.reaction_file_2);
        writeln("  reacting= ", Config.reacting);
        writeln("  T_frozen= ", Config.T_frozen);
        writeln("  axisymmetric= ", Config.axisymmetric);
        writeln("  max_x= ", Config.max_x);
        writeln("  max_step= ", Config.max_step);
        writeln("  dx= ", Config.dx);
        writeln("  cfl= ", Config.cfl);
        writeln("  cfl_count= ", Config.cfl_count);
        writeln("  print_count= ", Config.print_count);
        writeln("  plot_count= ", Config.plot_count);
        writeln("  x_order= ", Config.x_order);
        writeln("  n_streams= ", Config.n_streams);
    }
    foreach (i; 0 .. Config.n_streams) {
        streams ~= new StreamTube(i, configData);
        if (Config.verbosity_level >= 1) {
            writefln("  stream[%d]= %s", i, streams[i]);
        }
    }
    return 0;
}

int do_calculation()
{
    foreach (st; streams) {
        st.set_up_data_storage();
        st.set_up_inflow_boundary();
    }
    progress.step = 0;
    progress.x = 0.0;
    bool finished = false;
    while (!finished) {
        progress.step += 1;
        progress.x += Config.dx;
        foreach (st; streams) { st.set_up_slice(progress.x); }
        writeln("[TODO] compute to steady state.");
        writeln("[TODO] write flow data occasionally.");
        foreach (st; streams) { st.shuffle_data_west(); }
        finished = (progress.x > Config.max_x) || (progress.step >= Config.max_step);
    }
    writeln("TODO do_calculation properly");
    return 0;
}
