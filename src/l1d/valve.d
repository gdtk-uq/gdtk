// valve.d
// The information associated with a valve.
// PA Jacobs
// 2020-08-14

module valve;

import std.conv;
import std.stdio;
import std.file;
import std.string;
import std.json;
import std.format;

import nm.schedule;
import util.json_helper;
import config;

class Valve {
public:
    size_t indx;
    string label = "";
    double x;
    Schedule!double fopen_schedule;

    this(size_t indx, JSONValue jsonData)
    {
        if (L1dConfig.verbosity_level >= 3) {
            writeln("construct valve[", indx, "] from json=", jsonData);
        }
        this.indx = indx;
        label = getJSONstring(jsonData, "label", "");
        x = getJSONdouble(jsonData, "x", 0.0);
        double[] dummy_times = [0.0,];
        double[] times_values = getJSONdoublearray(jsonData, "times", dummy_times);
        double[] dummy_fopen = [1.0,];
        double[] fopen_values = getJSONdoublearray(jsonData, "fopen", dummy_fopen);
        fopen_schedule = new Schedule!double(times_values, fopen_values);
        if (L1dConfig.verbosity_level >= 1) {
            writeln("Valve[", indx, "]:");
            writeln("  x= ", x);
            writeln("  times= ", times_values);
            writeln("  fopen= ", fopen_values);
        }
        return;
    }

    @nogc
    double fopen(double t)
    {
        return fopen_schedule.interpolate_value(t);
    }
} // end class Valve
