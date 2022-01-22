// stream.d -- Part of the Puffin steady-flow calculator.
//
// PA Jacobs
// 2022-01-22
//
module stream;

import std.conv;
import std.stdio;
import std.string;
import std.json;
import std.format;
import std.range;
import std.math;
import std.algorithm;

import nm.schedule;
import json_helper;
import geom;
import gas;
import kinetics;
import config;
import flow;
import cell;

class StreamTube {
public:
    GasModel gmodel;
    GasState gs;
    int ncells;
    int n_bc;
    double dx_bc;
    double[] xs, y0s, y1s;
    int[] bc0s, bc1s;
    Schedule!double y_lower, y_upper;
    Schedule!int bc_lower, bc_upper;

    this(int indx, JSONValue configData)
    {
        ncells = getJSONint(configData, format("ncells_%d", indx), 0);
        gmodel = init_gas_model(Config.gas_model_file);
        auto inflowData = configData[format("inflow_%d", indx)];
        double p = getJSONdouble(inflowData, "p", 100.0e3);
        double T = getJSONdouble(inflowData, "T", 300.0);
        double[] default_massf = [1.0, ];
        foreach (i; 1 .. gmodel.n_species) { default_massf ~= 0.0; }
        double[] massf = getJSONdoublearray(inflowData, "massf", default_massf);
        gs = new GasState(gmodel);
        gs.p = p; gs.T = T; gs.massf[] = massf[];
        gmodel.update_thermo_from_pT(gs);
        //
        string fileName = format("%s/streamtube-%d.data", Config.job_name, indx);
        auto lines = File(fileName, "r").byLine();
        foreach (line; lines) {
            auto txt = line.strip();
            if (canFind(txt, "n_bc")) {
                n_bc = to!int(txt.split("=")[1]);
                continue;
            }
            if (canFind(txt, "dx_bc")) {
                dx_bc = to!double(txt.split("=")[1]);
                continue;
            }
            if (canFind(txt, "#")) continue;
            auto items = txt.split();
            if (items.length >= 5) {
                xs ~= to!double(items[0]);
                y0s ~= to!double(items[1]);
                y1s ~= to!double(items[2]);
                bc0s ~= to!int(items[3]);
                bc1s ~= to!int(items[4]);
            }
        }
        if (n_bc != xs.length) {
            writeln("WARNING: n_bc=", n_bc, " length=", xs.length);
        }
        if ((xs[1]-xs[0]-dx_bc)/dx_bc > 1.0e-9) {
            writeln("WARNING: dx=", xs[1]-xs[0], " dx_bc=", dx_bc);
        }
        y_lower = new Schedule!double(xs, y0s);
        y_upper = new Schedule!double(xs, y1s);
        bc_lower = new Schedule!int(xs, bc0s);
        bc_upper = new Schedule!int(xs, bc1s);
    } // end constructor

    override string toString()
    {
        string repr = "StreamTube(";
        repr ~= format("gmodel=%s, gs=%s", gmodel, gs);
        repr ~= format(", ncells=%d", ncells);
        repr ~= format(", n_bc=%d, dx_bc=%g", n_bc, dx_bc);
        repr ~= ")";
        return repr;
    }
} // end class StreamTube
