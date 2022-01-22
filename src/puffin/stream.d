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

    this(int indx, JSONValue configData)
    {
        gmodel = init_gas_model(Config.gas_model_file);
        auto myData = configData[format("inflow_%d", indx)];
        double p = getJSONdouble(myData, "p", 100.0e3);
        double T = getJSONdouble(myData, "T", 300.0);
        double[] default_massf = [1.0, ];
        foreach (i; 1 .. gmodel.n_species) { default_massf ~= 0.0; }
        double[] massf = getJSONdoublearray(myData, "massf", default_massf);
        gs = new GasState(gmodel);
        gs.p = p; gs.T = T; gs.massf[] = massf[];
        gmodel.update_thermo_from_pT(gs);
        //
        writeln("TODO -- read boundary data");
    }

    override string toString()
    {
        string repr = "StreamTube(";
        repr ~= format("gmodel=%s, gs=%s", gmodel, gs);
        repr ~= ")";
        return repr;
    }
} // end class StreamTube
