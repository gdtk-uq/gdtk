// piston.d for the Lagrangian 1D Gas Dynamics, also known as L1d4.
// PA Jacobs
// 2020-04-08
//
module piston;

import std.conv;
import std.stdio;
import std.string;
import std.json;

import json_helper;
import geom;
import gas;
import gasflow;
import config;
import endcondition;


class Piston {
public:
    size_t indx;
    string label = "";
    double mass;  // mass, kg
    double diam;  // diameter, m
    double L;     // length, m
    double x;     // position, m
    double vel;   // velocity, m/s
    double front_seal_f; // friction factor
    double front_seal_area; // area over which pressure acts
    double back_seal_f;
    double back_seal_area;
    double p_restrain;
    bool is_restrain;
    bool with_brakes;
    bool brakes_on;
    double x_buffer;
    bool hit_buffer;
    int ecL_id;
    int ecR_id;
    EndCondition ecL;
    EndCondition ecR;

    this(size_t indx, JSONValue jsonData)
    {
        if (L1dConfig.verbosity_level >= 3) {
            writeln("construct piston[", indx, "] from json=", jsonData);
        }
        this.indx = indx;
        label = getJSONstring(jsonData, "label", "");
        mass = getJSONdouble(jsonData, "mass", 0.0);
        diam = getJSONdouble(jsonData, "diam", 0.0);
        L = getJSONdouble(jsonData, "length", 0.0);
        front_seal_f = getJSONdouble(jsonData, "front_seal_f", 0.0);
        front_seal_area = getJSONdouble(jsonData, "front_seal_area", 0.0);
        back_seal_f = getJSONdouble(jsonData, "back_seal_f", 0.0);
        back_seal_area = getJSONdouble(jsonData, "back_seal_area", 0.0);
        p_restrain = getJSONdouble(jsonData, "p_restrain", 0.0);
        x_buffer = getJSONdouble(jsonData, "x_buffer", 0.0);
        with_brakes = getJSONbool(jsonData, "with_brakes", false);
        ecL_id = getJSONint(jsonData, "ecL_id", -1);
        ecR_id = getJSONint(jsonData, "ecR_id", -1);
        if (L1dConfig.verbosity_level >= 1) {
            writeln("Piston[", indx, "]:");
            writeln("  mass= ", mass);
            writeln("  diam= ", diam);
            writeln("  L= ", L);
            writeln("  front_seal_f= ", front_seal_f);
            writeln("  front_seal_area= ", front_seal_area);
            writeln("  back_seal_f= ", back_seal_f);
            writeln("  back_seal_area= ", back_seal_area);
            writeln("  p_restrain= ", p_restrain);
            writeln("  x_buffer= ", x_buffer);
            writeln("  with_brakes= ", with_brakes);
            writeln("  ecL_id= ", ecL_id);
            writeln("  ecR_id= ", ecR_id);
        }
    } // end Piston constructor

} // end class Piston
