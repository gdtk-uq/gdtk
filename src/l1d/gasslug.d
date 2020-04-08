// gasslug.d for the Lagrangian 1D Gas Dynamics, also known as L1d4.
// PA Jacobs
// 2020-04-08
//
module gasslug;

import std.conv;
import std.stdio;
import std.string;
import std.json;

import json_helper;
import geom;
import gas;
import gasflow;
import config;
import lcell;
import endcondition;
import simcore; // has the core data arrays

class GasSlug {
public:
    size_t indx;
    string label;
    size_t gmodel_id;
    GasModel gmodel;
    size_t ncells;
    int viscous_effects;
    bool adiabatic;
    int ecL_id;
    int ecR_id;
    EndCondition ecL;
    EndCondition ecR;
    size_t nhcells;
    size_t[] hcells;

    this(size_t indx, JSONValue jsonData)
    {
        if (L1dConfig.verbosity_level >= 3) {
            writeln("construct slug[", indx, "] from json=", jsonData);
        }
        this.indx = indx;
        label = getJSONstring(jsonData, "label", "");
        gmodel_id = getJSONint(jsonData, "gmodel_id", 0);
        gmodel = gmodels[gmodel_id];
        ncells = getJSONint(jsonData, "ncells", 0);
        viscous_effects = getJSONint(jsonData, "viscous_effects", 0);
        adiabatic = getJSONbool(jsonData, "adiabatic", false);
        ecL_id = getJSONint(jsonData, "ecL_id", -1);
        ecR_id = getJSONint(jsonData, "ecR_id", -1);
        nhcells = getJSONint(jsonData, "nhcells", 1);
        hcells = to!(size_t[])(getJSONintarray(jsonData, "hcells", [0,]));
        if (L1dConfig.verbosity_level >= 1) {
            writeln("GasSlug[", indx, "]:");
            writefln("  label= \"%s\"", label);
            writeln("  gmodel_id= ", gmodel_id);
            writeln("  ncells= ", ncells);
            writeln("  viscous_effects= ", viscous_effects);
            writeln("  adiabatic= ", adiabatic);
            writeln("  ecL_id= ", ecL_id);
            writeln("  ecR_id= ", ecR_id);
            writeln("  hcells= ", hcells);
        }        
    }
}
