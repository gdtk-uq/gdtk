// simcore.d for the Lagrangian 1D Gas Dynamics, also known as L1d4.
// PA Jacobs
// 2020-04-08
//
module simcore;

import std.conv;
import std.stdio;
import std.string;
//import std.typecons;
import std.json;
import std.file;

import json_helper;
import geom;
import gas;
import kinetics;
import gasflow;
import config;
import tube;
import gasslug;
import piston;
import endcondition;

__gshared static GasModel[] gmodels;
__gshared static ThermochemicalReactor[] reactors;
__gshared static Tube tube1;
__gshared static Piston[] pistons;
__gshared static GasSlug[] gasslugs;
__gshared static EndCondition[] ecs;


void init_simulation(int tindx_start)
{
    string dirName = L1dConfig.job_name;
    string configFileName = dirName~"/config.json";
    string content;
    try {
        content = readText(configFileName);
    } catch (Exception e) {
        string msg = text("Failed to read config file: ", configFileName);
        msg ~= text(" Message is: ", e.msg);
        throw new Exception(msg);
    }
    JSONValue jsonData;
    try {
        jsonData = parseJSON!string(content);
    } catch (Exception e) {
        string msg = text("Failed to parse JSON from config file: ", configFileName);
        msg ~= text(" Message is: ", e.msg);
        throw new Exception(msg);
    }
    // Now that we have parsed JSON data, proceed to update those config values.
    auto configData = jsonData["config"];
    L1dConfig.title = getJSONstring(configData, "title", "");
    L1dConfig.gas_model_files = getJSONstringarray(configData, "gas_model_files", []);
    L1dConfig.reaction_files_1 = getJSONstringarray(configData, "reaction_files_1", []);
    L1dConfig.reaction_files_2 = getJSONstringarray(configData, "reaction_files_2", []);
    L1dConfig.reacting = getJSONbool(configData, "reacting", false);
    if (L1dConfig.verbosity_level >= 1) {
        writeln("Config:");
        writefln("  title= \"%s\"", L1dConfig.title);
        writeln("  gas_model_files= ", L1dConfig.gas_model_files);
        writeln("  reaction_files_1= ", L1dConfig.reaction_files_1);
        writeln("  reaction_files_2= ", L1dConfig.reaction_files_2);
        writeln("  reacting= ", L1dConfig.reacting);
    }
    assert(L1dConfig.gas_model_files.length == L1dConfig.reaction_files_1.length &&
           L1dConfig.gas_model_files.length == L1dConfig.reaction_files_2.length,
           "Lengths of gas model and reaction file lists are inconsistent.");
    foreach (i, fileName; L1dConfig.gas_model_files) {
        auto gm = init_gas_model(fileName);
        gmodels ~= gm;
        auto fn1 = L1dConfig.reaction_files_1[i];
        auto fn2 = L1dConfig.reaction_files_2[i];
        if (fn1.length > 0) {
            reactors ~= init_thermochemical_reactor(gm, fn1, fn2);
        } else {
            reactors ~= null;
        }
    }
    L1dConfig.max_time = getJSONdouble(configData, "max_time", 0.0);
    L1dConfig.max_step = getJSONint(configData, "max_step", 0);
    L1dConfig.dt_init = getJSONdouble(configData, "dt_init", 0.0);
    L1dConfig.cfl_value = getJSONdouble(configData, "cfl_value", 0.0);
    L1dConfig.x_order = getJSONint(configData, "x_order", 0);
    L1dConfig.t_order = getJSONint(configData, "t_order", 0);
    L1dConfig.nslugs = getJSONint(configData, "nslugs", 0);
    L1dConfig.npistons = getJSONint(configData, "npistons", 0);
    L1dConfig.ndiaphragms = getJSONint(configData, "ndiaphragms", 0);
    L1dConfig.necs = getJSONint(configData, "necs", 0);
    if (L1dConfig.verbosity_level >= 1) {
        writeln("  max_time= ", L1dConfig.max_time);
        writeln("  max_step= ", L1dConfig.max_step);
        writeln("  dt_init= ", L1dConfig.dt_init);
        writeln("  cfl_value= ", L1dConfig.cfl_value);
        writeln("  x_order= ", L1dConfig.x_order);
        writeln("  t_order= ", L1dConfig.t_order);
        writeln("  nslugs= ", L1dConfig.nslugs);
        writeln("  npistons= ", L1dConfig.npistons);
        writeln("  ndiaphragms= ", L1dConfig.ndiaphragms);
        writeln("  necs= ", L1dConfig.necs);
    }
    //
    tube1 = new Tube();
    // [TODO] tube datails
    //
    foreach (i; 0 .. L1dConfig.nslugs) {
        auto myData = jsonData[format("slug_%d", i)];
        writeln("slug[", i, "] json=", myData);
        gasslugs ~= new GasSlug();
        // [TODO] config
        // [TODO] initial state
    }
    //
    foreach (i; 0 .. L1dConfig.npistons) {
        auto myData = jsonData[format("piston_%d", i)];
        writeln("piston[", i, "] json=", myData);
        pistons ~= new Piston();
        // [TODO] config
        // [TODO] initial state
    }
    //
    foreach (i; 0 .. L1dConfig.necs) {
        auto myData = jsonData[format("end_condition_%d", i)];
        writeln("end-condition[", i, "] json=", myData);
        // [TODO] need specific classes
        ecs ~= new EndCondition();
        // [TODO] if diaphragm, read initial state
    }
    return;
}
