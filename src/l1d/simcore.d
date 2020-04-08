// simcore.d for the Lagrangian 1D Gas Dynamics, also known as L1d4.
// PA Jacobs
// 2020-04-08
//
module simcore;

import std.conv;
import std.stdio;
import std.string;
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
import lcell;
import piston;
import endcondition;

__gshared static GasModel[] gmodels;
__gshared static ThermochemicalReactor[] reactors;
__gshared static Tube tube1;
__gshared static Piston[] pistons;
__gshared static GasSlug[] gasslugs;
__gshared static EndCondition[] ecs;
__gshared static Diaphragm[] diaphragms;


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
    tube1 = new Tube(L1dConfig.job_name~"/tube.data");
    //
    foreach (i; 0 .. L1dConfig.nslugs) {
        auto myData = jsonData[format("slug_%d", i)];
        size_t indx = gasslugs.length;
        gasslugs ~= new GasSlug(indx, myData);
    }
    //
    foreach (i; 0 .. L1dConfig.npistons) {
        auto myData = jsonData[format("piston_%d", i)];
        size_t indx = pistons.length;
        pistons ~= new Piston(indx, myData);
    }
    //
    foreach (i; 0 .. L1dConfig.necs) {
        auto myData = jsonData[format("end_condition_%d", i)];
        string ecClass = getJSONstring(myData, "class", "");
        size_t indx = ecs.length;
        switch (ecClass) {
        case "Diaphragm":
            auto myDia = new Diaphragm(indx, myData);
            ecs ~= myDia;
            diaphragms ~= myDia;
            break;
        case "GasInterface":
            ecs ~= new GasInterface(indx, myData);
            break;
        case "FreeEnd":
            ecs ~= new FreeEnd(indx, myData);
            break;
        case "VelocityEnd":
            ecs ~= new VelocityEnd(indx, myData);
            break;
        case "PistonFace":
            ecs ~= new PistonFace(indx, myData);
            break;
        default:
            string msg = text("Unknown EndCondition: ", ecClass);
            throw new Exception(msg);
        }
    }
    //
    // Work through dynamic components and read their initial state.
    foreach (i, s; gasslugs) {
        if (L1dConfig.verbosity_level >= 1) {
            writeln("Initial state of slug ", i);
        }
        string fileName = L1dConfig.job_name ~ format("/slug-%04d-faces.data", i);
        File fp = File(fileName, "r");
        s.read_face_data(fp, 0);
        fp.close();
        fileName = L1dConfig.job_name ~ format("/slug-%04d-cells.data", i);
        fp = File(fileName, "r");
        s.read_cell_data(fp, 0);
        fp.close();
        if (L1dConfig.verbosity_level >= 1) {
            LFace f = s.faces[$-1];
            LCell c = s.cells[$-1];
            writeln(format("  face[%d] x=%e area=%e", s.faces.length-1, f.x, f.area));
            writeln(format("  cell[%d] xmid=%e p=%e T=%e vel=%e",
                           s.cells.length-1, c.xmid, c.gas.p, c.gas.T, c.vel));
        }
    }
    foreach (i, p; pistons) {
        if (L1dConfig.verbosity_level >= 1) {
            writeln("Initial state of piston ", i);
        }
        string fileName = L1dConfig.job_name ~ format("/piston-%04d.data", i);
        File fp = File(fileName, "r");
        p.read_data(fp, 0);
        fp.close();
        if (L1dConfig.verbosity_level >= 1) {
            writeln(format("  x=%e, vel=%e is_restrain=%s brakes_on=%s hit_buffer=%s",
                           p.x, p.vel, p.is_restrain, p.brakes_on, p.hit_buffer));
        }
    }
    foreach (i, dia; diaphragms) {
        if (L1dConfig.verbosity_level >= 1) {
            writeln("Initial state of diaphragm (at EndCondition index)", dia.indx);
        }
        string fileName = L1dConfig.job_name ~ format("/diaphragm-%04d.data", i);
        File fp = File(fileName, "r");
        dia.read_data(fp, 0);
        fp.close();
        if (L1dConfig.verbosity_level >= 1) {
            writeln(format("  is_burst=%s", dia.is_burst));
        }
    }
    return;
} // end init_simulation()
