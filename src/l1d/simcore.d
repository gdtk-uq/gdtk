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
import gasflow;
import config;
import tube;
import gasslug;
import piston;
import endcondition;

__gshared static GasModel[] gmodels;
__gshared static Tube tube1;
__gshared static Piston[] pistons;
__gshared static GasSlug[] gasslugs;
__gshared static EndCondition[] ecs;


void init_simulation(int tindx_start)
{
    string dirName = L1dConfig.job_name;
    string fileName = dirName~"/config.json";
    string content;
    try {
        content = readText(fileName);
    } catch (Exception e) {
        string msg = text("Failed to read config file: ", fileName);
        msg ~= text(" Message is: ", e.msg);
        throw new Exception(msg);
    }
    JSONValue jsonData;
    try {
        jsonData = parseJSON!string(content);
    } catch (Exception e) {
        string msg = text("Failed to parse JSON from config file: ", fileName);
        msg ~= text(" Message is: ", e.msg);
        throw new Exception(msg);
    }
    // Now that we have parsed JSON data, proceed to update those config values.
    
    return;
}
