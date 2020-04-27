// postprocess.d for the Lagrangian 1D Gas Dynamics, also known as L1d4.
// PA Jacobs
// 2020-04-08
//
module postprocess;

import std.conv;
import std.stdio;
import std.string;
import std.json;
import std.file;

import json_helper;
import gas;
import config;
import tube;
import gasslug;
import lcell;
import piston;


void generate_xt_dataset(string varName, int tindxStart, int tindxEnd)
{
    return;
}


void extract_time_slice(int tindx)
{
    return;
}


void extract_piston_history(int pindx)
{
    return;
}
