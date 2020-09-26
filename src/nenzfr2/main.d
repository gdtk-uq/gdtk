// main.d for nenzfr2
// PJ 2020-09-26 Initial code built from Python prototype.

import std.stdio;
import std.array;
import std.string;
import std.conv;
import dyaml;

void main()
{
    writeln("NENZFR2: shock-tunnel with nonequilibrium nozzle flow.");
    writeln("Revision: PUT_REVISION_STRING_HERE");
    writeln("Compiler-name: PUT_COMPILER_NAME_HERE");
    // Read the input.
    auto config = dyaml.Loader.fromFile("input.yaml").load();
    // Extract our parameters and write them to stdout.
    writeln(config["title"].as!string);
    string gm1_filename = config["gas-model-1"].as!string;
    writeln("gas-model-1= ", gm1_filename);
    string gm2_filename = config["gas-model-2"].as!string;
    writeln("gas-model-2= ", gm2_filename);
    string reactions_filename = config["reactions"].as!string;
    writeln("reactions= ", reactions_filename);
    string[] species;
    foreach(string name; config["species"]) { species ~= name; }
    writeln("species= ", species);
    double T1 = to!double(config["T1"].as!string);
    double p1 = to!double(config["p1"].as!string);
    double Vs = to!double(config["Vs"].as!string);
    double pe = to!double(config["pe"].as!string);
    double ar = to!double(config["ar"].as!string);
    writeln("T1= ", T1);
    writeln("p1= ", p1);
    writeln("Vs= ", Vs);
    writeln("pe= ", pe);
    writeln("ar= ", ar);
}
