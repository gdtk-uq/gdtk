// cwrap.d
//
// Module for interfacing with foreign function libraries, specifically python's ffi
//
// @author: NNG
//

import core.runtime;
import core.stdc.string;
import std.string;
import std.stdio;
import std.conv;
import nenzf1d;

extern (C) int cwrap_init()
{
    writeln("Hello from d!");
    Runtime.initialize();
    writeln("Initialized!");
    return 0;
}

struct PlainConfig
{
    char* gm1_filename;
    char* gm2_filename;
    char* reactions_filename;
    char* reactions_filename2;
    int n_species;
    char** species;
    int n_molef;
    double* molef;
    double T1;
    double p1;
    double Vs;
    double pe;
    double meq_throat;
    double ar;
    double pp_ps;
    double C;
    int n_xi;
    double* xi;
    int n_di;
    double* di;
    double x_end;
    double t_final;
    double t_inc;
    double t_inc_factor;
    double t_inc_max;
}

Config plain_to_fancy(PlainConfig cfg){
/*
    Take a C style struct that is readable from python and convert it to a D-style
    config struct, of the kind that nenzf1d is expecting. 

    @author: Nick Gibbons
*/
    Config config;
    config.gm1_filename = to!string(cfg.gm1_filename);
    config.gm2_filename = to!string(cfg.gm2_filename);
    config.reactions_filename = to!string(cfg.reactions_filename);

    // There should be more defensive programming here.
    if (cfg.reactions_filename2 !is null){
        writeln("Found default for reactions_filename2");
        config.reactions_filename2 = to!string(cfg.reactions_filename2);
    }

    foreach(i; 0 .. cfg.n_species){
        string speciesname = to!string(cfg.species[i]);
        writeln("Got speciesname: ", speciesname);
        config.species ~= speciesname;
    }

    foreach(i; 0 .. cfg.n_molef) config.molef ~= cfg.molef[i];

    config.T1         = cfg.T1;
    config.p1         = cfg.p1;
    config.Vs         = cfg.Vs;
    if (cfg.pe!=0.0)         config.pe         = cfg.pe;
    if (cfg.meq_throat!=0.0) config.meq_throat = cfg.meq_throat;
    if (cfg.ar!=0.0)         config.ar         = cfg.ar;
    if (cfg.pp_ps!=0.0)      config.pp_ps      = cfg.pp_ps;
    if (cfg.C!=0.0)          config.C          = cfg.C;

    foreach(i; 0 .. cfg.n_xi) config.xi ~= cfg.xi[i];
    foreach(i; 0 .. cfg.n_di) config.di ~= cfg.di[i];

    if (cfg.x_end==0.0)        config.x_end        = config.xi[$-1];
    if (cfg.t_final!=0.0)      config.t_final      = cfg.t_final;
    if (cfg.t_inc!=0.0)        config.t_inc        = cfg.t_inc;
    if (cfg.t_inc_factor!=0.0) config.t_inc_factor = cfg.t_inc_factor;
    if (cfg.t_inc_max!=0.0)    config.t_inc_max    = cfg.t_inc_max;

    return config;
}

extern (C)
{
    extern PlainConfig cfg;
}

extern (C) int run(int verbosityLevel, PlainConfig cfg)

{
    Config config = plain_to_fancy(cfg);
    //writeln("Produced fancy config parameters:");
    //writeln("  gas-model-1= ", config.gm1_filename);
    //writeln("  gas-model-2= ", config.gm2_filename);
    //writeln("  reactions-file= ", config.reactions_filename);
    //writeln("  reactions_file2= ", config.reactions_filename2);
    //writeln("  species= ", config.species);
    //writeln("  molef= ", config.molef);
    //writeln("  T1= ", config.T1);
    //writeln("  p1= ", config.p1);
    //writeln("  Vs= ", config.Vs);
    //writeln("  pe= ", config.pe);
    //writeln("  meq_throat= ", config.meq_throat);
    //writeln("  ar= ", config.ar);
    //writeln("  pp_ps= ", config.pp_ps);
    //writeln("  C= ", config.C);
    //writeln("  xi= ", config.xi);
    //writeln("  di= ", config.di);
    //writeln("  x_end= ", config.x_end);
    //writeln("  t_final= ", config.t_final);
    int exitCode=0;
    try{
        exitCode = nenzf1d.run(verbosityLevel, config);
    }
    catch(Exception e) {
        writeln("Caught exception in cwrap.run, message was:");
        writeln(e.msg);
        exitCode = 1;
    }
    return exitCode;
}

