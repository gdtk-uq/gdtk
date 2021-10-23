/*
Wrapper functions for interfacing nenzf1d with python, using the ffi library.

Notes:
 - This file is compiled when invoking the target $ make libnenzf1d.so

@author: NNG
*/

import core.runtime;
import core.stdc.string;
import std.string;
import std.stdio;
import std.conv;
import configuration;
import shock_tube_nozzle;

Config config;
Result result;

extern (C) int cwrap_init()
{
    Runtime.initialize();
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
        config.reactions_filename2 = to!string(cfg.reactions_filename2);
    }

    foreach(i; 0 .. cfg.n_species){
        string speciesname = to!string(cfg.species[i]);
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

struct PlainResult{
/*
    Store the calculated gas state at the end of the facility nozzle as a packaged struct of numbers.

    @author: Nick Gibbons
*/
    double x;
    double area_ratio;
    double velocity;
    double Mach_number;
    double p_pitot;
    double rayleigh_pitot;
    double pressure, density, temperature, viscosity;
    int n_T_modes;
    double* T_modes;
    int n_massf;
    double* massf;
    double massflux_rel_err, enthalpy_rel_err, pitot_rel_err;
}

PlainResult d_result_to_c(ref Result result){
/*
    Take the D style struct returned from nenzf1d and convert into a plain C one
    that we can pass back out to python.

    @author: Nick Gibbons
*/
    PlainResult rst;

    rst.x = result.x;
    rst.area_ratio = result.area_ratio;
    rst.velocity = result.velocity;
    rst.Mach_number = result.Mach_number;

    rst.p_pitot = result.p_pitot;
    rst.rayleigh_pitot = result.rayleigh_pitot;
    rst.pressure = result.pressure;
    rst.density = result.density;
    rst.temperature = result.temperature;
    rst.viscosity = result.viscosity;

    rst.n_T_modes = to!int(result.T_modes.length);
    rst.T_modes = result.T_modes.ptr;

    rst.n_massf = to!int(result.massf.length);
    rst.massf = result.massf.ptr;

    rst.massflux_rel_err = result.massflux_rel_err;
    rst.enthalpy_rel_err = result.enthalpy_rel_err;
    rst.pitot_rel_err    = result.pitot_rel_err;

    return rst;
}

extern (C) PlainResult run(int verbosityLevel, PlainConfig cfg)

{
    int exitCode = 0;
    config = plain_to_fancy(cfg);

    PlainResult rst;
    try{
        result = shock_tube_nozzle.analyse(verbosityLevel, config);
        rst = d_result_to_c(result);
    }
    catch(Exception e) {
        writeln("Caught exception in cwrap.run, message was:");
        writeln(e.msg);
        exitCode = 1;
    }

    stdout.flush();
    // Set the thing exitCodePointer is pointing at to be equal to exitCode
    //*exitCodePointer = exitCode;
    return rst;
}

