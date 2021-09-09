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

struct PlainConfig
{
    int n_gm1_filename;
    char* gm1_filename;
    int n_gm2_filename;
    char* gm2_filename;
    int n_reactions_filename;
    char* reactions_filename;
    int n_reactions_filename2;
    char* reactions_filename2;
    int nn_species;
    int* n_species;
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
    double t_final;  // fixme, this used to default to 2.0 * x_end / v
    double t_inc;
    double t_inc_factor;
    double t_inc_max;
}

extern (C)
{
    extern PlainConfig cfg;
}

extern (C) void pass_a_struct(PlainConfig cfg)
{
    writeln("Got struct with T1: ", cfg.T1);
    write("Got molef: ");
    foreach(i; 0 .. cfg.n_molef) writef("%f ", cfg.molef[i]);
    writeln("");
    write("Got gm1_filename: ");
    foreach(i; 0 .. cfg.n_gm1_filename) writef("%s", cfg.gm1_filename[i]);
    writeln("");
    writeln("Got species: ");
    foreach(i; 0 .. cfg.nn_species){
        foreach(j; 0 .. cfg.n_species[i]) writef("%s", cfg.species[i][j]);
        writeln("");
    }
}

extern (C) int cwrap_run()
{
    writeln("Hello from d!");
    return 0;
}
