/*
Wrapper functions for interfacing slf with python, using the ffi library.

Notes:
 - This file is compiled when invoking the target $ make libslf.so

@author: NNG
*/

import core.runtime;
import core.stdc.string;
import std.string;
import std.stdio;
import std.conv;

import slf;

Flame flame;

extern (C) int cwrap_init()
{
    Runtime.initialize();
    return 0;
}


extern (C) int init_slf(const char* file_name)
{
    try {
        flame = new Flame(to!string(file_name));
        return 0;
    } catch (Exception e) {
        stderr.writeln("Failed to construct a new flame.");
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

extern (C) int set_initial_condition()
{
    if (flame !is null) {
        flame.set_initial_condition();
        return 0;
    } else {
        stderr.writeln("Flame not initialised.");
        return -1;
    }
}

extern (C) int run()
{
    int exitCode = 0;

    try{
        exitCode = flame.run();
    }
    catch(Exception e) {
        writeln("Caught exception in cwrap.run, message was:");
        writeln(e.toString());
        exitCode = 1;
    }

    stdout.flush();
    return exitCode;
}

extern (C) int save_solution()
{
    int exitCode = 0;
    try{
        flame.save_solution();
    }
    catch(Exception e) {
        writeln("Caught exception in cwrap.save_solution, message was:");
        writeln(e.toString());
        exitCode = 1;
    }

    stdout.flush();
    return exitCode;
}

extern (C) int save_log()
{
    int exitCode = 0;
    try{
        flame.save_log();
    }
    catch(Exception e) {
        writeln("Caught exception in cwrap.save_log, message was:");
        writeln(e.msg);
        exitCode = 1;
    }

    stdout.flush();
    return exitCode;
}

extern (C) int back_out_scalar_dissipation()
{
    int exitCode = 0;
    try{
        flame.back_out_scalar_dissipation();
    }
    catch(Exception e) {
        writeln("Caught exception in cwrap.back_out_scalar_dissipation, message was:");
        writeln(e.toString());
        exitCode = 1;
    }

    stdout.flush();
    return exitCode;
}

extern (C) int extract_residual(double* R)
{
    int exitCode = 0;
    try{
        second_derivs_from_cent_diffs(flame.pm, flame.U, flame.U2nd);
        compute_residual(flame.gm, flame.reactor, flame.gs, flame.omegaMi, flame.pm,
                         flame.U, flame.U2nd, flame.R);
        foreach(i; 0 .. flame.n) {
            R[i] = flame.R[i].re;
        }
    }
    catch(Exception e) {
        writeln("Caught exception in cwrap.extract_residual, message was:");
        writeln(e.toString());
        exitCode = 1;
    }

    stdout.flush();
    return exitCode;
}

extern (C) int get_nsp() { return to!int(flame.nsp); }
extern (C) int get_neq() { return to!int(flame.neq); }
extern (C) int get_N()   { return to!int(flame.N); }
extern (C) int get_n()   { return to!int(flame.n); }

extern (C) double get_D()    { return to!double(flame.pm.D.re); }
extern (C) double get_p()    { return to!double(flame.pm.p.re); }
extern (C) double get_dZ()   { return to!double(flame.pm.dZ.re); }
extern (C) double get_T0()   { return to!double(flame.pm.T0.re); }
extern (C) double get_T1()   { return to!double(flame.pm.T1.re); }
extern (C) double get_targetGRR()   { return to!double(flame.config.targetGRR); }

extern (C) void set_D(double D) { flame.pm.D.re = D; }
extern (C) void set_p(double p) { flame.pm.p.re = p; }
extern (C) void set_targetGRR(double targetGRR) { flame.config.targetGRR = targetGRR; }

immutable string get_array_function_def = "
extern (C) int get_%s(double* %s) {
    try { foreach(i; 0 .. flame.%s.length) %s[i] = flame.%s[i].re; return 0; }
    catch (Exception e) { stderr.writeln(\"Failed to get %s because reason: \", e.msg); return -1; }
}";

mixin(format(get_array_function_def, "chi","chi","pm.chi","chi","pm.chi","chi"));
mixin(format(get_array_function_def, "Z", "Z",   "pm.Z",  "Z",  "pm.Z",  "Z"));
mixin(format(get_array_function_def, "Y0", "Y0", "pm.Y0", "Y0", "pm.Y0", "Y0"));
mixin(format(get_array_function_def, "Y1", "Y1", "pm.Y1", "Y1", "pm.Y1", "Y1"));
mixin(format(get_array_function_def, "U0", "U0", "pm.U0", "U0", "pm.U0", "U0"));
mixin(format(get_array_function_def, "U1", "U1", "pm.U1", "U1", "pm.U1", "U1"));
mixin(format(get_array_function_def, "U", "U", "U", "U", "U", "U"));
mixin(format(get_array_function_def, "R", "R", "R", "R", "R", "R"));

extern (C) int get_T(double* T) {
    size_t nsp = flame.nsp;
    size_t neq = flame.neq;
    try {
        foreach(i; 0 .. flame.N) T[i] = flame.U[i*neq + nsp].re;
        return 0;
    }
    catch (Exception e) {
        stderr.writeln("Failed to get temperature because reason: ", e.msg);
        return -1;
    }
}

extern (C) int get_Y(double* Y) {
    size_t nsp = flame.nsp;
    size_t neq = flame.neq;
    size_t idx = 0;
    try {
        foreach(i; 0 .. flame.N) {
            foreach(j; 0 .. flame.nsp) {
                Y[idx] = flame.U[i*neq + j].re;
                idx += 1;
            }
        }
        return 0;
    }
    catch (Exception e) {
        stderr.writeln("Failed to get temperature because reason: ", e.msg);
        return -1;
    }
}

extern (C) int set_U(double* U) {
    try {
        foreach(i; 0 .. flame.n) {
            flame.U[i].re = U[i];
            flame.U[i].im = 0.0;
        }
        return 0;
    }
    catch (Exception e) {
        stderr.writeln("Failed to set conserved quantities U because reason: ", e.msg);
        return -1;
    }
}

extern (C) int set_chi(double* chi) {
    try {
        foreach(i; 0 .. flame.N) flame.pm.chi[i].re = chi[i];
        foreach(i; 0 .. flame.N) flame.pm.chi[i].im = 0.0;
        return 0;
    }
    catch (Exception e) {
        stderr.writeln("Failed to set scalar dissipation because reason: ", e.msg);
        return -1;
    }
}
