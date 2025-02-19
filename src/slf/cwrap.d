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
        writeln(e.msg);
        exitCode = 1;
    }

    stdout.flush();
    return exitCode;
}

extern (C) int get_size_of_solution()
{
    return to!int(flame.n);
}


extern (C) int get_solution(double* U)
{
    // Return the flame field as a numpy array
    try {
        foreach(i; 0 .. flame.U.length) U[i] = flame.U[i].re;
        return 0;
    } catch (Exception e) {
        stderr.writeln("Exception message: ", e.msg);
        return -1;
    }
}

