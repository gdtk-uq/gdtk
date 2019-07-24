// cwrap_gasmodule.d
//
// This particular module will be compiled into a loadable library and
// provides a C-level API for the GasModel and GasState classes
// that may be called from any language with a C-foreign-function-interface.
//
// PJ 2019-07-24: just enough to try integration with the gas makefile.
//

import std.stdio;
import std.format;
import std.conv;
import core.runtime;

extern (C) int cwrap_gas_module_init()
{
    writeln("cwrap_gas_module_init() start of call");
    Runtime.initialize();
    double[] b;
    foreach(i; 0 .. 5) { b ~= to!double(i); }
    writeln("b= ", b);
    return 0;
}

shared static this()
{
    writeln("libgasmodule.so shared static this");
}

shared static ~this()
{
    writeln("libgasmodule.so shared static ~this");
}
