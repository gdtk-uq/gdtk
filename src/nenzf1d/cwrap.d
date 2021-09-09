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

// You may need some code here to convert an Extern Config which is all pointers and stuff to a d config
// with it's own copies of the data packed into nice d arrays. This shouldn't be too difficult
struct Foo 
{   int a;
    double b;
    int n_c;
    double* c;
    int n_name;
    char* name;
}
extern (C)
{
    extern Foo bar;
}

extern (C) void pass_a_struct(Foo bar)
{
    writeln("Got struct with a: ", bar.a, " and b: ", bar.b, " and c: ", *(bar.c));
    writeln("c: ");
    foreach(i; 0 .. bar.n_c) writeln(*(bar.c+i));
    writeln("Name: ");
    foreach(i; 0 .. bar.n_name) writeln(*(bar.name+i));
}

extern (C) int cwrap_run()
{
    writeln("Hello from d!");
    return 0;
}
