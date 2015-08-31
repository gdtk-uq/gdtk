// newtoncotes_demo.d
// Peter J.
// 2014-06-20, Demo code adapted from the Python version.
//
// Example transcript:
// peterj@helmholtz ~/cfcfd3/dlang/nm $ ./newtoncotes_demo 
// Begin Newton-Cotes integration demonstration...
// Estimates of pi/4: 0.785398 0.785398
// errors: 8.47579e-09 -1.0154e-08
// number of function calls: 315 55
// Done.

import std.stdio;
import std.math;
import newtoncotes;

static uint count1 = 0;
static uint count2 = 0;

double fun1(double x) {
    count1++;
    return abs(x) < 1.0 ? sqrt(1.0 - x*x): 0.0;
}

double fun2(double x) {
    count2++;
    return 1.0 / (1.0 + x * x);
}

void main() {
    writeln("Begin Newton-Cotes integration demonstration...");
    double a = 0.0; 
    double b = 1.0;
    double pi4_1 = integrate!fun1(a, b);
    double pi4_2 = integrate!fun2(a, b);
    writefln("Estimates of pi/4: %s %s", pi4_1, pi4_2);
    writefln("errors: %s %s", PI/4 - pi4_1, PI/4 - pi4_2);
    writefln("number of function calls: %s %s", count1, count2);
    writeln("Done.");
}
