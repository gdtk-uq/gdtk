/**
 * linesearch_demo.d
 * Example of using the one-dimensional line search for minimization. 
 *
 * Peter J. 2014-06-14
 *
 * Example transcript:
 *
 * peterj@laval ~/cfcfd3/dlang/nm $ ./linesearch_demo 
 * Begin demonstration of line search.
 * x=-1.4721 f(-1.4721)=2.13091
 * bracket=(-0.588534,-0.588533)
 * Done.
 */

import std.stdio;
import std.math;
import linesearch;

void main() {
    writeln("Begin demonstration of line search.");
    double fdemo(double x) {
        return exp(x) + 2.0 - cos(x);
    }
    double x = -1.4721;
    writeln("x=", x, " f(", x, ")=", fdemo(x));

    double a = -3;
    double b = 1;
    minimize!fdemo(a, b, 1.0e-6);
    writeln("bracket=(", a, ",", b, ")");
    writeln("Done.");
}
