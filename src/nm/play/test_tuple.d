import std.typecons;
import std.math;
import std.stdio;

Tuple!(double, double) sincos(double x) 
{
    return tuple(cast(double)sin(x), cast(double)cos(x));
}

void main(string[] argv)
{
    double a, b;
    TypeTuple!(a, b) = sincos(1.0);
    writeln("a= ", a, " b= ", b);
    return;
}
