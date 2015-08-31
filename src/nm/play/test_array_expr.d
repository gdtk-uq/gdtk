// test_array_expr.d

import std.stdio;

double[] f1(double[] x)
{
    double[] result = x.dup;
    result[] = x[] * 3.0;
    return result;
}

void main () 
{
    double[] a = [1.0, 2.0];
    writeln("a= ", a);
    version(broken) {
    double[] b = f1(2.0*a[]);
    }
    version(works) {
    double[] tmp = a.dup;
    tmp[] = 2.0 * a[];
    double[] b = f1(tmp);
    }
    writeln("b= ", b);
}


/+ transcript...
peterj@helmholtz ~/cfcfd3/dlang/nm $ dmd -version=broken test_array_expr.d
Internal error: e2ir.c 2077
peterj@helmholtz ~/cfcfd3/dlang/nm $ dmd -version=works test_array_expr.d 
peterj@helmholtz ~/cfcfd3/dlang/nm $ ./test_array_expr 
a= [1, 2]
b= [6, 12]
peterj@helmholtz ~/cfcfd3/dlang/nm $ dmd -v
DMD32 D Compiler v2.065
Copyright (c) 1999-2013 by Digital Mars written by Walter Bright
Documentation: http://dlang.org/
Usage:
+/
