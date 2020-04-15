// tube.d for the Lagrangian 1D Gas Dynamics, also known as L1d4.
// PA Jacobs
// 2020-04-08
//
module tube;

import std.stdio;
import std.string;
import std.format;
import std.conv;
import std.math;

import geom;
import config;

class Tube {
public:
    int n; // number of small segments
    double[] xs;
    double[] ds;
    double[] areas;
    double[] K_over_Ls;
    double[] Ts;

    this(string fileName)
    {
        auto fp = File(fileName, "r");
        string text = fp.readln().chomp();
        text.formattedRead!"# n= %d"(n);
        if (L1dConfig.verbosity_level >= 1) {
            writefln("Tube:\n  n= %d", n);
        }
        text = fp.readln(); // read and ignore line with variable names
        xs.length = n+1;
        ds.length = n+1;
        areas.length = n+1;
        K_over_Ls.length = n+1;
        Ts.length = n+1;
        foreach (i; 0 .. n+1) {
            text = fp.readln().chomp();
            text.formattedRead!"%e %e %e %e %e"(xs[i], ds[i], areas[i], K_over_Ls[i], Ts[i]);
        }
        if (L1dConfig.verbosity_level >= 1) {
            int i = 0; writefln("  xs[%d]= %e ds[%d]= %e", i, xs[i], i, ds[i]);
            i = n; writefln("  xs[%d]= %e ds[%d]= %e", i, xs[i], i, ds[i]);
        }
    } // end constructor

    @nogc
    double[4] eval(double x)
    {
        double d, K_over_L, Twall;
        if (x <= xs[0]) {
            d = ds[0];
            K_over_L = K_over_Ls[0];
            Twall = Ts[0];
        } else if (x >= xs[$-1]) {
            d = ds[$-1];
            K_over_L = K_over_Ls[$-1];
            Twall = Ts[$-1];
        } else {
            double dx = xs[1] - xs[0];
            int i = cast(int)floor((x - xs[0])/dx);
            double frac = (x - xs[i])/dx;
            d = (1.0-frac)*ds[i] + frac*ds[i+1];
            K_over_L = (1.0-frac)*K_over_Ls[i] + frac*K_over_Ls[i+1];
            Twall = (1.0-frac)*Ts[i] + frac*Ts[i+1];
        }
        double area = std.math.PI*(d^^2)/4;
        double[4] daKT; daKT[0] = d; daKT[1] = area; daKT[2] = K_over_L; daKT[3] = Twall;
        return daKT;
    } // end eval()

} // end class Tube
