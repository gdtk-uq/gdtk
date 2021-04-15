// tube.d for the Lagrangian 1D Gas Dynamics, also known as L1d4.
// PA Jacobs
// 2020-04-08 Static tube definition.
// 2021-04-16 Matt McGilvray's heat-transfer-augmentation factor.
//
module tube;

import std.stdio;
import std.string;
import std.format;
import std.conv;
import std.math;
import std.algorithm;

import geom;
import config;

class Tube {
public:
    int n; // The number of small segments that define the tube.
    // The tube properties ar assumed to be linearly distributed
    // between the sample points.
    double[] xs;         // Sample point positions.
    double[] ds;         // Tube diameters.
    double[] areas;      // Cross-section areas.
    // 2021-04-16: So far, we have let the cross-section areas be computed
    // from the diameter.  The areas in the file are just for plotting.
    double[] K_over_Ls;  // Head-loss factor, K, divided by finite distance L.
    double[] Ts;         // Tube-wall temperature.
    double[] vfs;        // Viscous factor, presumably 1.0,
    // but it can be reduced to 0.0 to eliminate viscous wall effects locally.
    double[] htcfs;      // Heat-transfer-coefficient factor.

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
        vfs.length = n+1;
        htcfs.length = n+1;
        foreach (i; 0 .. n+1) {
            text = fp.readln().chomp();
            text.formattedRead!("%e %e %e %e %e %e %e")
                (xs[i], ds[i], areas[i], K_over_Ls[i], Ts[i], vfs[i], htcfs[i]);
        }
        if (L1dConfig.verbosity_level >= 1) {
            int i = 0; writefln("  xs[%d]= %e ds[%d]= %e", i, xs[i], i, ds[i]);
            i = n; writefln("  xs[%d]= %e ds[%d]= %e", i, xs[i], i, ds[i]);
        }
    } // end constructor

    @nogc
    double[6] eval(double x)
    {
        double d, K_over_L, Twall, vf, htcf;
        if (x <= xs[0]) {
            d = ds[0];
            K_over_L = K_over_Ls[0];
            Twall = Ts[0];
            vf = vfs[0];
            htcf = htcfs[0];
        } else if (x >= xs[$-1]) {
            d = ds[$-1];
            K_over_L = K_over_Ls[$-1];
            Twall = Ts[$-1];
            vf = vfs[$-1];
            htcf = htcfs[$-1];
        } else {
            double dx = xs[1] - xs[0];
            int i = cast(int)floor((x - xs[0])/dx);
            double frac = (x - xs[i])/dx;
            i = max(0, min(i, ds.length-2)); // deal with potential round-off error
            d = (1.0-frac)*ds[i] + frac*ds[i+1];
            K_over_L = (1.0-frac)*K_over_Ls[i] + frac*K_over_Ls[i+1];
            Twall = (1.0-frac)*Ts[i] + frac*Ts[i+1];
            vf = (1.0-frac)*vfs[i] + frac*vfs[i+1];
            htcf = (1.0-frac)*htcfs[i] + frac*htcfs[i+1];
        }
        double area = std.math.PI*(d^^2)/4;
        return [d, area, K_over_L, Twall, vf, htcf];
    } // end eval()

} // end class Tube
