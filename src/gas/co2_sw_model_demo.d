/**
 * CO2SW_model_demo.d
 *
 * Author: Jonathan H.
 * Version: 2015-07-22
 */

import std.stdio;
import gas.gas_model;
import std.datetime;

void main() {
    writeln("Begin demonstration of using the gasmodel and Gas_data classes using CO2 Span Wagner...");
    auto gm = init_gas_model("sample-data/co2sw-gas-model.lua");
    foreach(i; 0 .. gm.n_species) {
        writeln("species[", i, "] name=", gm.species_name(i));
    }
    auto gd = new GasState(gm, 6.0e6, 300.0);
    writefln("R= %s, pressure= %.8f, temperature= %.8f", gm.R(gd), gd.p, gd.T);
    gm.update_thermo_from_pT(gd);
    gm.update_sound_speed(gd);
    double s = gm.entropy(gd);
    writefln("rho= %.8f, e= %.8f, a= %.8f, s=%.8f", gd.rho, gd.u, gd.a, s); 
    writeln("------------------------------------------");
    writeln("Update again from p, s");
    writeln("------------------------------------------");
    int ncycles = 100;
    StopWatch sw;
    sw.start();
    for(int i = 0; i != ncycles; i++) gm.update_thermo_from_ps(gd,s);
    sw.stop();
    gm.update_sound_speed(gd);
    long t_eval_ps = sw.peek().usecs;
    writefln("R= %.8f, pressure= %.16f, temperature= %.8f", gm.R(gd), gd.p, gd.T);
    writefln("rho= %.8f, e= %.8f, a= %.8f, s=%.8f", gd.rho, gd.u, gd.a, s); 
    writeln("------------------------------------------");
    writeln("Using LUT to update again from rho, e");
    writeln("------------------------------------------");
    sw.start();
    for(int i = 0; i != ncycles; i++) gm.update_thermo_from_rhou(gd);
    sw.stop();
    long t_lut = sw.peek().usecs - t_eval_ps;
    gm.update_sound_speed(gd);
    writefln("R= %.8f, pressure= %.16f, temperature= %.8f", gm.R(gd), gd.p, gd.T);
    writefln("rho= %.8f, e= %.8f, a= %.8f", gd.rho, gd.u, gd.a); 
    writeln("gd= ", gd);
    writeln("-------------------------------");
    writefln("Time taken for %s evaluations:", ncycles);
    writefln("Update_ps: %s usecs; LUT_rhoe: %s usecs;", t_eval_ps, t_lut);
    writefln("Update_ps is %s times slower", t_eval_ps/t_lut);
    writeln("-------------------------------");
    auto gd2 = new GasState(gm, 200.0e3, 400.0);
    writeln("gd2=", gd2);
    auto gd3 = new GasState(gm, 100.0e3, 300.0);
    gd3.copy_average_values_from([gd2, gd], gm);
    writeln("after average gd3=", gd3);
    gd2.copy_values_from(gd);
    writeln("after copy gd2=", gd2);
    writeln("Done.");
}
