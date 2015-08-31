/**
 * idealgas_demo.d
 *
 * Author: Peter J.
 * Version: 2014-06-22
 */

import std.stdio;
import gas.gas_model;
import gas.ideal_gas;

void main() {
    writeln("Begin demonstration of using the IdealGas and GasState classes...");
    auto gm = new IdealGas();
    writeln("species name=", gm.species_name(0));
    writeln("gm=", gm);
    auto gd = new GasState(gm, 100.0e3, 300.0);
    writefln("R= %s, pressure= %s, temperature= %s", gm.R(gd), gd.p, gd.T[0]);
    gm.update_thermo_from_pT(gd);
    gm.update_sound_speed(gd);
    writefln("rho= %s, e= %s, a= %s", gd.rho, gd.e[0], gd.a);
    gm.update_trans_coeffs(gd);
    writefln("mu= %s, k= %s", gd.mu, gd.k[0]);
    writeln("Done.");
}
