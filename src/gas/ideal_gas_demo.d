/**
 * idealgas_demo.d
 *
 * Author: Peter J.
 * Version: 2014-06-22
 */

import std.stdio;
import gas.gas_model;
import gas.gas_state;
import gas.ideal_gas;
import util.lua;
import util.lua_service;

void main() {
    writeln("Begin demonstration of using the IdealGas and GasState classes...");
    lua_State* L = init_lua_State();
    doLuaFile(L, "sample-data/ideal-air-gas-model.lua");
    auto gm = new IdealGas(L);
    writeln("species name=", gm.species_name(0));
    writeln("gm=", gm);
    auto gd = GasState(gm, 100.0e3, 300.0);
    writefln("R= %s, pressure= %s, temperature= %s", gm.R(gd), gd.p, gd.T);
    gm.update_thermo_from_pT(gd);
    gm.update_sound_speed(gd);
    writefln("rho= %s, e= %s, a= %s", gd.rho, gd.u, gd.a);
    gm.update_trans_coeffs(gd);
    writefln("mu= %s, k= %s", gd.mu, gd.k);
    writeln("Done.");
}
