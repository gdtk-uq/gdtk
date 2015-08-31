// flowstate_demo.d
// Peter J. 
// 2014-07-17: get serious about building the CFD code infrastructure...

import std.stdio;
import geom;
import gas;
import flowstate;

void main()
{
    writeln("test flowstate module");
    auto gm = init_gas_model("sample-data/ideal-air-gas-model.lua");
    auto gd = new GasState(gm, 100.0e3, 300.0);
    writefln("R= %s, pressure= %s, temperature= %s", gm.R(gd), gd.p, gd.T[0]);
    gm.update_thermo_from_pT(gd);
    gm.update_sound_speed(gd);
    writefln("rho= %s, e= %s, a= %s", gd.rho, gd.e[0], gd.a); 

    auto flow = new FlowState(gm, 100.0e3, [300.0,], Vector3(1.0,0.0,0.0));
    writefln("flow rho= %s, e= %s, a= %s", flow.gas.rho, flow.gas.e[0], flow.gas.a); 
    writeln("flow=", flow);
    writeln("done.");
}

