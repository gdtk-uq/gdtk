// fv_demo.d
// Exerside the finite-volume classes.
// PJ, 2014-07-17

import std.stdio;
import geom;
import gas;
import flowstate;
import conservedquantities;
import fvcore;
import fvinterface;
import fvvertex;
import fvcell;
import globalconfig;
import fluxcalc;

void main()
{
    writeln("test fv modules");
    auto gm = init_gas_model("sample-data/ideal-air-gas-model.lua");
    auto flow = new FlowState(gm, 100.0e3, [300.0,], Vector3(1.0,0.0,0.0));
    writeln("flow=", flow);
    auto flow2 = new FlowState(gm, 120.0e3, [300.0,], Vector3(0.0,1.0,0.0));
    auto flow3 = new FlowState(gm, 120.0e3, [300.0,], Vector3(0.0,1.0,0.0));
    flow3.copy_average_values_from([flow, flow2], gm);
    writeln("flow3=", flow3);

    writeln("-----------------------");
    auto Q = new ConservedQuantities(gm);
    Q.mass = 99.0;
    Q.energies[0] = 301.0;
    writeln("Q=", Q);
    Q.clear_values();
    writeln("cleared Q=", Q);

    writeln("-----------------------");
    writeln("face_name[face_index(\"south\")]=", face_name[face_index("south")]);
    writeln("update_scheme_name[update_scheme_index[\"classic_rk3\"]]=",
	    gasdynamic_update_scheme_name(update_scheme_from_name("classic_rk3")));

    writeln("-----------------------");
    auto iface = new FVInterface(gm);
    writeln("iface=", iface);

    writeln("-----------------------");
    auto vtx = new FVVertex(gm);
    writeln("vtx=", vtx);

    writeln("-----------------------");
    auto cell = new FVCell(gm);
    GlobalConfig.gmodel = gm; // The following call needs the gas model in place.
    writeln("variable_list_for_cell=", variable_list_for_cell());
    string sample = "1.0 2.0 3.0 0.000999 0.1 1.1 1.2 1.3 100.0e3 345.0 1.8e-5 "
	~ "0.0123 0.999 0.0888 1 0.05 1.009 1.0 2.65e5 311";
    cell.scan_values_from_string(sample);
    cell.encode_conserved(0, 0, 0.0);
    cell.decode_conserved(0, 0, 0.0);
    cell.check_flow_data();
    writeln("cell=", cell);
    writeln("string written=", cell.write_values_to_string());
    writeln("sample=", sample);

    writeln("------------------------");
    adaptive_flux(flow, flow2, iface);
    writeln("iface=", iface);

    writeln("done.");
}
