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
    auto Q = new ConservedQuantities(gm.n_species, gm.n_modes);
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
    auto iface = new FVInterface(gm, false, false);
    writeln("iface=", iface);

    writeln("-----------------------");
    auto vtx = new FVVertex(gm, false);
    writeln("vtx=", vtx);

    writeln("-----------------------");
    // The following calls need the gas model in place.
    GlobalConfig.gas_model_file = "sample-data/ideal-air-gas-model.lua"; 
    LocalConfig myConfig = new LocalConfig();
    auto cell = new FVCell(myConfig);
    writeln("variable_list_for_cell=", variable_list_for_cell(myConfig.gmodel, myConfig.include_quality,
                                                              myConfig.MHD, myConfig.divergence_cleaning,
                                                              myConfig.radiation));
    string sample = "1.0 2.0 3.0 0.000999 0.1 1.1 1.2 1.3 100.0e3 345.0 1.8e-5 "
        ~ "0.0123 0.999 0.0888 1 0.05 1.009 1.0 2.65e5 311";
    cell.scan_values_from_string(sample, true);
    cell.encode_conserved(0, 0, 0.0);
    cell.decode_conserved(0, 0, 0.0);
    cell.check_flow_data();
    writeln("cell=", cell);
    writeln("string written=", cell.write_values_to_string());
    writeln("sample=", sample);

    writeln("------------------------");
    adaptive_flux(flow, flow2, iface, myConfig.gmodel);
    writeln("iface=", iface);

    writeln("done.");
}
