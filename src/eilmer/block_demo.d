// block_demo.d
// Exercise the Structured-grid blocks as we build them.
// PJ 2014-07-20, 2015-05-05 gmodel now part of Block class.

import std.stdio;
import geom;
import gas;
import globalconfig;
import flowstate;
import sblock;

void main()
{
    writeln("Block demo, structured blocks only for now...");
    GlobalConfig.verbosity_level = 2;
    GlobalConfig.gas_model_file = "sample-data/ideal-air-gas-model.lua";
    auto myConfig = new LocalConfig();

    throw new Error("Bit rot has set in -- don't use this demo.");

    version(old_stuff) {
        auto blk = new SBlock(1, 10, 40, 1, "test");
        blk.assemble_arrays();
        blk.bind_interfaces_and_vertices_to_cells();
        writeln("blk=", blk);
        blk.grid = new StructuredGrid("sample-data/cone20.grid.b0000.t0000.gz", "gziptext");
        blk.sync_vertices_from_underlying_grid(0);
        blk.sync_vertices_to_underlying_grid(0);
        blk.write_underlying_grid("test-grid.txt.gz");
        auto sim_time = blk.read_solution("sample-data/cone20.flow.b0000.t0000.gz");
        blk.write_solution("test-flow.txt.gz", 1.0);

        auto flow = new FlowState(myConfig.gmodel, 100.0e3, [300.0,], Vector3(1.0,0.0,0.0));
        writeln("flow=", flow);

        writeln("Done.");
    }
}
