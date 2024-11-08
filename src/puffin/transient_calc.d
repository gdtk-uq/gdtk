// transient_calc.d -- Part of the Lorikeet transient-flow calculator.
//
// PA Jacobs
// 2022-12-12: Adapt from Puffin and Chicken codes.
//
module transient_calc;

import std.conv;
import std.stdio;
import std.string;
import std.json;
import std.file;
import std.datetime;
import std.format;
import std.range;
import std.math;
import std.algorithm;
import core.stdc.math: HUGE_VAL;
import std.parallelism;

import util.json_helper;
import geom;
import gas;
import kinetics;
import config;
import flow;
import fluidblock;

// We use __gshared so that several threads may access
// the following array concurrently.
__gshared static FluidBlock[] fluidBlocks;

struct ProgressData {
    int step = 0; // steps so far
    double t = 0.0; // time at start of gas-dynamic update
    double dt = 0.0; // increment of time for the gas-dynamic update
    int tindx = 0; // index of the flow data just read or written
    double plot_at_t = 0.0; // time at which to write another blob of flow data
    double[] dt_values; // a place to store the allowable dt for each block
    int steps_since_last_plot_write = 0;
    SysTime wall_clock_start;
}

__gshared static ProgressData progress;

void init_simulation(int tindx)
// Set up configuration and read block data for a given tindx.
{
    string dirName = Config.job_name;
    JSONValue configData = readJSONfile(dirName~"/config.json");
    parse_config_data_for_transient_solver(configData);
    //
    // Parallel runtime will have main thread plus extras.
    int extraThreadsInPool = min(Config.maxCPUs-1, Config.n_fluid_blocks-1);
    defaultPoolThreads(extraThreadsInPool);
    //
    JSONValue blockData = configData["fluid_blocks"];
    foreach (i; 0 .. Config.n_fluid_blocks) {
        fluidBlocks ~= new FluidBlock(i, blockData[i]);
        if (Config.verbosity_level > 1) {
            writefln("  fluidBlocks[%d]= %s", i, fluidBlocks[i]);
        }
    }
    //
    // We set up all blocks fully, whether they are active or not.
    // This is because we need a full complement of block partitions
    // for the VTK structured grid that is built for Paraview.
    foreach (b; fluidBlocks) {
        b.set_up_data_storage();
        b.read_grid_data();
        b.set_up_geometry();
        b.read_flow_data(tindx);
        // [TODO] Should shift the initialization of the kinetics into fluidblock module.
        if (Config.reacting) {
            b.thermochemUpdate =
                init_thermochemical_reactor(b.gmodel, Config.reaction_file_1, Config.reaction_file_2);
        }
    }
    // Read times file to get our starting time for this simulation.
    // This allows for restarted simulations, with tindx > 0.
    double start_t = 0.0;
    auto tf = File(format("%s/times.data", Config.job_name), "r");
    bool tindxFound = false;
    foreach (string line; lines(tf)) {
        if (line.startsWith("#")) continue;
        auto items = line.split();
        int index = to!int(items[0]);
        if (index == tindx) {
            tindxFound = true;
            start_t = to!double(items[1]);
            break;
        }
    }
    if (!tindxFound) {
        throw new Exception("Did not find requested starting tindx in times file.");
    }
    writefln("Start simulation at t= %g", start_t);
    progress.step = 0;
    progress.t = start_t;
    progress.dt = Config.dt_init;
    progress.tindx = tindx;
    progress.plot_at_t = Config.plot_dt;
    progress.steps_since_last_plot_write = 0;
    progress.dt_values.length = fluidBlocks.length;
    return;
} // end init_calculation()


void do_time_integration()
{
    progress.wall_clock_start = Clock.currTime();
    foreach (b; fluidBlocks) {
        if (b.active) { b.encode_conserved(0); }
    }
    while (progress.t < Config.max_t && progress.step < Config.max_step) {
        //
        // 1. Occasionally set size of time step.
        if (progress.step > 0 && (progress.step % Config.cfl_count)==0) {
            foreach (j, b; parallel(fluidBlocks, 1)) {
                progress.dt_values[j] = (b.active) ? b.estimate_allowable_dt() : HUGE_VAL;
            }
            double smallest_dt = progress.dt_values.minElement();
            // Make the transition to larger allowable time step not so sudden.
            progress.dt = fmin(1.5*progress.dt, smallest_dt);
        }
        //
        // 2. Take a step.
        int attempt_number = 0;
        bool step_failed;
        do {
            ++attempt_number;
            step_failed = false;
            try {
                // 1. Predictor (Euler) stage.
                apply_boundary_conditions();
                foreach (b; parallel(fluidBlocks, 1)) {
                    if (!b.active) continue;
                    b.mark_shock_cells();
                    b.update_conserved_for_stage(1, progress.dt);
                }
                // 2. Corrector stage.
                if (Config.t_order > 1) {
                    apply_boundary_conditions();
                    foreach (b; parallel(fluidBlocks, 1)) {
                        if (!b.active) continue;
                        b.update_conserved_for_stage(2, progress.dt);
                    }
                }
            } catch (Exception e) {
                writefln("Step failed e.msg=%s", e.msg);
                step_failed = true;
                progress.dt *= 0.2;
                // We need to restore the flow-field.
                foreach (b; parallel(fluidBlocks, 1)) {
                    if (b.active) { b.decode_conserved(0); }
                }
            }
        } while (step_failed && (attempt_number <= 3));
        if (step_failed) {
            throw new Exception("Gasdynamic step failed after 3 attempts.");
        }
        //
        // 3. Prepare for next time step.
        foreach (b; parallel(fluidBlocks, 1)) {
            if (b.active) { b.transfer_conserved_quantities(1, 0); }
        }
        //
        // 4. Chemistry step is loosely coupled to the gas dynamics.
        if (Config.reacting) {
            foreach (b; parallel(fluidBlocks, 1)) {
                if (b.active) { b.thermochemical_increment(progress.dt); }
            }
        }
        //
        progress.t += progress.dt;
        progress.step++;
        //
        // 5. Occasional console output.
        if (Config.verbosity_level >= 1 &&
            ((progress.step % Config.print_count) == 0)) {
            // For reporting wall-clock time, convert with precision of milliseconds.
            auto elapsed_ms = (Clock.currTime() - progress.wall_clock_start).total!"msecs"();
            double elapsed_s = to!double(elapsed_ms)/1000;
            double WCtFT = ((progress.t > 0.0) && (progress.step > 0)) ?
                elapsed_s*(Config.max_t-progress.t)/progress.dt/progress.step : 0.0;
            writefln("Step=%d t=%.3e dt=%.3e WC=%.3f WCtFT=%.3f",
                     progress.step, progress.t, progress.dt, elapsed_s, WCtFT);
            stdout.flush();
        }
        //
        // 5. Write a flow solution (maybe).
        if (progress.t >= progress.plot_at_t) {
            int tindx = progress.tindx + 1;
            foreach (b; fluidBlocks) { b.write_flow_data(tindx); }
            append(format("%s/times.data", Config.job_name), format("%d %g\n", tindx, progress.t));
            progress.steps_since_last_plot_write = 0;
            progress.plot_at_t += Config.plot_dt;
            progress.tindx = tindx;
        } else {
            progress.steps_since_last_plot_write++;
        }
    } // end while
    //
    // Write the final slice, maybe.
    if (progress.steps_since_last_plot_write > 0) {
        int tindx = progress.tindx + 1;
        foreach (b; fluidBlocks) { b.write_flow_data(tindx); }
        append(format("%s/times.data", Config.job_name), format("%d %g\n", tindx, progress.t));
        progress.tindx = tindx;
    }
    return;
} // end do_time integration()


void apply_boundary_conditions()
// Application of the boundary conditions is essentially filling the
// ghost-cell flow states with suitable data.
// We do this in a context that has a view of all blocks so that exchange
// of flow data between blocks is easy.
{
    foreach (b; parallel(fluidBlocks, 1)) {
        if (!b.active) continue;
        //
        final switch (b.bc_west.code) {
        case BCCode.wall_with_slip: {
            // The slip-wall condition is implemented by filling the ghost cells
            // with reflected-normal-velocity flow.
            foreach (j; 0 .. b.njc) {
                auto face = b.ifaces[b.iface_index(0,j)];
                FlowState2D* fstate = &(face.left_cells[0].fs);
                fstate.copy_values_from(face.right_cells[0].fs);
                fstate.vel.transform_to_local_frame(face.n, face.t1);
                fstate.vel.x = -(fstate.vel.x);
                fstate.vel.transform_to_global_frame(face.n, face.t1);
                //
                fstate = &(face.left_cells[1].fs);
                fstate.copy_values_from(face.right_cells[1].fs);
                fstate.vel.transform_to_local_frame(face.n, face.t1);
                fstate.vel.x = -(fstate.vel.x);
                fstate.vel.transform_to_global_frame(face.n, face.t1);
            }
        } break;
        case BCCode.exchange: {
            // We fill the ghost-cells with data from the neighbour fluid block.
            int other_i = b.i - 1;
            if (other_i < 0) { other_i = Config.nib-1; } // Wrap around.
            int other_j = b.j;
            int other_id = Config.blk_ids[other_i][other_j];
            auto other_b = fluidBlocks[other_id];
            foreach (j; 0 .. b.njc) {
                auto face = b.ifaces[b.iface_index(0,j)];
                auto other_face = other_b.ifaces[other_b.iface_index(other_b.nic,j)];
                face.left_cells[0].fs.copy_values_from(other_face.left_cells[0].fs);
                face.left_cells[1].fs.copy_values_from(other_face.left_cells[1].fs);
            }
        } break;
        case BCCode.inflow: {
            foreach (j; 0 .. b.njc) {
                auto face = b.ifaces[b.iface_index(0,j)];
                face.left_cells[0].fs.copy_values_from(*(b.bc_west.fs));
                face.left_cells[1].fs.copy_values_from(*(b.bc_west.fs));
            }
        } break;
        case BCCode.outflow: {
            foreach (j; 0 .. b.njc) {
                auto face = b.ifaces[b.iface_index(0,j)];
                face.left_cells[0].fs.copy_values_from(face.right_cells[0].fs);
                face.left_cells[1].fs.copy_values_from(face.right_cells[0].fs);
            }
        }
        } // end switch bc_west
        //
        final switch (b.bc_east.code) {
        case BCCode.wall_with_slip: {
            // The slip-wall condition is implemented by filling the ghost cells
            // with reflected-normal-velocity flow.
            foreach (j; 0 .. b.njc) {
                auto face = b.ifaces[b.iface_index(b.nic,j)];
                FlowState2D* fstate = &(face.right_cells[0].fs);
                fstate.copy_values_from(face.left_cells[0].fs);
                fstate.vel.transform_to_local_frame(face.n, face.t1);
                fstate.vel.x = -(fstate.vel.x);
                fstate.vel.transform_to_global_frame(face.n, face.t1);
                //
                fstate = &(face.right_cells[1].fs);
                fstate.copy_values_from(face.left_cells[1].fs);
                fstate.vel.transform_to_local_frame(face.n, face.t1);
                fstate.vel.x = -(fstate.vel.x);
                fstate.vel.transform_to_global_frame(face.n, face.t1);
            }
        } break;
        case BCCode.exchange: {
            // We fill the ghost-cells with data from the neighbour fluid block.
            int other_i = b.i + 1;
            if (other_i >= Config.nib) { other_i = 0; } // Wrap around.
            int other_j = b.j;
            int other_id = Config.blk_ids[other_i][other_j];
            auto other_b = fluidBlocks[other_id];
            foreach (j; 0 .. b.njc) {
                auto face = b.ifaces[b.iface_index(b.nic,j)];
                auto other_face = other_b.ifaces[other_b.iface_index(0,j)];
                face.right_cells[0].fs.copy_values_from(other_face.right_cells[0].fs);
                face.right_cells[1].fs.copy_values_from(other_face.right_cells[1].fs);
            }
        } break;
        case BCCode.inflow: {
            foreach (j; 0 .. b.njc) {
                auto face = b.ifaces[b.iface_index(b.nic,j)];
                face.right_cells[0].fs.copy_values_from(*(b.bc_east.fs));
                face.right_cells[1].fs.copy_values_from(*(b.bc_east.fs));
            }
        } break;
        case BCCode.outflow: {
            foreach (j; 0 .. b.njc) {
                auto face = b.ifaces[b.iface_index(b.nic,j)];
                face.right_cells[0].fs.copy_values_from(face.left_cells[0].fs);
                face.right_cells[1].fs.copy_values_from(face.left_cells[0].fs);
            }
        }
        } // end switch bc_east
        //
        final switch (b.bc_south.code) {
        case BCCode.wall_with_slip: {
            // The slip-wall condition is implemented by filling the ghost cells
            // with reflected-normal-velocity flow.
            foreach (i; 0 .. b.nic) {
                auto face = b.jfaces[b.jface_index(i,0)];
                FlowState2D* fstate = &(face.left_cells[0].fs);
                fstate.copy_values_from(face.right_cells[0].fs);
                fstate.vel.transform_to_local_frame(face.n, face.t1);
                fstate.vel.x = -(fstate.vel.x);
                fstate.vel.transform_to_global_frame(face.n, face.t1);
                //
                fstate = &(face.left_cells[1].fs);
                fstate.copy_values_from(face.right_cells[1].fs);
                fstate.vel.transform_to_local_frame(face.n, face.t1);
                fstate.vel.x = -(fstate.vel.x);
                fstate.vel.transform_to_global_frame(face.n, face.t1);
            }
        } break;
        case BCCode.exchange: {
            // We fill the ghost-cells with data from the neighbour fluid block.
            int other_i = b.i;
            int other_j = b.j - 1;
            if (other_j < 0) { other_j = Config.njb-1; } // Wrap around.
            int other_id = Config.blk_ids[other_i][other_j];
            auto other_b = fluidBlocks[other_id];
            foreach (i; 0 .. b.nic) {
                auto face = b.jfaces[b.jface_index(i,0)];
                auto other_face = other_b.jfaces[other_b.jface_index(i,other_b.njc)];
                face.left_cells[0].fs.copy_values_from(other_face.left_cells[0].fs);
                face.left_cells[1].fs.copy_values_from(other_face.left_cells[1].fs);
            }
        } break;
        case BCCode.inflow: {
            foreach (i; 0 .. b.nic) {
                auto face = b.jfaces[b.jface_index(i,0)];
                face.left_cells[0].fs.copy_values_from(*(b.bc_south.fs));
                face.left_cells[1].fs.copy_values_from(*(b.bc_south.fs));
            }
        } break;
        case BCCode.outflow: {
            foreach (i; 0 .. b.nic) {
                auto face = b.jfaces[b.jface_index(i,0)];
                face.left_cells[0].fs.copy_values_from(face.right_cells[0].fs);
                face.left_cells[1].fs.copy_values_from(face.right_cells[0].fs);
            }
        }
        } // end switch bc_south
        //
        final switch (b.bc_north.code) {
        case BCCode.wall_with_slip: {
            // The slip-wall condition is implemented by filling the ghost cells
            // with reflected-normal-velocity flow.
            foreach (i; 0 .. b.nic) {
                auto face = b.jfaces[b.jface_index(i,b.njc)];
                FlowState2D* fstate = &(face.right_cells[0].fs);
                fstate.copy_values_from(face.left_cells[0].fs);
                fstate.vel.transform_to_local_frame(face.n, face.t1);
                fstate.vel.x = -(fstate.vel.x);
                fstate.vel.transform_to_global_frame(face.n, face.t1);
                //
                fstate = &(face.right_cells[1].fs);
                fstate.copy_values_from(face.left_cells[1].fs);
                fstate.vel.transform_to_local_frame(face.n, face.t1);
                fstate.vel.x = -(fstate.vel.x);
                fstate.vel.transform_to_global_frame(face.n, face.t1);
            }
        } break;
        case BCCode.exchange: {
            // We fill the ghost-cells with data from the neighbour fluid block.
            int other_i = b.i;
            int other_j = b.j + 1;
            if (other_j >= Config.njb) { other_j = 0; } // Wrap around.
            int other_id = Config.blk_ids[other_i][other_j];
            auto other_b = fluidBlocks[other_id];
            foreach (i; 0 .. b.nic) {
                auto face = b.jfaces[b.jface_index(i,b.njc)];
                auto other_face = other_b.jfaces[other_b.jface_index(i,0)];
                face.right_cells[0].fs.copy_values_from(other_face.right_cells[0].fs);
                face.right_cells[1].fs.copy_values_from(other_face.right_cells[1].fs);
            }
        } break;
        case BCCode.inflow: {
            foreach (i; 0 .. b.nic) {
                auto face = b.jfaces[b.jface_index(i,b.njc)];
                face.right_cells[0].fs.copy_values_from(*(b.bc_north.fs));
                face.right_cells[1].fs.copy_values_from(*(b.bc_north.fs));
            }
        } break;
        case BCCode.outflow: {
            foreach (i; 0 .. b.nic) {
                auto face = b.jfaces[b.jface_index(i,b.njc)];
                face.right_cells[0].fs.copy_values_from(face.left_cells[0].fs);
                face.right_cells[1].fs.copy_values_from(face.left_cells[0].fs);
            }
        }
        } // end switch bc_north
    } // end foreach b
    return;
} // end apply_boundary_conditions()
