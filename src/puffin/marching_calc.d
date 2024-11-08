// marching_calc.d -- Part of the Puffin steady-flow calculator.
//
// PA Jacobs
// 2022-01-22
//
module marching_calc;

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
import std.parallelism;

import util.json_helper;
import geom;
import config;
import streamtube;
import kinetics.init_thermochemical_reactor;

// We use __gshared so that several threads may access
// the following array concurrently.
__gshared static StreamTube[] streams;

struct ProgressData {
    int step = 0;
    double x = 0.0;
    double dx = 0.0;
    double plot_at_x = 0.0;
    double[] dt_values;
    int steps_since_last_plot_write = 0;
    SysTime wall_clock_start;
}

__gshared static ProgressData progress;

void init_calculation()
{
    string dirName = Config.job_name;
    JSONValue configData = readJSONfile(dirName~"/config.json");
    parse_config_data_for_marching_solver(configData);
    foreach (i; 0 .. Config.n_streams) {
        streams ~= new StreamTube(i, configData);
        if (Config.verbosity_level > 1) {
            writefln("  stream[%d]= %s", i, streams[i]);
        }
    }
    //
    // Parallel runtime will have main thread plus extras.
    int extraThreadsInPool = min(Config.maxCPUs-1, Config.n_streams-1);
    defaultPoolThreads(extraThreadsInPool);
    //
    foreach (st; streams) {
        st.set_up_data_storage();
        st.set_up_inflow_boundary();
        if (st.active.get_value(0.0) == 1) {
            st.write_flow_data(true, true);
        } else {
            st.write_flow_data(true, false);
        }
        // [TODO] Should shift the initialization of the kinetics into stream module.
        if (Config.reacting) {
            st.thermochemUpdate = init_thermochemical_reactor(st.gmodel,
                                                              Config.reaction_file_1,
                                                              Config.reaction_file_2);
        }
    }
    progress.step = 0;
    progress.x = 0.0;
    progress.dx = Config.dx;
    progress.plot_at_x = Config.plot_dx;
    progress.steps_since_last_plot_write = 0;
    progress.dt_values.length = streams.length;
    return;
} // end init_calculation()

void do_space_marching_calculation()
{
    progress.wall_clock_start = Clock.currTime();
    while (progress.x < Config.max_x && progress.step < Config.max_step) {
        // 1. Set size of space step.
        if (progress.dx < Config.dx) { progress.dx *= 1.2; }
        progress.dx = min(Config.dx, progress.dx);
        // 2. Take a step.
        int attempt_number = 0;
        bool step_failed;
        do {
            ++attempt_number;
            step_failed = false;
            try {
                foreach (st; parallel(streams, 1)) { st.set_up_slice(progress.x + progress.dx); }
                relax_slice_to_steady_flow(progress.x + 0.5*progress.dx);
            } catch (Exception e) {
                writefln("Step failed e.msg=%s", e.msg);
                step_failed = true;
                progress.dx *= 0.2;
            }
        } while (step_failed && (attempt_number <= 3));
        if (step_failed) {
            throw new Exception("Step failed after 3 attempts.");
        }
        //
        // 3. Prepare for next spatial step.
        foreach (st; parallel(streams, 1)) { st.shuffle_data_west(); }
        progress.x += progress.dx;
        progress.step++;
        //
        // 4. Occasional console output.
        if (Config.verbosity_level >= 1 &&
            ((progress.step % Config.print_count) == 0)) {
            // For reporting wall-clock time, convert with precision of milliseconds.
            auto elapsed_ms = (Clock.currTime() - progress.wall_clock_start).total!"msecs"();
            double elapsed_s = to!double(elapsed_ms)/1000;
            double WCtFT = ((progress.x > 0.0) && (progress.step > 0)) ?
                elapsed_s*(Config.max_x-progress.x)/progress.dx/progress.step : 0.0;
            writefln("Step=%d x=%.3e dx=%.3e WC=%.1f WCtFT=%.1f",
                     progress.step, progress.x, progress.dx, elapsed_s, WCtFT);
            stdout.flush();
        }
        //
        // 5. Write a flow slice (maybe).
        if (progress.x >= progress.plot_at_x) {
            foreach (st; streams) {
                if (st.active.get_value(progress.x) == 1) {
                    st.write_flow_data(false, true);
                }
            }
            progress.steps_since_last_plot_write = 0;
            progress.plot_at_x += Config.plot_dx;
        } else {
            progress.steps_since_last_plot_write++;
        }
    } // end while
    //
    // Write the final slice, maybe.
    if (progress.steps_since_last_plot_write > 0) {
        foreach (st; streams) { st.write_flow_data(false, true); }
    }
    return;
} // end do_space_marching_calculation()


void relax_slice_to_steady_flow(double xmid)
// We are operating on a slice with its mid-point being at x=xmid.
// There may be more than one StreamTube in the slice and any exchange
// boundary condition will necessarily involve two StreamTubes.
{
    foreach (j, st; parallel(streams, 1)) {
        st.encode_conserved(0);
        progress.dt_values[j] = st.estimate_allowable_dt();
    }
    double dt = progress.dt_values.minElement();
    //
    foreach (k; 0 .. Config.max_step_relax) {
        // 1. Predictor (Euler) step..
        apply_boundary_conditions(xmid);
        foreach (st; parallel(streams, 1)) {
            if (st.active.get_value(xmid) == 1) {
                st.mark_shock_cells();
                st.predictor_step(dt, xmid);
            }
        }
        // 2. (possible) corrector step.
        if (Config.t_order > 1) {
            apply_boundary_conditions(xmid);
            foreach (st; parallel(streams, 1)) {
                if (st.active.get_value(xmid) == 1) {
                    st.corrector_step(dt, xmid);
                }
            }
        }
        // 3. Chemistry step is loosely coupled to the gas dynamics.
        if (Config.reacting) {
            foreach (st; parallel(streams, 1)) {
                if (st.active.get_value(xmid) == 1) {
                    st.thermochemical_increment(dt);
                }
            }
        }
        // 4. Prepare for next step.
        foreach (st; parallel(streams, 1)) {
            st.transfer_conserved_quantities(1, 0);
        }
        // 5. [TODO] measure residuals overall
        // break relaxation loop early, if residuals are small enough
    }
    return;
} // end relax_slice_to_steady_flow()

void apply_boundary_conditions(double xmid)
// Look up boundary conditions at xmid and apply boundary conditions.
// Application of the boundary conditions is essentially filling the
// ghost-cell flow states with suitable data.
{
    // We could make the following loop parallel, but it seems that the
    // amount of work per thread is too small to benefit.
    foreach (i, st; streams) {
        int bc0 = st.bc_lower.get_value(xmid);
        switch (bc0) {
        case BCCode.wall:
            // The slip-wall condition is implemented by filling the ghost cells
            // with reflected-normal-velocity flow.
            auto fstate = &(st.ghost_cells_left[0].fs);
            auto face = st.jfaces[0];
            fstate.copy_values_from(st.cells[0].fs);
            fstate.vel.transform_to_local_frame(face.n, face.t1);
            fstate.vel.x = -(fstate.vel.x);
            fstate.vel.transform_to_global_frame(face.n, face.t1);
            //
            fstate = &(st.ghost_cells_left[1].fs);
            fstate.copy_values_from(st.cells[1].fs);
            fstate.vel.transform_to_local_frame(face.n, face.t1);
            fstate.vel.x = -(fstate.vel.x);
            fstate.vel.transform_to_global_frame(face.n, face.t1);
            break;
        case BCCode.exchange:
            // We fill the ghost-cells with data drom the neighbour streamtube.
            auto nst = streams[i-1];
            st.ghost_cells_left[0].fs.copy_values_from(nst.cells[nst.ncells-1].fs);
            st.ghost_cells_left[1].fs.copy_values_from(nst.cells[nst.ncells-2].fs);
            break;
        default:
            throw new Exception("Unknown BCCode.");
        }
        int bc1 = st.bc_upper.get_value(xmid);
        switch (bc1) {
        case BCCode.wall:
            auto fstate = &(st.ghost_cells_right[0].fs);
            auto face = st.jfaces[st.ncells];
            fstate.copy_values_from(st.cells[st.ncells-1].fs);
            fstate.vel.transform_to_local_frame(face.n, face.t1);
            fstate.vel.x = -(fstate.vel.x);
            fstate.vel.transform_to_global_frame(face.n, face.t1);
            //
            fstate = &(st.ghost_cells_right[1].fs);
            fstate.copy_values_from(st.cells[st.ncells-2].fs);
            fstate.vel.transform_to_local_frame(face.n, face.t1);
            fstate.vel.x = -(fstate.vel.x);
            fstate.vel.transform_to_global_frame(face.n, face.t1);
            break;
        case BCCode.exchange:
            // We fill the ghost-cells with data from the neighbour streamtube.
            auto nst = streams[i+1];
            st.ghost_cells_right[0].fs.copy_values_from(nst.cells[0].fs);
            st.ghost_cells_right[1].fs.copy_values_from(nst.cells[1].fs);
            break;
        default:
            throw new Exception("Unknown BCCode.");
        }
    } // end foreach st
    return;
} // end apply_boundary_conditions()
