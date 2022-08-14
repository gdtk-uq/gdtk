// marching_calc.d
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

import json_helper;
import geom;
import config;
import streamtube;

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
    Config.title = getJSONstring(configData, "title", "");
    Config.gas_model_file = getJSONstring(configData, "gas_model_file", "");
    Config.reaction_file_1 = getJSONstring(configData, "reaction_files_1", "");
    Config.reaction_file_2 = getJSONstring(configData, "reaction_file_2", "");
    Config.reacting = getJSONbool(configData, "reacting", false);
    Config.T_frozen = getJSONdouble(configData, "T_frozen", 300.0);
    Config.axisymmetric = getJSONbool(configData, "axisymmetric", false);
    Config.max_x = getJSONdouble(configData, "max_x", 0.0);
    Config.max_step = getJSONint(configData, "max_step", 0);
    Config.dx = getJSONdouble(configData, "dx", 0.0);
    Config.max_step_relax = getJSONint(configData, "max_step_relax", 100);
    Config.cfl = getJSONdouble(configData, "cfl", 0.5);
    Config.print_count = getJSONint(configData, "print_count", 50);
    Config.plot_dx = getJSONdouble(configData, "plot_dx", 1.0e-2);
    Config.x_order = getJSONint(configData, "x_order", 2);
    Config.t_order = getJSONint(configData, "t_order", 2);
    Config.flux_calc = to!FluxCalcCode(getJSONint(configData, "flux_calc", 0));
    Config.compression_tol = getJSONdouble(configData, "compression_tol", -0.3);
    Config.shear_tol = getJSONdouble(configData, "shear_tol", 0.2);
    Config.n_streams = getJSONint(configData, "n_streams", 1);
    if (Config.verbosity_level >= 1) {
        writeln("Config:");
        writefln("  title= \"%s\"", Config.title);
        writeln("  gas_model_files= ", Config.gas_model_file);
        writeln("  reaction_files_1= ", Config.reaction_file_1);
        writeln("  reaction_files_2= ", Config.reaction_file_2);
        writeln("  reacting= ", Config.reacting);
        writeln("  T_frozen= ", Config.T_frozen);
        writeln("  axisymmetric= ", Config.axisymmetric);
        writeln("  max_x= ", Config.max_x);
        writeln("  max_step= ", Config.max_step);
        writeln("  dx= ", Config.dx);
        writeln("  max_step_relax= ", Config.max_step_relax);
        writeln("  cfl= ", Config.cfl);
        writeln("  print_count= ", Config.print_count);
        writeln("  plot_dx= ", Config.plot_dx);
        writeln("  x_order= ", Config.x_order);
        writeln("  t_order= ", Config.t_order);
        writeln("  flux_calc= ", Config.flux_calc);
        writeln("  compression_tol= ", Config.compression_tol);
        writeln("  shear_tol= ", Config.shear_tol);
        writeln("  n_streams= ", Config.n_streams);
    }
    foreach (i; 0 .. Config.n_streams) {
        streams ~= new StreamTube(i, configData);
        if (Config.verbosity_level >= 1) {
            writefln("  stream[%d]= %s", i, streams[i]);
        }
    }
    //
    foreach (st; streams) {
        st.set_up_data_storage();
        st.set_up_inflow_boundary();
        st.write_flow_data(true);
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
    while (progress.x < Config.max_x || progress.step < Config.max_step) {
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
                foreach (st; streams) { st.set_up_slice(progress.x + progress.dx); }
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
        foreach (st; streams) { st.shuffle_data_west(); }
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
            foreach (st; streams) { st.write_flow_data(false); }
            progress.steps_since_last_plot_write = 0;
            progress.plot_at_x += Config.plot_dx;
        } else {
            progress.steps_since_last_plot_write++;
        }
    } // end while
    //
    // Write the final slice, maybe.
    if (progress.steps_since_last_plot_write > 0) {
        foreach (st; streams) { st.write_flow_data(false); }
    }
    return;
} // end do_space_marching_calculation()

@nogc
void relax_slice_to_steady_flow(double xmid)
// We are operating on a slice with its mid-point being at x=xmid.
// There may be more than one StreamTube in the slice and any exchange
// boundary condition will necessarily involve two StreamTubes.
{
    foreach (j, st; streams) {
        st.encode_conserved(0);
        progress.dt_values[j] = st.estimate_allowable_dt();
    }
    double dt = progress.dt_values[0];
    foreach (j; 1 .. progress.dt_values.length) {
        dt = fmin(dt, progress.dt_values[j]);
    }
    //
    foreach (k; 0 .. Config.max_step_relax) {
        // 1. Predictor (Euler) step..
        apply_boundary_conditions(xmid);
        foreach (st; streams) { st.mark_shock_cells(); }
        foreach (st; streams) { st.predictor_step(dt); }
        if (Config.t_order > 1) {
            apply_boundary_conditions(xmid);
            foreach (st; streams) {
                st.corrector_step(dt);
                st.transfer_conserved_quantities(2, 0);
            }
        } else {
            // Clean-up after Euler step.
            foreach (st; streams) {
                st.transfer_conserved_quantities(1, 0);
            }
        }
        // 3. [TODO] measure residuals overall
        // break early, if residuals are small enough
    }
    return;
} // end relax_slice_to_steady_flow()

@nogc
void apply_boundary_conditions(double xmid)
// Look up boundary conditions at xmid and apply boundary conditions.
// Application of the boundary conditions is essentially filling the
// ghost-cell flow states with suitable data.
{
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
            // We fill the ghost-cells with data drom the neighbour streamtube.
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
