// simcore.d for the Lagrangian 1D Gas Dynamics, also known as L1d4.
// PA Jacobs
// 2020-04-08
//
module simcore;

import std.conv;
import std.stdio;
import std.string;
import std.json;
import std.file;
import std.datetime;
import std.algorithm;
import std.range;

import nm.schedule;
import util.json_helper;
import geom;
import gas;
import kinetics;
import gasdyn.gasflow;
import config;
import tube;
import gasslug;
import lcell;
import piston;
import valve;
import endcondition;
import misc;

__gshared static GasModel[] gmodels;
__gshared static uint overall_species_count;
__gshared static uint[][] overall_species_index;
__gshared static uint overall_modes_count;
__gshared static uint[][] overall_modes_index;
__gshared static ThermochemicalReactor[] reactors;
__gshared static Tube tube1;
__gshared static Piston[] pistons;
__gshared static Valve[] valves;
__gshared static GasSlug[] gasslugs;
__gshared static EndCondition[] ecs;
__gshared static Diaphragm[] diaphragms;

struct SimulationData {
    int step = 0;
    int halt_now = 0;
    double sim_time = 0.0;
    double dt_global;
    double cfl;
    double t_plot;
    double t_hist;
    int tindx;
    int steps_since_last_plot_write;
    int steps_since_last_hist_write;
    SysTime wall_clock_start;
}

__gshared static SimulationData sim_data;


void init_simulation(int tindx_start)
{
    sim_data.tindx = tindx_start;
    sim_data.wall_clock_start = Clock.currTime();
    //
    string dirName = L1dConfig.job_name;
    JSONValue jsonData = readJSONfile(dirName~"/config.json");
    // Now that we have parsed JSON data, proceed to update those config values.
    auto configData = jsonData["config"];
    L1dConfig.title = getJSONstring(configData, "title", "");
    L1dConfig.gas_model_files = getJSONstringarray(configData, "gas_model_files", []);
    L1dConfig.reaction_files_1 = getJSONstringarray(configData, "reaction_files_1", []);
    L1dConfig.reaction_files_2 = getJSONstringarray(configData, "reaction_files_2", []);
    L1dConfig.reacting = getJSONbool(configData, "reacting", false);
    L1dConfig.T_frozen = getJSONdouble(configData, "T_frozen", 300.0);
    if (L1dConfig.verbosity_level >= 1) {
        writeln("Config:");
        writefln("  title= \"%s\"", L1dConfig.title);
        writeln("  gas_model_files= ", L1dConfig.gas_model_files);
        writeln("  reaction_files_1= ", L1dConfig.reaction_files_1);
        writeln("  reaction_files_2= ", L1dConfig.reaction_files_2);
        writeln("  reacting= ", L1dConfig.reacting);
        writeln("  T_frozen= ", L1dConfig.T_frozen);
    }
    assert(L1dConfig.gas_model_files.length == L1dConfig.reaction_files_1.length &&
           L1dConfig.gas_model_files.length == L1dConfig.reaction_files_2.length,
           "Lengths of gas model and reaction file lists are inconsistent.");
    overall_species_count = 0;
    foreach (gi, fileName; L1dConfig.gas_model_files) {
        auto gm = init_gas_model(fileName);
        gmodels ~= gm;
        overall_species_index ~= array(iota(overall_species_count, overall_species_count+gm.n_species));
        overall_species_count += gm.n_species;
        overall_modes_index ~= array(iota(overall_modes_count, overall_modes_count+gm.n_modes));
        overall_modes_count += gm.n_modes;
        auto fn1 = L1dConfig.reaction_files_1[gi];
        auto fn2 = L1dConfig.reaction_files_2[gi];
        if (fn1.length > 0) {
            reactors ~= init_thermochemical_reactor(gm, fn1, fn2);
        } else {
            reactors ~= null;
        }
    }
    L1dConfig.max_time = getJSONdouble(configData, "max_time", 0.0);
    L1dConfig.max_step = getJSONint(configData, "max_step", 0);
    L1dConfig.dt_init = getJSONdouble(configData, "dt_init", 0.0);
    int n_cfl = getJSONint(configData, "n_cfl", 0);
    double[] dummy; foreach (i; 0 .. n_cfl) { dummy ~= 0.0; }
    double[] cfl_times = getJSONdoublearray(configData, "cfl_times", dummy);
    double[] cfl_values = getJSONdoublearray(configData, "cfl_values", dummy);
    L1dConfig.cfl_schedule = new Schedule!double(cfl_times, cfl_values);
    L1dConfig.cfl_count = getJSONint(configData, "cfl_count", 10);
    L1dConfig.print_count = getJSONint(configData, "print_count", 50);
    L1dConfig.x_order = getJSONint(configData, "x_order", 0);
    L1dConfig.t_order = getJSONint(configData, "t_order", 0);
    int n_dt_plot = getJSONint(configData, "n_dt_plot", 0);
    dummy.length = n_dt_plot; foreach (ref v; dummy) { v = 0.0; }
    double[] t_changes = getJSONdoublearray(configData, "t_change", dummy);
    double[] values = getJSONdoublearray(configData, "dt_plot", dummy);
    L1dConfig.dt_plot = new Schedule!double(t_changes, values);
    values[] = getJSONdoublearray(configData, "dt_hist", dummy);
    L1dConfig.dt_hist = new Schedule!double(t_changes, values);
    L1dConfig.hloc_n = getJSONint(configData, "hloc_n", 0);
    L1dConfig.hloc_x.length = L1dConfig.hloc_n;
    dummy.length = L1dConfig.hloc_n; foreach (ref v; dummy) { v = 0.0; }
    L1dConfig.hloc_x[] = getJSONdoublearray(configData, "hloc_x", dummy);
    L1dConfig.nslugs = getJSONint(configData, "nslugs", 0);
    L1dConfig.npistons = getJSONint(configData, "npistons", 0);
    L1dConfig.nvalves = getJSONint(configData, "nvalves", 0);
    L1dConfig.ndiaphragms = getJSONint(configData, "ndiaphragms", 0);
    L1dConfig.necs = getJSONint(configData, "necs", 0);
    if (L1dConfig.verbosity_level >= 1) {
        writeln("  max_time= ", L1dConfig.max_time);
        writeln("  max_step= ", L1dConfig.max_step);
        writeln("  dt_init= ", L1dConfig.dt_init);
        writeln("  cfl_schedule= ", L1dConfig.cfl_schedule);
        writeln("  cfl_count= ", L1dConfig.cfl_count);
        writeln("  print_count= ", L1dConfig.print_count);
        writeln("  x_order= ", L1dConfig.x_order);
        writeln("  t_order= ", L1dConfig.t_order);
        writeln("  dt_plot= ", L1dConfig.dt_plot);
        writeln("  dt_hist= ", L1dConfig.dt_hist);
        writeln("  hloc_x= ", L1dConfig.hloc_x);
        writeln("  nslugs= ", L1dConfig.nslugs);
        writeln("  npistons= ", L1dConfig.npistons);
        writeln("  nvalves= ", L1dConfig.nvalves);
        writeln("  ndiaphragms= ", L1dConfig.ndiaphragms);
        writeln("  necs= ", L1dConfig.necs);
    }
    //
    tube1 = new Tube(L1dConfig.job_name~"/tube.data");
    //
    foreach (i; 0 .. L1dConfig.nslugs) {
        auto myData = jsonData[format("slug_%d", i)];
        size_t indx = gasslugs.length;
        gasslugs ~= new GasSlug(indx, myData);
    }
    //
    foreach (i; 0 .. L1dConfig.npistons) {
        auto myData = jsonData[format("piston_%d", i)];
        size_t indx = pistons.length;
        pistons ~= new Piston(indx, myData);
    }
    //
    foreach (i; 0 .. L1dConfig.nvalves) {
        auto myData = jsonData[format("valve_%d", i)];
        size_t indx = valves.length;
        valves ~= new Valve(indx, myData);
    }
    //
    foreach (i; 0 .. L1dConfig.necs) {
        auto myData = jsonData[format("end_condition_%d", i)];
        string ecClass = getJSONstring(myData, "class", "");
        size_t indx = ecs.length;
        switch (ecClass) {
        case "Diaphragm":
            auto myDia = new Diaphragm(indx, myData);
            ecs ~= myDia;
            diaphragms ~= myDia;
            break;
        case "GasInterface":
            ecs ~= new GasInterface(indx, myData);
            break;
        case "FreeEnd":
            ecs ~= new FreeEnd(indx, myData);
            break;
        case "VelocityEnd":
            ecs ~= new VelocityEnd(indx, myData);
            break;
        case "PistonFace":
            ecs ~= new PistonFace(indx, myData);
            break;
        default:
            string msg = text("Unknown EndCondition: ", ecClass);
            throw new Exception(msg);
        }
    }
    // Work through the end conditions and make connections from the
    // connected gas slugs and pistons back to the end conditions.
    foreach (ec; ecs) {
        if (ec.slugL) {
            if (ec.slugL_end == End.L) {
                ec.slugL.ecL = ec;
            } else {
                ec.slugL.ecR = ec;
            }
        }
        if (ec.slugR) {
            if (ec.slugR_end == End.L) {
                ec.slugR.ecL = ec;
            } else {
                ec.slugR.ecR = ec;
            }
        }
        if (ec.pistonL) {
            if (ec.pistonL_face == End.L) {
                ec.pistonL.ecL = ec;
            } else {
                ec.pistonL.ecR = ec;
            }
        }
        if (ec.pistonR) {
            if (ec.pistonR_face == End.L) {
                ec.pistonR.ecL = ec;
            } else {
                ec.pistonR.ecR = ec;
            }
        }
    }
    // Work through dynamic components and read their initial state.
    foreach (i, s; gasslugs) {
        if (L1dConfig.verbosity_level >= 1) {
            writeln("Initial state of slug ", i);
        }
        string fileName = L1dConfig.job_name ~ format("/slug-%04d-faces.data", i);
        File fp = File(fileName, "r");
        s.read_face_data(fp, tindx_start);
        fp.close();
        fileName = L1dConfig.job_name ~ format("/slug-%04d-cells.data", i);
        fp = File(fileName, "r");
        s.read_cell_data(fp, tindx_start);
        fp.close();
        if (L1dConfig.verbosity_level >= 1) {
            LFace f = s.faces[$-1];
            LCell c = s.cells[$-1];
            writeln(format("  face[%d] x=%e area=%e", s.faces.length-1, f.x, f.area));
            writeln(format("  cell[%d] xmid=%e p=%e T=%e vel=%e",
                           s.cells.length-1, c.xmid, c.gas.p, c.gas.T, c.vel));
        }
    }
    foreach (i, p; pistons) {
        if (L1dConfig.verbosity_level >= 1) {
            writeln("Initial state of piston ", i);
        }
        string fileName = L1dConfig.job_name ~ format("/piston-%04d.data", i);
        File fp = File(fileName, "r");
        p.read_data(fp, tindx_start);
        fp.close();
        if (L1dConfig.verbosity_level >= 1) {
            writeln(format("  x=%e, vel=%e is_restrain=%s brakes_on=%s on_buffer=%s",
                           p.x, p.vel, p.is_restrain, p.brakes_on, p.on_buffer));
        }
    }
    foreach (i, dia; diaphragms) {
        if (L1dConfig.verbosity_level >= 1) {
            writeln("Initial state of diaphragm (at EndCondition index)", dia.indx);
        }
        string fileName = L1dConfig.job_name ~ format("/diaphragm-%04d.data", i);
        File fp = File(fileName, "r");
        dia.read_data(fp, tindx_start);
        fp.close();
        if (L1dConfig.verbosity_level >= 1) {
            writeln(format("  state=%s", dia.state));
        }
    }
    // Note that, for a restart, sim_time will generally be nonzero
    sim_data.dt_global = L1dConfig.dt_init;
    sim_data.sim_time = get_time_from_times_file(tindx_start);
    sim_data.cfl = L1dConfig.cfl_schedule.get_value(sim_data.sim_time);
    sim_data.t_plot = sim_data.sim_time + L1dConfig.dt_plot.get_value(sim_data.sim_time);
    sim_data.t_hist = sim_data.sim_time + L1dConfig.dt_hist.get_value(sim_data.sim_time);
    sim_data.steps_since_last_plot_write = 0;
    sim_data.steps_since_last_hist_write = 0;
    if (L1dConfig.verbosity_level >= 1) {
        // For reporting wall-clock time, convert with precision of milliseconds.
        auto elapsed_ms = (Clock.currTime() - sim_data.wall_clock_start).total!"msecs"();
        double elapsed_s = to!double(elapsed_ms)/1000;
        writefln("Finished initialization at WC=%.3f seconds", elapsed_s);
    }
    return;
} // end init_simulation()


void integrate_in_time()
{
    if (L1dConfig.verbosity_level >= 1) {
        writeln("Begin time integration to max_time=", L1dConfig.max_time);
    }
    sim_data.step = 0;
    append(L1dConfig.job_name~"/events.txt", format("t=%e Begin stepping\n", sim_data.sim_time));
    //
    // Main time loop.
    while (sim_data.sim_time <= L1dConfig.max_time &&
           sim_data.step <= L1dConfig.max_step &&
           sim_data.halt_now == 0) {
        if (L1dConfig.verbosity_level >= 2) {
            writeln("Begin time step ", sim_data.step+1);
        }
        // 1. Set the size of the time step.
        if ((sim_data.step % L1dConfig.cfl_count) == 0) {
            sim_data.cfl = L1dConfig.cfl_schedule.get_value(sim_data.sim_time);
            double dt_allowed = gasslugs[0].suggested_time_step(sim_data.cfl);
            foreach (i; 1 .. gasslugs.length) {
                dt_allowed = min(dt_allowed, gasslugs[i].suggested_time_step(sim_data.cfl));
            }
            if (dt_allowed < sim_data.dt_global) {
                // Reduce immediately.
                sim_data.dt_global = dt_allowed;
            } else {
                // Cautious increase, only if we have taken some steps.
                if (sim_data.step > 0) {
                    sim_data.dt_global += 0.5*(dt_allowed - sim_data.dt_global);
                }
            }
        }
        // 2.1 Check Piston and buffer.
        foreach (p; pistons) { p.check_for_buffer_interaction(sim_data.sim_time); }
        // 2.2 Update state of end conditions.
        foreach (ec; ecs) {
            auto dia = cast(Diaphragm) ec;
            if (dia) { dia.update_state(sim_data.sim_time); }
        }
        // 3. Record current state of dynamic components.
        foreach (p; pistons) { p.record_state(); }
        foreach (s; gasslugs) {
            s.compute_areas_and_volumes();
            s.encode_conserved();
            s.record_state();
        }
        int attempt_number = 0;
        bool step_failed;
        do {
            ++attempt_number;
            step_failed = false;
            try {
                // 4. Predictor update.
                // 4.1 Update the end conditions.
                foreach (ec; ecs) {
                    auto dia = cast(Diaphragm) ec;
                    if (dia) { dia.update_state(sim_data.sim_time); }
                }
                // 4.2 Update dynamic elements.
                foreach (s; gasslugs) {
                    s.time_derivatives(0, sim_data.sim_time);
                    s.predictor_step(sim_data.dt_global);
                    if (s.bad_cells() > 0) { throw new Exception("Bad cells"); }
                }
                foreach (p; pistons) {
                    p.time_derivatives(0);
                    p.predictor_step(sim_data.dt_global);
                }
            } catch (Exception e) {
                writefln("Predictor step failed e.msg=%s", e.msg);
                step_failed = true;
                foreach (p; pistons) { p.restore_state(); }
                foreach (s; gasslugs) { s.restore_state(); }
                sim_data.dt_global *= 0.2;
            }
            if (L1dConfig.t_order == 2 && !step_failed) {
                try {
                    // 5. Corrector update.
                    // 5.1 Update the end conditions.
                    foreach (ec; ecs) {
                        // [TODO] Diaphragms
                    }
                    // 5.2 Update dynamic elements.
                    foreach (s; gasslugs) {
                        s.time_derivatives(1, sim_data.sim_time+sim_data.dt_global);
                        s.corrector_step(sim_data.dt_global);
                        if (s.bad_cells() > 0) { throw new Exception("Bad cells"); }
                    }
                    foreach (p; pistons) {
                        p.time_derivatives(1);
                        p.corrector_step(sim_data.dt_global);
                    }
                } catch (Exception e) {
                    writefln("Corrector step failed e.msg=%s", e.msg);
                    step_failed = true;
                    foreach (p; pistons) { p.restore_state(); }
                    foreach (s; gasslugs) { s.restore_state(); }
                    sim_data.dt_global *= 0.2;
                    continue;
                }
            }
            if (L1dConfig.reacting && !step_failed) {
                foreach (s; gasslugs) { s.chemical_increment(sim_data.dt_global); }
            }
        } while (step_failed && (attempt_number <= 3));
        if (step_failed) {
            throw new Exception("Step failed after 3 attempts.");
        }
        // 6. Occasional console output.
        if (L1dConfig.verbosity_level >= 1 &&
            ((sim_data.step % L1dConfig.print_count) == 0)) {
            // For reporting wall-clock time, convert with precision of milliseconds.
            auto elapsed_ms = (Clock.currTime() - sim_data.wall_clock_start).total!"msecs"();
            double elapsed_s = to!double(elapsed_ms)/1000;
            double WCtFT = ((sim_data.sim_time > 0.0) && (sim_data.step > 0)) ?
                elapsed_s*(L1dConfig.max_time-sim_data.sim_time)/sim_data.dt_global/sim_data.step : 0.0;
            double WCtMS = (sim_data.step > 0) ?
                (elapsed_s*(L1dConfig.max_step-sim_data.step))/sim_data.step : 0.0;
            writefln("Step=%d t=%.3e dt=%.3e cfl=%.3f WC=%.1f WCtFT=%.1f WCtMS=%.1f",
                     sim_data.step, sim_data.sim_time, sim_data.dt_global,
                     sim_data.cfl, elapsed_s, WCtFT, WCtMS);
            stdout.flush();
        }
        // 7. Update time and (maybe) write solution.
        sim_data.step += 1;
        sim_data.sim_time += sim_data.dt_global;
        if (sim_data.sim_time >= sim_data.t_plot) {
            write_state_gasslugs_pistons_diaphragms();
            sim_data.t_plot += L1dConfig.dt_plot.get_value(sim_data.sim_time);
            sim_data.steps_since_last_plot_write = 0;
        } else {
            sim_data.steps_since_last_plot_write++;
        }
        if (sim_data.sim_time >= sim_data.t_hist) {
            write_data_at_history_locations(sim_data.sim_time);
            write_energies(sim_data.sim_time);
            sim_data.t_hist += L1dConfig.dt_hist.get_value(sim_data.sim_time);
            sim_data.steps_since_last_hist_write = 0;
        } else {
            sim_data.steps_since_last_hist_write++;
        }
    } // End main time loop.
    //
    append(L1dConfig.job_name~"/events.txt", format("t=%e Finished stepping\n", sim_data.sim_time));
    // Write a final time solution.
    if (sim_data.steps_since_last_plot_write > 0) {
        write_state_gasslugs_pistons_diaphragms();
    }
    if (sim_data.steps_since_last_hist_write > 0) {
        write_data_at_history_locations(sim_data.sim_time);
        write_energies(sim_data.sim_time);
    }
    return;
} // end integrate_in_time()


void write_state_gasslugs_pistons_diaphragms()
{
    sim_data.tindx += 1;
    if (L1dConfig.verbosity_level >= 1) {
        writeln("Write state data at tindx=", sim_data.tindx);
    }
    string fileName = L1dConfig.job_name ~ "/times.data";
    File fp = File(fileName, "a");
    fp.writefln("%d %e", sim_data.tindx, sim_data.sim_time);
    fp.close();
    foreach (i, s; gasslugs) {
        if (L1dConfig.verbosity_level >= 2) { writeln("  Writing state data for slug ", i); }
        fileName = L1dConfig.job_name ~ format("/slug-%04d-faces.data", i);
        fp = File(fileName, "a");
        s.write_face_data(fp, sim_data.tindx, false);
        fp.close();
        fileName = L1dConfig.job_name ~ format("/slug-%04d-cells.data", i);
        fp = File(fileName, "a");
        s.write_cell_data(fp, sim_data.tindx, false);
        fp.close();
    }
    foreach (i, p; pistons) {
        if (L1dConfig.verbosity_level >= 2) { writeln("  Writing state of piston ", i); }
        fileName = L1dConfig.job_name ~ format("/piston-%04d.data", i);
        fp = File(fileName, "a");
        p.write_data(fp, sim_data.tindx, false);
        fp.close();
    }
    foreach (i, dia; diaphragms) {
        if (L1dConfig.verbosity_level >= 2) {
            writeln("  Writing state of diaphragm (at EndCondition index)", dia.indx);
        }
        fileName = L1dConfig.job_name ~ format("/diaphragm-%04d.data", i);
        fp = File(fileName, "a");
        dia.write_data(fp, sim_data.tindx, false);
        fp.close();
    }
    return;
} // end write_state_gasslugs_pistons_diaphragms()


void write_data_at_history_locations(double t)
{
    foreach (i; 0 .. L1dConfig.hloc_n) {
        string fileName = L1dConfig.job_name ~ format("/history-loc-%04d.data", i);
        double x = L1dConfig.hloc_x[i];
        File fp = File(fileName, "a");
        foreach (s; gasslugs) { s.write_history_loc_data(fp, t, x); }
        fp.close();
    }
    return;
} // end write_data_at_history_locations()


void write_energies(double t)
{
    string fileName = L1dConfig.job_name ~ format("/energies.data");
    File fp = File(fileName, "a");
    fp.write(format("%e", t));
    double e_total = 0.0;
    foreach (s; gasslugs) {
        double e = s.energy;
        fp.write(format(" %e", e));
        e_total += e;
    }
    foreach (p; pistons) {
        double e = p.energy;
        fp.write(format(" %e", e));
        e_total += e;
    }
    fp.write(format(" %e\n", e_total));
    fp.close();
} // end write_energies()
