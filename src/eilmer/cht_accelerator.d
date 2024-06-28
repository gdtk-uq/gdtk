/**
 *
 * This program loosely couples Eilmer's fluid domain and solid domain solvers.
 * The coupling is via a Flux Forwad Temperature Back (FFTB) boundary condition.
 * Steady-state fluid domain updates are via the Newton-Krylov solver, and
 * transient (time-accurate) solid domain updates are via Super-Time-Stepping.
 *
 * Author: Kyle A. Damm
 * Date:   2022-10-02
 *
 * Note: This program is currently under development and is not intended for wider use yet.
 *
 */

import core.thread;
import core.runtime;

import std.stdio;
import std.format;
import std.conv;
import std.parallelism;
import std.algorithm;
import std.getopt;
import std.json;
import std.file;
import std.math;
import std.array;

import fileutil;
import json_helper;
import simcore_exchange;
import bc;
import flowstate;
import geom;
import steadystate_core;
import special_block_init;
import globaldata;
import globalconfig;
import simcore : init_simulation;
import postprocess : readTimesFile;
import loads;
import fluidblockio;
import fluidblockio_old;
import fluidblockio_new;
import solid_loose_coupling_update;

version(mpi_parallel) {
    import mpi;
}

int main(string[] args)
{
    int exitFlag;
    version(mpi_parallel) {
        // This preamble copied directly from the OpenMPI hello-world example.
        auto c_args = Runtime.cArgs;
        MPI_Init(&(c_args.argc), &(c_args.argv));
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        GlobalConfig.mpi_rank_for_local_task = rank;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        GlobalConfig.mpi_size = size;
        scope(success) { MPI_Finalize(); }
        // Make note that we are in the context of an MPI task, presumably, one of many.
        GlobalConfig.in_mpi_context = true;
        GlobalConfig.is_master_task = (GlobalConfig.mpi_rank_for_local_task == 0);
    } else {
        // We are NOT in the context of an MPI task.
        GlobalConfig.in_mpi_context = false;
        GlobalConfig.is_master_task = true;
    }

    if (GlobalConfig.is_master_task) {
        writeln("Eilmer 4.0 compressible-flow simulation code -- Conjugate Heat Transfer (CHT) accelerator.");
        writeln("Revision-id: PUT_REVISION_STRING_HERE");
        writeln("Revision-date: PUT_REVISION_DATE_HERE");
        writeln("Compiler-name: PUT_COMPILER_NAME_HERE");
        writeln("Build-date: PUT_BUILD_DATE_HERE");
        write("Build-flavour: ");
        version(flavour_debug) { writeln("debug"); }
        version(flavour_profile) { writeln("profile"); }
        version(flavour_fast) { writeln("fast"); }
    }

    string msg;
    version(mpi_parallel) {
        msg ~= "Usage: e4-cht-mpi [OPTIONS]\n";
    } else {
        msg ~= "Usage: e4-cht-shared [OPTIONS]\n";
    }
    msg ~= "OPTIONS include the following:\n";
    msg ~= "Option:                             Comment:\n";
    msg ~= "\n";
    msg ~= "  --job=<string>                    file names built from this string\n";
    msg ~= "  --verbosity=<int>                 defaults to 0\n";
    msg ~= "  --snapshot-start=<int>|last       defaults to 0\n";
    version(mpi_parallel) {
        msg ~= "  --threads-per-mpi-task=<int>      defaults to 1\n";
    } else {
        msg ~= "  --max-cpus=<int>                  defaults to ";
        msg ~= to!string(totalCPUs) ~" on this machine\n";
    }
    msg ~= "  --max-wall-clock=<int>            in seconds\n";
    msg ~= "  --help                            writes this message\n";

    if (args.length < 2) {
        if (GlobalConfig.is_master_task) {
            writeln("Too few arguments.");
            write(msg);
            stdout.flush();
        }
        exitFlag = 1;
        return exitFlag;
    }
    
    string jobName = "";
    int verbosityLevel = 0;
    string snapshotStartStr = "0";
    int snapshotStart = 0;
    int maxCPUs = totalCPUs;
    int threadsPerMPITask = 1;
    string maxWallClock = "432000"; // 5 days default
    bool helpWanted = false;
    try {
        getopt(args,
               "job", &jobName,
               "verbosity", &verbosityLevel,
               "snapshot-start", &snapshotStartStr,
               "max-cpus", &maxCPUs,
               "threads-per-mpi-task", &threadsPerMPITask,
               "max-wall-clock", &maxWallClock,
               "help", &helpWanted
               );
    } catch (Exception e) {
        if (GlobalConfig.is_master_task) {
            writeln("Problem parsing command-line options.");
            writeln("Arguments not processed: ");
            args = args[1 .. $]; // Dispose of program name in first argument
            foreach (arg; args) writeln("   arg: ", arg);
            write(msg);
            stdout.flush();
        }
        exitFlag = 1;
        return exitFlag;
    }
    if (verbosityLevel > 0) {
        if (GlobalConfig.is_master_task) {
            writeln("Begin simulation with command-line arguments.");
            writeln("  jobName: ", jobName);
            writeln("  snapshotStart: ", snapshotStartStr);
            writeln("  maxWallClock: ", maxWallClock);
            writeln("  verbosityLevel: ", verbosityLevel);
        }
    }
    version(mpi_parallel) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (GlobalConfig.is_master_task) {
            writefln("Parallelism: Distributed memory with message passing (MPI), number of tasks %d", GlobalConfig.mpi_size);
        }
        stdout.flush();
        // Give master_task a chance to be seen first.
        Thread.sleep(dur!("msecs")(100));
        MPI_Barrier(MPI_COMM_WORLD);
        // Now, get all tasks to report.
        char[256] processor_name;
        int len;
        MPI_Get_processor_name(processor_name.ptr, &len);
        if (verbosityLevel>1){
            writefln("MPI-parallel, start task %d on processor %s",
                     GlobalConfig.mpi_rank_for_local_task,
                     to!string(processor_name[0..len]));
            stdout.flush();
            Thread.sleep(dur!("msecs")(100));
            MPI_Barrier(MPI_COMM_WORLD);
        }
    } else {
        writeln("Parallelism: Shared memory");
    }
    if (helpWanted) {
        if (GlobalConfig.is_master_task) {
            write(msg);
            stdout.flush();
        }
        exitFlag = 0;
        return exitFlag;
    }

    if (jobName.length == 0) {
        if (GlobalConfig.is_master_task) {
            writeln("Need to specify a job name.");
            write(msg);
            stdout.flush();
        }
        exitFlag = 1;
        return exitFlag;
    }

    GlobalConfig.base_file_name = jobName;
    GlobalConfig.verbosity_level = verbosityLevel;
    maxCPUs = min(max(maxCPUs, 1), totalCPUs); // don't ask for more than available
    threadsPerMPITask = min(threadsPerMPITask, totalCPUs);

    switch (snapshotStartStr) {
    case "last":
        auto timesDict = readTimesFile(jobName);
        auto tindxList = timesDict.keys;
        sort(tindxList);
        snapshotStart = tindxList[$-1];
        break;
    default:
        snapshotStart = to!int(snapshotStartStr);
    }

    if (GlobalConfig.is_master_task) { writefln("Initialising simulation from snapshot: %d", snapshotStart); }
    init_simulation(snapshotStart, -1, maxCPUs, threadsPerMPITask, maxWallClock);

    // Check that the user has supplied a structured grid
    int goodToProceed = 1;
    foreach (blk; localFluidBlocks) {
        if (blk.grid_type == Grid_t.unstructured_grid) { goodToProceed = 0; }
    }
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &(goodToProceed), 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    }
    if (!goodToProceed) {
        if (GlobalConfig.is_master_task) {
            writeln("The CHT accelerator is not yet compatible with unstructured grids.");
            writeln("Bailing out!");
            stdout.flush();
        }
        exitFlag = 1;
        return exitFlag;
    }


    // verify_jacobian();

    // set some local config parameters
    bool init_precondition_matrix;
    bool solid_domain_only = GlobalConfig.sdluOptions.solidDomainOnly;
    SolidTimeIntegrationScheme solid_time_integration_scheme = GlobalConfig.sdluOptions.solidTimeIntegrationScheme;

    // warn users about setting explicit time integration scheme
    GasdynamicUpdate explicit_scheme = GlobalConfig.gasdynamic_update_scheme;
    if (solid_time_integration_scheme == SolidTimeIntegrationScheme.explicit &&
        (explicit_scheme != GasdynamicUpdate.rkl1 && explicit_scheme != GasdynamicUpdate.rkl2)) {
        if (GlobalConfig.is_master_task) {
            writeln("ERROR: the user has selected to solve the solid domain with an incompatible explicit update scheme.");
            writeln("       Please set either rkl1 or rkl2 as the config.gasdynamic_update_scheme parameter.");
            writeln("Bailing out!");
            stdout.flush();
        }
        exitFlag = 1;
        return exitFlag;
    }

    // Additional memory allocation specific to steady-state solver
    allocate_global_fluid_workspace();
    foreach (blk; localFluidBlocks) { blk.allocate_GMRES_workspace(); }
    if (solid_time_integration_scheme == SolidTimeIntegrationScheme.implicit) {
        if (GlobalConfig.is_master_task) {
            writeln("Allocating memory for Newton-Krylov solver in solid domain...");
        }
        allocate_global_solid_workspace();
        foreach (sblk; localSolidBlocks) { sblk.allocate_GMRES_workspace(); }
    }

    if (solid_domain_only) {
        // we will only run the solid domain solver
        init_precondition_matrix = true;
        double target_time = GlobalConfig.max_time;
        final switch (solid_time_integration_scheme) {
            case SolidTimeIntegrationScheme.explicit:
                integrate_solid_in_time_explicit(target_time);
                break;
            case SolidTimeIntegrationScheme.implicit:
                integrate_solid_in_time_implicit(target_time, init_precondition_matrix);
                break;
        } // end switch

        // write out a flow/solid solution
        write_cht_solution(jobName, target_time, 0, 0.0);

        if (GlobalConfig.is_master_task) {
            writeln("Done simulation.");
            stdout.flush();
        }

        exitFlag = 0;
        return exitFlag;
    }

    // else we continue with the fluid-solid coupled solver...
    // read .cht JSON file
    int npoints, n_startup_steps;
    bool constant_freestream, warm_start_fluid_solve;
    double dt_couple_max, dt_couple_init;
    JSONValue jsonData     = readJSONfile(jobName~".cht");
    npoints                = to!int(jsonData["npoints"].integer);
    constant_freestream    = jsonData["constant_freestream"].boolean;
    warm_start_fluid_solve = jsonData["warm_start_fluid_solve"].boolean;
    dt_couple_max          = jsonData["dt_couple_max"].floating;
    dt_couple_init         = jsonData["dt_couple_init"].floating;
    n_startup_steps        = to!int(jsonData["n_startup_steps"].integer);

    if (n_startup_steps > npoints) {
        if (GlobalConfig.is_master_task) {
            writeln("The user has requested more startup steps than points in the trajectory.");
            writeln("Bailing out!");
            stdout.flush();
        }
        exitFlag = 1;
        return exitFlag;
    }

    // bootstrap the coupled boundary condition by sending the initial solid temperature to the fluid domain
    send_solid_domain_boundary_temperature_data_to_gas_domain();

    // run CHT simulation
    double time = 0.0;
    double dt_couple = dt_couple_init;
    double dt_factor;
    if (n_startup_steps == 0) {
        dt_factor = 1.0;
        dt_couple_max = dt_couple_init;
    }
    else {
        dt_factor = pow(dt_couple_max/dt_couple_init, 1.0/n_startup_steps);
    }

    // we alternate between fluid and solid solves...
    foreach (idx; 0..npoints) {

        int io_idx = idx;

        // we only need to initialise the precondition matrix for the steady-state solver on the first iteration
        init_precondition_matrix = (idx == 0) ? true: false;

        // if we have a constant freestream then we only have one entry (point_0) in the .cht file
        if (constant_freestream) { idx = 0; }

        if (GlobalConfig.is_master_task) {
            writefln("#################################################");
            writefln("### Simulating point %d at t = %.5f seconds ###", idx, time);
            writefln("#################################################");
        }
        time += dt_couple;

        // set the inflow condition for this point
        if (!constant_freestream) {
            FlowState inflow = FlowState(jsonData["point_"~to!string(idx)]["flowstate"], GlobalConfig.gmodel_master);
            set_inflow_condition(inflow);
        }

        // set the initial condition for this point if not warm starting with previous flow solution
        if (!warm_start_fluid_solve) {
            FlowState initial = FlowState(jsonData["point_"~to!string(idx)]["flowstate"], GlobalConfig.gmodel_master);
            set_initial_condition(initial);
        }

        // fluid domain solver
        iterate_to_steady_state(snapshotStart, maxCPUs, threadsPerMPITask, init_precondition_matrix);
        if (GlobalConfig.is_master_task) {
            // save a copy of the steady-state solver diagnostics file
            string file_name = "e4-nk.diagnostics.dat";
            if (exists(file_name)) {
                rename(file_name, to!string(io_idx)~"_"~file_name);
            } else {
                writef("WARNING: Could not find %s, no copy has been saved. \n", file_name);
            }
        }

        // write loads file
        write_loads(io_idx);

        // transfer heat flux data from fluid domain to solid domain (Flux Forward)
        send_gas_domain_boundary_heat_flux_data_to_solid_domain();

        // solid domain solver
        final switch (solid_time_integration_scheme) {
            case SolidTimeIntegrationScheme.explicit:
                integrate_solid_in_time_explicit(dt_couple);
                break;
            case SolidTimeIntegrationScheme.implicit:
                integrate_solid_in_time_implicit(dt_couple, init_precondition_matrix);
                break;
        } // end switch

        // transfer temperature data from solid domain to fluid domain (Temperature Back)
        send_solid_domain_boundary_temperature_data_to_gas_domain();

        // write out a flow/solid solution
        write_cht_solution(jobName, time, io_idx, dt_couple);

        // increase dt_couple
        dt_couple = fmin(dt_couple_max, dt_couple*dt_factor);
    }

    // Perform RHS evaluation to fill the fluid/solid BC with the latest solid temperature...
    steadystate_core.evalRHS(0.0, 0);

    // ... and write out the final loads file
    // NOTE: the heat flux in this loads file will not be physical
    write_loads(npoints);

    if (GlobalConfig.is_master_task) {
        writeln("Done simulation.");
        stdout.flush();
    }

    exitFlag = 0;
    return exitFlag;
}

void write_loads(int idx) {
    // helper function that writes out the loads file
    int dir_id = 1_000_000 + idx; // TODO: temporary hack, we need a more elegant solution. KAD 2022-11-08
    if (GlobalConfig.is_master_task) { init_current_loads_tindx_dir(dir_id); }
    version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }
    wait_for_current_tindx_dir(dir_id);
    write_boundary_loads_to_file(SimState.time, dir_id);
    if (GlobalConfig.is_master_task) { update_loads_times_file(SimState.time, dir_id); }
}


void set_initial_condition(FlowState initial) {
    // helper function that copies the provided FlowState into all of the cell flowstates
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach (cell; blk.cells) {
            cell.fs.copy_values_from(initial);
            blk.myConfig.gmodel.update_thermo_from_pT(cell.fs.gas);
            cell.encode_conserved(0, 0, blk.omegaz);
            cell.decode_conserved(0, 0, blk.omegaz);
        }
    }

    // we need to apply the special blending if required
    if (GlobalConfig.diffuseWallBCsOnInit) {
        foreach (myblk; parallel(localFluidBlocks,1)) {
            diffuseWallBCsIntoBlock(myblk, GlobalConfig.nInitPasses, GlobalConfig.initTWall);
        }
    }

    // we also switch off any flags that may have been turned on during the previous flow simulation here
    GlobalConfig.frozen_limiter = false; // this is a formality, since it won't have any effect for structured grids
    GlobalConfig.frozen_shock_detector = false;
    return;
}

void set_inflow_condition(FlowState inflow) {
    // helper function that copies the provided FlowState into all inflow_supersonic BCs
    foreach (blk; parallel(localFluidBlocks,1)) {
        foreach (boundary; blk.bc) {
            if (boundary.type == "inflow_supersonic") {
                auto inflowBC0 = cast(GhostCellFlowStateCopy) boundary.preReconAction[0];
                assert(inflowBC0 !is null, "Oops, this should be a GhostCellFlowStateCopy object.");
                inflowBC0.fstate.copy_values_from(inflow);
                auto inflowBC1 = cast(BIE_FlowStateCopy) boundary.preSpatialDerivActionAtBndryFaces[0];
                assert(inflowBC1 !is null, "Oops, this should be a BIE_FlowStateCopy object.");
                inflowBC1.fstate.copy_values_from(inflow);
            }
        }
    }

    return;
}

void write_cht_solution(string jobName, double time, int idx, double dt) {
    // helper function to write a flow and solid solution
    int tindx = 1_000_000 + idx; // TODO: temporary hack, we need a more elegant solution. KAD 2022-11-18
    FluidBlockIO[] io_list = localFluidBlocks[0].block_io;
    bool legacy = is_legacy_format(GlobalConfig.flow_format);
    if (GlobalConfig.is_master_task){
        if (legacy) {
            ensure_directory_is_present(make_path_name!"flow"(tindx));
        } else {
            foreach(io; io_list) {
                string path = "CellData/"~io.tag;
                if (io.do_save()) ensure_directory_is_present(make_path_name(path, tindx));
            }
        }
        ensure_directory_is_present(make_path_name!"solid"(tindx));
    }
    version(mpi_parallel) {
        MPI_Barrier(MPI_COMM_WORLD);
    }
    foreach (blk; localFluidBlocks) {
        if (legacy) {
            auto fileName = make_file_name!"flow"(jobName, blk.id, tindx, GlobalConfig.flowFileExt);
            blk.write_solution(fileName, time);
        } else {
            foreach(io; blk.block_io) {
                auto fileName = make_file_name("CellData", io.tag, jobName, blk.id, tindx, GlobalConfig.flowFileExt);
                if (io.do_save()) io.save_to_file(fileName, time);
            }
        }
    }
    foreach (sblk; localSolidBlocks) {
        auto fileName = make_file_name!"solid"(jobName, sblk.id, tindx, "gz");
        sblk.writeSolution(fileName, time);
    }

    // Update times file
    if (GlobalConfig.is_master_task) {
        auto writer = appender!string();
        formattedWrite(writer, "%04d %.18e %.18e\n", tindx, time, dt);
        append("config/" ~ GlobalConfig.base_file_name ~ ".times", writer.data);
    }

}
