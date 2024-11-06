/**
 * Core set of functions and structures used to advance the flow field in time.
 *
 * Authors: PAJ and RJG
 * Date: 2024-02-07
 * History:
 *   2024-02-07
 */

module timemarching;

import std.stdio : write, writeln, writefln, stdout, File;
import std.datetime : DateTime, Clock;
import core.memory : GC;
import std.conv : to;
import std.math : FloatingPointControl;
import std.parallelism : parallel;
import std.algorithm : min;
import std.file;
import std.math;
import std.format : format, formattedWrite;
import std.array : appender, split;
import std.string : splitLines;
import dyaml;

import ntypes.complex;
import util.time_utils : timeStringToSeconds;
import nm.number : number;
import conservedquantities : new_ConservedQuantities;
import simcore : check_run_time_configuration,
       synchronize_corner_coords_for_all_blocks,
       call_UDF_at_timestep_start,
       call_UDF_at_timestep_end,
       call_UDF_at_write_to_file,
       update_ch_for_divergence_cleaning,
       set_mu_and_k,
       chemistry_step,
       compute_Linf_residuals,
       compute_L2_residual,
       compute_mass_balance;
import simcore_gasdynamic_step : sts_gasdynamic_explicit_increment_with_fixed_grid,
       gasdynamic_explicit_increment_with_fixed_grid,
       gasdynamic_implicit_increment_with_fixed_grid,
       gasdynamic_explicit_increment_with_moving_grid,
       gasdynamic_implicit_increment_with_moving_grid,
       determine_time_step_size,
       local_dt_allow,
       local_dt_allow_parab,
       local_cfl_max,
       local_invalid_cell_count;
import simcore_solid_step : determine_solid_time_step_size, solid_step;

import lmrexceptions;
import lmrconfig;

import globalconfig;
import globaldata;
import init;
import fileutil : ensure_directory_is_present;
import lmr.loads : computeRunTimeLoads,
       init_current_loads_indx_dir,
       wait_for_current_indx_dir,
       writeLoadsToFile,
       count_written_loads,
       update_loads_metadata_file;
import lmr.fvcell : FVCell;
import lmr.history : initHistoryCells, writeHistoryCellsToFiles;

version(mpi_parallel) {
    import mpi;
}

void initTimeMarchingSimulation(int snapshotStart, int maxCPUs, int threadsPerMPITask, string maxWallClock)
{
    alias cfg = GlobalConfig;
    if (cfg.verbosity_level > 0 && cfg.is_master_task) {
        writeln("lmr run: Begin initTimeMarchingSimulation()...");
    }

    SimState.is_restart = (snapshotStart > 0);
    SimState.maxWallClockSeconds = timeStringToSeconds(maxWallClock);
    SimState.wall_clock_start = Clock.currTime();

    // Initialise baseline configuration
    initConfiguration();
    if (cfg.nFluidBlocks == 0 && cfg.is_master_task) {
        throw new NewtonKrylovException("No FluidBlocks; no point in continuing with simulation initialisation.");
    }
    // and read the control because it sets up part of initial configuration
    readControl();
    // [TODO] RJG, 2024-02-07
    // Implement opening of progress file.

    if (GlobalConfig.writeTransientResiduals && GlobalConfig.is_master_task && (snapshotStart == 0)) {
        std.file.write(lmrCfg.transResidFile, "# step time wall-clock mass x-mom y-mom z-mom energy L2 mass-balance\n");
    }

    SimState.current_tindx = snapshotStart;

    if (SimState.current_tindx == 0) {
        // On a fresh start, set loads tindx to 0
        SimState.current_loads_tindx = 0;
    }
    else {
        // We'll need to work a little harder and find out how many loads we've already written.
        SimState.current_loads_tindx = count_written_loads();
    }

    initLocalBlocks();
    initThreadPool(maxCPUs, threadsPerMPITask);
    initFluidBlocksBasic(true);
    initFluidBlocksMemoryAllocation();
    // [TODO] RJG, 2024-04-07
    // Add Lachlan's FSI initialisation here.
    // Well maybe not here exactly, but think on where to do this.
    initFluidBlocksGlobalCellIDStarts();
    initFluidBlocksZones();
    initFluidBlocksFlowField(snapshotStart);

    initFullFaceDataExchange();
    initMappedCellDataExchange();
    initGhostCellGeometry();
    initLeastSquaresStencils();

    if ((cfg.interpolation_order > 1) &&
        ((cfg.unstructured_limiter == UnstructuredLimiter.hvenkat) ||
         (cfg.unstructured_limiter == UnstructuredLimiter.venkat) ||
         (cfg.unstructured_limiter == UnstructuredLimiter.hvenkat_mlp) ||
         (cfg.unstructured_limiter == UnstructuredLimiter.venkat_mlp))) {
        initUSGlimiters();
    }

    // [TODO] RJG, 2024-04-07
    // Here is where initialise GPU chemistry, if we are going to continue
    // with that particular implementation.

    if (cfg.nSolidBlocks > 0) {
        initSolidBlocks();
        initFluidSolidExchangeBoundaries();
    }

    initHistoryCells();

    if (GlobalConfig.write_loads && (SimState.current_loads_tindx == 0)) {
        if (GlobalConfig.is_master_task) {
            initLoadsFiles();
        }
    }

    // For a simulation with shock fitting, the files defining the rails for
    // vertex motion and the scaling of vertex velocities throughout the blocks
    // will have been written by prep.lua + output.lua.
    if (GlobalConfig.grid_motion == GlobalConfig.grid_motion.shock_fitting) {
        foreach (i, fba; fluidBlockArrays) {
            if (fba.shock_fitting) {
                version(mpi_parallel) {
                    // The MPI tasks associated with this FBArray will have
                    // their own communicator for synchronizing the content
                    // of their shock-fitting arrays.
                    // We don't care about the rank of each task within
                    // that communicator.
                    MPI_Comm_split(MPI_COMM_WORLD, to!int(i), 0, &(fba.mpicomm));
                }
                // FIX-ME 2024-02-28 PJ Make use of Rowan's config information to find files.
                fba.read_rails_file(format("lmrsim/grid/gridarray-%04d.rails", i));
                fba.read_velocity_weights(format("lmrsim/grid/gridarray-%04d.weights", i));
            }
        }
    }

    // [TODO] RJG, 2024-02-07
    // When we have solid domains in place, we'll need to implement the initialisation of
    // the cell mapping between fluid and solid domains.

    orderBlocksBySize();
    // Flags to indicate that the saved output is fresh.
    // On startup or restart, it is assumed to be so.
    SimState.output_just_written = true;
    SimState.history_just_written = true;
    SimState.loads_just_written = true;
    // When starting a new calculation,
    // set the global time step to the initial value.
    // If we have elected to run with a variable time step,
    // this value may get revised on the very first step.
    SimState.dt_global = GlobalConfig.dt_init;
    // dt_global will be updated if restart branch is selected below.

    if (!SimState.is_restart) {
        SimState.time = 0.0;
        SimState.step = 0;
        if (cfg.is_master_task) {
            // Clean out any existing times file.
            if (lmrCfg.timesFile.exists) lmrCfg.timesFile.remove;
            // We can put an initial entry in the times file now.
            addToTimesFile();
        }
    }
    else {
        // SimState.time, SimState.step and SimState.dt_global set in prepareTimesFileOnRestart
        // It's convenient to do that while we have the times file loaded.
        if (cfg.is_master_task) {
            // RJG, 2024-11-06
            // For large-scale work (many processes on parallel filesystems),
            // we need to be careful about attempting to read one small file from many processes.
            // The solution is to read only from master, then broadcast to other ranks.
            // This is the same solution NNG implemented in the steady-state codepath for
            // broadcasting restart information.
            //
            // This issue was noted by Sebastiaan van Oeveren when working on the UQ Bunya cluster.
            prepareTimesFileOnRestart(SimState.current_tindx);
            writeln("*** RESTARTING SIMULATION ***");
            writefln("RESTART-SNAPSHOT: %d", SimState.current_tindx);
            writefln("RESTART-STEP: %d", SimState.step);
            writefln("RESTART-TIME: %8.3e", SimState.time);
        }
        version(mpi_parallel) {
            // After setting SimState on master, broadcast to other ranks
            // 0. Make non-shared temporaries
            double simTime = SimState.time;
            int step = SimState.step;
            double dtGlobal = SimState.dt_global;
            // 1. Broadcast
            MPI_Bcast(&simTime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Bcast(&step, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&dtGlobal, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            // 2. Write to shared variables
            SimState.time = simTime;
            SimState.step = step;
            SimState.dt_global = dtGlobal;
        }
    }

    initMasterLuaState();
    initCornerCoordinates();
    if (cfg.turb_model.needs_dwall) initWallDistances();

    if (cfg.compute_run_time_loads) initRunTimeLoads();

    // [TODO] RJG, 2024-02-07
    // Configure Nick's electric field solver.

    // Keep our memory foot-print small.
    if (GlobalConfig.verbosity_level >= 2) { writeln("Before GC.collect."); }
    GC.collect();
    GC.minimize();
    if (GlobalConfig.verbosity_level >= 2) { writeln("After GC.collect."); }
    //
    version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }

    // Close out by computing some memory usage and reporting
    auto myStats = GC.stats();
    double heapUsed = to!double(myStats.usedSize)/(2^^20);
    double heapFree = to!double(myStats.freeSize)/(2^^20);
    double minTotal = heapUsed+heapFree;
    double maxTotal = heapUsed+heapFree;

    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE,&heapUsed,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&heapFree,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&minTotal,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE,&maxTotal,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    }

    if (GlobalConfig.is_master_task) {
        writefln("Heap memory used: %.0f MB, unused: %.0f MB, total: %.0f MB (%.0f-%.0f MB per task)",
                 heapUsed, heapFree, heapUsed+heapFree, minTotal, maxTotal);
        stdout.flush();
    }
    version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }
    if (GlobalConfig.verbosity_level > 0 && GlobalConfig.is_master_task) {
        // For reporting wall-clock time, convert to seconds with precision of milliseconds.
        double wall_clock_elapsed = to!double((Clock.currTime() - SimState.wall_clock_start).total!"msecs"())/1000.0;
        writefln("Done init_simulation() at wall-clock(WC)= %.1f sec", wall_clock_elapsed);
        stdout.flush();
    }
} // end initTimeMarchinSimulation()

/**
 * Large coordinating function to advance the flow field via time integration.
 *
 * Authors: PAJ and RJG
 * Date: 2024-02-07
 * History:
 *   2024-02-07 - brought across from Eilmer4, only I/O changed.
 */
int integrateInTime(double targetTimeAsRequested)
{
    number mass_balance = to!number(0.0);
    number L2_residual = to!number(0.0);
    auto Linf_residuals = new_ConservedQuantities(GlobalConfig.cqi.n);
    version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }
    if (GlobalConfig.verbosity_level > 0 && GlobalConfig.is_master_task) {
        writeln("Integrate in time.");
        stdout.flush();
    }
    SimState.target_time = (GlobalConfig.block_marching) ? targetTimeAsRequested : GlobalConfig.max_time;
    // The next time for output...
    SimState.t_plot = SimState.time + GlobalConfig.dt_plot_schedule.get_value(SimState.time);
    SimState.t_history = SimState.time + GlobalConfig.dt_history;
    SimState.t_loads = SimState.time + GlobalConfig.dt_loads;
    //
    if (GlobalConfig.viscous) {
        // We have requested viscous effects but their application may be delayed
        // until the flow otherwise starts.
        if (GlobalConfig.viscous_delay > 0.0 && SimState.time < GlobalConfig.viscous_delay) {
            // We will initially turn-down viscous effects and
            // only turn them up when the delay time is exceeded.
            GlobalConfig.viscous_factor = 0.0;
        } else {
            // No delay in applying full viscous effects.
            GlobalConfig.viscous_factor = 1.0;
        }
    } else {
        // We haven't requested viscous effects at all.
        GlobalConfig.viscous_factor = 0.0;
    }
    foreach (myblk; localFluidBlocksBySize) {
        myblk.myConfig.viscous_factor = GlobalConfig.viscous_factor;
    }
    //
    // Select the actual gasdynamic update function.
    // These functions are sitting in module simcore_gasdynamic_step.
    void function() gasdynamic_step;
    if (GlobalConfig.grid_motion == GridMotion.none) {
        // Fixed grid
        if (GlobalConfig.with_super_time_stepping) {
            gasdynamic_step = &sts_gasdynamic_explicit_increment_with_fixed_grid;
        } else if (is_explicit_update_scheme(GlobalConfig.gasdynamic_update_scheme)) {
            gasdynamic_step = &gasdynamic_explicit_increment_with_fixed_grid;
        } else {
            gasdynamic_step = &gasdynamic_implicit_increment_with_fixed_grid;
        }
    } else {
        // Moving Grid
        if (is_explicit_update_scheme(GlobalConfig.gasdynamic_update_scheme)) {
            gasdynamic_step = &gasdynamic_explicit_increment_with_moving_grid;
        } else {
            gasdynamic_step = &gasdynamic_implicit_increment_with_moving_grid;
        }
    }
    if (!gasdynamic_step) {
        throw new Error("Did not set a valid gasdynamic_step function.");
    }
    //
    local_dt_allow.length = localFluidBlocks.length+localSolidBlocks.length; // prepare array for use later (we use this for solid blocks as well)
    local_cfl_max.length = localFluidBlocks.length; // prepare array for use later
    local_dt_allow_parab.length = localFluidBlocks.length;
    local_invalid_cell_count.length = localFluidBlocks.length;
    //
    // Normally, we can terminate upon either reaching
    // a maximum time or upon reaching a maximum iteration count.
    shared bool finished_time_stepping = (SimState.time >= SimState.target_time) ||
        (SimState.step >= GlobalConfig.max_step);
    // There are occasions where things go wrong and an exception may be caught.
    bool caughtException = false;
    //----------------------------------------------------------------
    //                 Top of main time-stepping loop
    //----------------------------------------------------------------
    // When SolidBlocks are in play,
    // we let the fluid time step settle before loosely coupling the solid domain;
    // during this period the fluid and solid are tightly coupled.
    // 1000 iterations appears to be sufficient.
    int solid_domain_loose_coupling_delay = 1000;
    int update_solid_domain_on_step = solid_domain_loose_coupling_delay;
    auto coupling_with_solid_domains_save = GlobalConfig.coupling_with_solid_domains;
    if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.loose) { GlobalConfig.coupling_with_solid_domains = SolidDomainCoupling.tight; }
    //
    SimState.wall_clock_start = Clock.currTime();
    while ( !finished_time_stepping ) {
        try {
            if (SimState.step == solid_domain_loose_coupling_delay &&
                coupling_with_solid_domains_save == SolidDomainCoupling.loose) {
                // switch to loose coupling
                GlobalConfig.coupling_with_solid_domains = SolidDomainCoupling.loose;
            }
            //
            // 0.0 Run-time configuration may change, a halt may be called, etc.
            check_run_time_configuration(targetTimeAsRequested);
            if (GlobalConfig.grid_motion != GridMotion.none) { synchronize_corner_coords_for_all_blocks(); }
            //
            // 1.0 Maintain a stable time step size, and other maintenance, as required.
            // The user may also have some start-of-time-step maintenance to do
            // via their Lua script file.  Let them have first go.
            if (GlobalConfig.udf_supervisor_file.length > 0) {
                // Note that the following call allows the user to do almost anything
                // at the start of the time step, including changing flow states in cells.
                call_UDF_at_timestep_start();
                // If the user has adjusted any of the flow states via the Lua functions,
                // we will beed to re-encode the conserved quantities, so that they have
                // consistent data.
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    foreach (cell; blk.cells) { cell.encode_conserved(0, 0, blk.omegaz); }
                }
            }
            if (!GlobalConfig.fixed_time_step) { determine_time_step_size(); }
            if (GlobalConfig.divergence_cleaning && !GlobalConfig.MHD_static_field) { update_ch_for_divergence_cleaning(); }
            // If using k-omega, we need to set mu_t and k_t BEFORE we call convective_update
            // because the convective update is where the interface values of mu_t and k_t are set.
            // only needs to be done on initial step, subsequent steps take care of setting these values
            if ((SimState.step == 0) && GlobalConfig.turb_model.isTurbulent) {
                set_mu_and_k();
            }
            //
            // 2.0 Attempt a time step.
            // 2.1 Chemistry 1/2 step, if appropriate.
            if (GlobalConfig.strangSplitting == StrangSplittingMode.half_R_full_T_half_R &&
                GlobalConfig.with_local_time_stepping) {
                throw new Error("StrangSplitting.half_R_full_T_half_R and LTS aren't currently compatible");
            }
            if (GlobalConfig.reacting &&
                GlobalConfig.chemistry_update == ChemistryUpdateMode.split &&
                (GlobalConfig.strangSplitting == StrangSplittingMode.half_R_full_T_half_R) &&
                (SimState.time > GlobalConfig.reaction_time_delay)) {
                double dt_chem = 0.5*SimState.dt_global;
                dt_chem *= GlobalConfig.reaction_fraction_schedule.interpolate_value(SimState.time);
                chemistry_step(dt_chem);
            }
            //
            // 2.2 Step the gasdynamic processes.
            if (SimState.step >= GlobalConfig.freeze_limiter_on_step && !(GlobalConfig.frozen_limiter)) {
                // Freeze the limiter at this point to help alleviate any ringing of the residuals.
                GlobalConfig.frozen_limiter = true;
            }
            gasdynamic_step();
            //
            // 2.3 Solid domain update (if loosely coupled)
            // If tight coupling, then this has already been performed in the gasdynamic_step().
            if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.loose &&
                SimState.step == update_solid_domain_on_step) {
                // Determine stable time step in solid domain.
                double dt_solid_stable = determine_solid_time_step_size();
                int n_solid_coupling = to!int(min((floor(dt_solid_stable/SimState.dt_global)), int.max));
                double dt_solid = n_solid_coupling*SimState.dt_global;
                update_solid_domain_on_step = SimState.step + n_solid_coupling;
                if (GlobalConfig.is_master_task) {
                    writeln("-- SOLID DOMAIN UPDATE: cfl=", GlobalConfig.solid_domain_cfl,
                            ", dt_solid=", dt_solid,
                            ", next update on step=", update_solid_domain_on_step);
                }
                solid_step(dt_solid);
            }
            //
            // 2.4 Chemistry step or 1/2 step (if appropriate).
            if (GlobalConfig.reacting &&
                GlobalConfig.chemistry_update == ChemistryUpdateMode.split &&
                (SimState.time > GlobalConfig.reaction_time_delay)) {
                double dt_chem = (GlobalConfig.strangSplitting == StrangSplittingMode.full_T_full_R) ?
                    SimState.dt_global : 0.5*SimState.dt_global;
                dt_chem *= GlobalConfig.reaction_fraction_schedule.interpolate_value(SimState.time);
                chemistry_step(dt_chem);
            }
            // 2.5 Update electric field solution (if needed)
            if (GlobalConfig.solve_electric_field && ((SimState.step+1)%GlobalConfig.electric_field_count==0)){
                if (GlobalConfig.is_master_task) writeln("Called field.solve_efield(): ...");
                eField.solve_efield(localFluidBlocks, GlobalConfig.is_master_task);
                eField.compute_electric_field_vector(localFluidBlocks);

                double current_in, current_out;
                eField.compute_boundary_current(localFluidBlocks, current_in, current_out);
                if (GlobalConfig.is_master_task) {
                    writeln("Called field.compute_boundary_current() ...");
                    writefln("    Current in:  %f (A/m)", current_in);
                    writefln("    Current out: %f (A/m)", current_out);
                }
            }

            // 3.0 Update the time record and (occasionally) print status.
            SimState.step = SimState.step + 1;
            if (GlobalConfig.is_master_task) {
                try {
                    std.file.write(lmrCfg.progFile, format("%d\n", SimState.step));
                } catch (Exception e) {
                    // do nothing
                }
            }
            //
            SimState.output_just_written = false;
            SimState.history_just_written = false;
            SimState.loads_just_written = false;
            if ((SimState.step % GlobalConfig.print_count) == 0) {
                // Print the current time-stepping status.
                auto writer = appender!string();
                if (!GlobalConfig.fixed_time_step) {
                    formattedWrite(writer, "Step=%7d t=%10.3e dt=%10.3e cfl=%.2f ",
                               SimState.step, SimState.time, SimState.dt_global, SimState.cfl_max);
                } else {
                    formattedWrite(writer, "Step=%7d t=%10.3e dt=%10.3e cfl=N/A ",
                               SimState.step, SimState.time, SimState.dt_global);
                }

                // For reporting wall-clock time, convert to seconds with precision of milliseconds.
                double wall_clock_elapsed = to!double((Clock.currTime()-SimState.wall_clock_start).total!"msecs"())/1000.0;
                double wall_clock_per_step = wall_clock_elapsed / SimState.step;
                double WCtFT = (GlobalConfig.max_time - SimState.time) / SimState.dt_global * wall_clock_per_step;
                double WCtMS = (GlobalConfig.max_step - SimState.step) * wall_clock_per_step;
                formattedWrite(writer, "WC=%.1f WCtFT=%.1f WCtMS=%.1f",
                               wall_clock_elapsed, WCtFT, WCtMS);
                if (GlobalConfig.verbosity_level >= 0 && GlobalConfig.is_master_task) {
                    writeln(writer.data);
                    stdout.flush();
                }
                version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }

                if (GlobalConfig.writeTransientResiduals) {
                    // We also compute the residual information and write to residuals file.
                    // These data can be used to monitor the progress of a steady-state calculation.
                    compute_mass_balance(mass_balance);
	                compute_L2_residual(L2_residual);
                    compute_Linf_residuals(Linf_residuals);
                    auto cqi = GlobalConfig.cqi;
                    version(mpi_parallel) {
                        // Reduce residual values across MPI tasks.
                        double my_local_value;
                        foreach (i; 0 .. cqi.n) {
                            my_local_value = Linf_residuals[i].re;
                            MPI_Allreduce(MPI_IN_PLACE, &my_local_value, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                            Linf_residuals[i] = to!number(my_local_value);
                        }
                        my_local_value = mass_balance.re;
                        MPI_Allreduce(MPI_IN_PLACE, &my_local_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                        mass_balance.re = my_local_value;
                        my_local_value = L2_residual.re;
                        MPI_Allreduce(MPI_IN_PLACE, &my_local_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                        L2_residual.re = my_local_value;
                    }
                    L2_residual = sqrt(L2_residual);
                    if (GlobalConfig.is_master_task) {
                        string txt = format("%7d %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e\n",
                                            SimState.step, SimState.time, wall_clock_elapsed,
                                            Linf_residuals[cqi.mass].re,
                                            Linf_residuals[cqi.xMom].re,
                                            Linf_residuals[cqi.yMom].re,
                                            (cqi.threeD) ? Linf_residuals[cqi.zMom].re : 0.0,
                                            Linf_residuals[cqi.totEnergy].re,
                                            fabs(L2_residual.re), fabs(mass_balance.re));
                        std.file.append(lmrCfg.transResidFile, txt);
                    }
                } // end if writeTransienResiduals
                version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }
                if (GlobalConfig.with_super_time_stepping) {
                    if (GlobalConfig.is_master_task) {
                        auto writer2 = appender!string();
                        formattedWrite(writer2, "SUPER-TIME-STEP: ");
                        formattedWrite(writer2, "S= %d dtHYPER= %10.6e dtPARAB= %10.6e ",
                                       SimState.s_RKL, SimState.dt_global, SimState.dt_global_parab);
                        writeln(writer2.data);
                        stdout.flush();
                    }
                } // end if with_super_time_stepping
                version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }
            } // end if (step...

            // [TODO] RJG, 2024-02-07
            // This bit of code will no longer work in manner below.
            // We'll think on an arrangement for auxiliary variables that works with
            // our new new I/O.
            /*
            // update any auxiliary values
            if (GlobalConfig.new_flow_format) {
                foreach (myblk; parallel(localFluidBlocksBySize,1)) {
                    myblk.update_aux(SimState.dt_global, SimState.time, SimState.step);
                }
            }
            */

            //
            // 4.0 (Occasionally) Write out an intermediate solution
            if ( SimState.step == GlobalConfig.write_flow_solution_at_step ) {
                writeSnapshotFiles_timemarching();
                GC.collect();
                GC.minimize();
            }
            if ((SimState.time >= SimState.t_plot) && !SimState.output_just_written) {
                writeSnapshotFiles_timemarching();
                if (GlobalConfig.udf_supervisor_file.length > 0) { call_UDF_at_write_to_file(); }
                SimState.output_just_written = true;
                SimState.t_plot = SimState.t_plot + GlobalConfig.dt_plot_schedule.get_value(SimState.time);
                GC.collect();
                GC.minimize();
            }
            // [TODO] RJG, 2024-02-07
            // This type of snapshot below is now an overloaded term.
            // What's more the use case for this might be vanishingly small.
            // This is another bit of code to consider if it's time to hit the super fund.
            //
            /*
            // 4.1 (Occasionally) Write out a snapshot
            if ((SimState.step % GlobalConfig.snapshotCount) == 0) {
                if (GlobalConfig.is_master_task) {
                    writeln("***");
                    writefln("*** Writing a snapshot: step= %4d, t= %6.3e", SimState.step, SimState.time);
                    writeln("***");
                }
                write_snapshot_files();
                GC.collect();
                GC.minimize();
            }
            */
            // 4.2 (Occasionally) Write out the cell history data and loads on boundary groups data
            if ((SimState.time >= SimState.t_history) && !SimState.history_just_written) {
                version(FSI) {
                    foreach (FEMModel; FEMModels) {
                        FEMModel.WriteFSIToHistory(SimState.time);
                    }
                }
                writeHistoryCellsToFiles(SimState.time);
                SimState.history_just_written = true;
                SimState.t_history = SimState.t_history + GlobalConfig.dt_history;
                GC.collect();
                GC.minimize();
            }

            if (GlobalConfig.write_loads &&
                ( ((SimState.time >= SimState.t_loads) && !SimState.loads_just_written) ||
                  SimState.step == GlobalConfig.write_loads_at_step )) {
                writeLoadsFiles_timemarching();
                SimState.loads_just_written = true;
                SimState.current_loads_tindx = SimState.current_loads_tindx + 1;
                SimState.t_loads = SimState.t_loads + GlobalConfig.dt_loads;
                GC.collect();
                GC.minimize();
            }

            /* RJG, 2024-02-12  commented out during initial testing.
            //
            // 4.3 Increment the DFT in each cell
            if ((SimState.step % GlobalConfig.DFT_step_interval == 0) && GlobalConfig.do_temporal_DFT) {
                foreach (blk; localFluidBlocks) {
                    blk.increment_DFT(SimState.step / GlobalConfig.DFT_step_interval - 1);
                }
            }
            */
            //
            // 5.0 Update the run-time loads calculation, if required
            if (GlobalConfig.compute_run_time_loads) {
                if ((SimState.step % GlobalConfig.run_time_loads_count) == 0) {
                    if (GlobalConfig.viscous) {
                        // Prep the viscous fluxes by recomputing.
                        // We'll wipe out everything else in the fluxes vector
                        // and recall the the viscous_flux calc so that we get
                        // only the viscous contribution.
                        foreach (blk; parallel(localFluidBlocksBySize,1)) {
                            blk.clear_fluxes_of_conserved_quantities();
                            if (blk.active) {
                                blk.average_turbulent_transprops_to_faces();
                                blk.viscous_flux();
                            }
                        }
                    }
                    computeRunTimeLoads();
                }
            }
            //
            // 6.0 Allow the user to do special actions at the end of a timestep..
            if (GlobalConfig.udf_supervisor_file.length > 0) { call_UDF_at_timestep_end(); }
            //
        } catch(Exception e) {
            writefln("Exception caught while trying to take step %d.", SimState.step);
            writeln("----- Begin exception message -----");
            writeln(e);
            writeln("----- End exception message -----");
            caughtException = true;
            finished_time_stepping = true;
        }
        //
        // Loop termination criteria:
        // (a) Exception caught.
        // (b) Reaching a maximum simulation time or target time.
        // (c) Reaching a maximum number of steps.
        // (d) Finding that the "halt_now" parameter has been set
        //     in the control-parameter file.
        //     This provides a semi-interactive way to terminate the
        //     simulation and save the data.
        // (e) Exceeding a maximum number of wall-clock seconds.
        //
        // Note that the max_time and max_step control parameters can also
        // be found in the control-parameter file (which may be edited
        // while the code is running).
        //
        if (SimState.time >= SimState.target_time) { finished_time_stepping = true; }
        if (SimState.step >= GlobalConfig.max_step) { finished_time_stepping = true; }
        if (GlobalConfig.halt_now == 1) { finished_time_stepping = true; }
        auto wall_clock_elapsed = (Clock.currTime() - SimState.wall_clock_start).total!"seconds"();
        if (SimState.maxWallClockSeconds > 0 && (wall_clock_elapsed > SimState.maxWallClockSeconds)) {
            finished_time_stepping = true;
        }
        version(mpi_parallel) {
            // If one task is finished time-stepping, all tasks have to finish.
            int localFinishFlag = to!int(finished_time_stepping);
            MPI_Allreduce(MPI_IN_PLACE, &localFinishFlag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            finished_time_stepping = to!bool(localFinishFlag);
        }
        if(finished_time_stepping && GlobalConfig.is_master_task) {
            // Make an announcement about why we are finishing time-stepping.
            write("STOP-REASON: ");
            if (caughtException) {
                writeln("An exception was caught during the time step.");
            }
            if (SimState.time >= SimState.target_time) {
                writefln("maximum-time");
                writefln("Reached target simulation time of %g seconds.", SimState.target_time);
            }
            if (SimState.step >= GlobalConfig.max_step) {
                writefln("maximum-steps");
                writefln("Reached maximum number of steps with step=%d.", SimState.step);
            }
            if (GlobalConfig.halt_now == 1) { writeln("Halt set in control file."); }
            if (SimState.maxWallClockSeconds > 0 && (wall_clock_elapsed > SimState.maxWallClockSeconds)) {
                writefln("maximum-wall-clock");
                writefln("Reached maximum wall-clock time with elapsed time %s.", to!string(wall_clock_elapsed));
            }
            writefln("FINAL-STEP: %d", SimState.step);
            writefln("FINAL-TIME: %g", SimState.time);
            stdout.flush();
        }
    } // end while !finished_time_stepping
    // [TODO] RJG, 2024-02-07
    // Disabled for a little bit while we get the new time-marching code going.
    //
    /*
    if (GlobalConfig.solve_electric_field && caughtException==false){
        if (GlobalConfig.is_master_task) writeln("Called field.solve_efield(): ...");
        eField.solve_efield(localFluidBlocks, GlobalConfig.is_master_task);
        eField.compute_electric_field_vector(localFluidBlocks);

        double current_in, current_out;
        eField.compute_boundary_current(localFluidBlocks, current_in, current_out);
        if (GlobalConfig.is_master_task) {
            writeln("Called field.compute_boundary_current() ...");
            writefln("    Current in:  %f (A/m)", current_in);
            writefln("    Current out: %f (A/m)", current_out);
        }
    }
    */
    if (GlobalConfig.verbosity_level > 0 && GlobalConfig.is_master_task) {
        writeln("Done integrate_in_time().");
        stdout.flush();
    }
    // We use a strange mix of exception handling and error flags
    // because we could not think of an easier way to deal with
    // the MPI issues of shutting down cleanly.
    if (caughtException) {
        return -1; // failed integration
    } else {
        return 0; // success
    }
}

/**
 * This function writes the loads files
 *
 * Authors: KAD
 * Date: 2024-03-15
 */
void writeLoadsFiles_timemarching()
{
    alias cfg = GlobalConfig;
    if (cfg.is_master_task) {
        writefln("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        writefln("+   Writing loads at step = %6d; t = %8.3e s    +", SimState.step, SimState.time);
        writefln("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
    }

    if (cfg.is_master_task) {
        init_current_loads_indx_dir(SimState.current_loads_tindx);
    }
    version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }
    wait_for_current_indx_dir(SimState.current_loads_tindx);
    writeLoadsToFile(SimState.time, SimState.current_loads_tindx);
    if (cfg.is_master_task) {
        update_loads_metadata_file(SimState.time, SimState.step, SimState.current_loads_tindx);
    }
}

/**
 * This function writes the flow field files, but also takes care of accompanying step <-> time mapping.
 *
 * Authors: RJG
 * Date: 2024-02-07
 */
void writeSnapshotFiles_timemarching()
{
    alias cfg = GlobalConfig;
    if (cfg.is_master_task) {
        writefln("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
        writefln("+   Writing snapshot at step = %6d; t = %8.3e s +", SimState.step, SimState.time);
        writefln("++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
    }

    // We're writing a snapshot, so time to bump the tindx.
    SimState.current_tindx = SimState.current_tindx + 1;

    // Write snapshot
    auto dirName = snapshotDirectory(SimState.current_tindx);
    if (cfg.is_master_task) {
        ensure_directory_is_present(dirName);
    }
    // Wait for master to complete building the directory for the next snapshot
    version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }

    foreach (blk; localFluidBlocksBySize) { // 2024-03-05 PJ remove parallel
        auto fileName = fluidFilename(SimState.current_tindx, blk.id);
        FVCell[] cells;
        cells.length = blk.cells.length;
        foreach (i, ref c; cells) c = blk.cells[i];
        fluidBlkIO.writeVariablesToFile(fileName, cells);
    }

    if (cfg.grid_motion != GridMotion.none) {
        foreach (blk; localFluidBlocksBySize) { // 2024-03-05 PJ remove parallel
            blk.sync_vertices_to_underlying_grid(0);
            auto gridName = gridFilename(SimState.current_tindx, blk.id);
            blk.write_underlying_grid(gridName);
        }
    }

    foreach (blk; localSolidBlocks) { // 2024-03-05 PJ remove parallel
        auto filename = solidFilename(SimState.current_tindx, blk.id);
        FVCell[] cells;
        cells.length = blk.cells.length;
        foreach (i, ref c; cells) c = blk.cells[i];
        solidBlkIO.writeVariablesToFile(filename, cells);
    }

    if (cfg.is_master_task) {
        addToTimesFile();
    }
}

/**
 * Add an entry in the times file.
 *
 * NOTE:
 * What happens if this is a restart from some earlier tindx?
 * Perhaps there is already an entry for this particular tindx (from a previous run).
 * Using YAML, we just stick the duplicate key at the end.
 * The parsers should hold onto the *last* duplicate key they find, which is what we want on restarted sims.
 * We might need to think about properly scrubbing the file during initialisation.
 *
 * Author: Rowan J. Gollan
 * Date: 2024-02-24
 */

void addToTimesFile()
{
    // Add entry in times file
    auto f = File(lmrCfg.timesFile, "a");
    string key = format(lmrCfg.snapshotIdxFmt, SimState.current_tindx);
    f.writefln("'%s':", key);
    f.writefln("   step: %d", SimState.step);
    f.writefln("   time: %.18e", SimState.time);
    f.writefln("   dt:   %.18e", SimState.dt_global);
    if (isNaN(SimState.cfl_max)) {
        // Then it has never been set.
        // This might occur on initialisation, or if we are running with fixed timestepping.
        // Put a value of -1.0 as the signal to outside world.
        f.writeln("   cfl:  -1.0");
    }
    else {
        f.writefln("   cfl:  %.18e", SimState.cfl_max);
    }
    double wall_clock_elapsed = to!double((Clock.currTime() - SimState.wall_clock_start).total!"msecs"())/1000.0;
    f.writefln("   wall-clock-elapsed: %.3f", wall_clock_elapsed);
    f.close();
}

/**
 * Prepare times file on restart.
 *
 * When we restart, we want to keep all times file entries up
 * to and including restart index. For those beyond,
 * we want to remove.
 *
 * Authors: RJG
 * Date: 2024-03-10
 */
void prepareTimesFileOnRestart(int restartIdx)
{
    Node timesData = dyaml.Loader.fromFile(lmrCfg.timesFile).load();
    string snap;
    Node snapData;
    auto f = File(lmrCfg.timesFile, "w");
    foreach (i; 0 .. restartIdx+1) {
        snap = format(lmrCfg.snapshotIdxFmt, i);
        snapData = timesData[snap];
        f.writefln("'%s':", snap);
        f.writefln("   step: %d", snapData["step"].as!int);
        f.writefln("   time: %.18e", snapData["time"].as!double);
        f.writefln("   dt:   %.18e", snapData["dt"].as!double);
        f.writefln("   cfl:  %.18e", snapData["cfl"].as!double);
        f.writefln("   wall-clock-elapsed: %.3f", snapData["wall-clock-elapsed"].as!double);
    }
    f.close();
    SimState.step = snapData["step"].as!int;
    SimState.time = snapData["time"].as!double;
    SimState.dt_global = snapData["dt"].as!double;
}

/**
 * This function is resposible for cleaning up (typically open files and recording state)
 * at the end of a transient simulation.
 *
 * Authors: PAJ and RJG
 * Date: 2024-02-12
 */

void finalizeSimulation_timemarching()
{
    if (GlobalConfig.verbosity_level > 0 && GlobalConfig.is_master_task) {
        writeln("Finalize the simulation.");
    }
    if (!SimState.output_just_written) {
        writeSnapshotFiles_timemarching();
        if (GlobalConfig.udf_supervisor_file.length > 0) { call_UDF_at_write_to_file(); }
    }

    /* RJG, 2024-02-12 commented out during initial testing.
     * Let's temporarily disable this auxiliary I/O while we get lmr5 going for transient mode.

    if (!SimState.history_just_written) {
        write_history_cells_to_files(SimState.time);
        version(FSI) {
            foreach (FEMModel; FEMModels) {
                FEMModel.WriteFSIToHistory(SimState.time);
            }
        }
    }
    */
    
    if (GlobalConfig.write_loads && !SimState.loads_just_written) {
        writeLoadsFiles_timemarching();
    }
    /*
    if (!GlobalConfig.new_flow_format && GlobalConfig.do_temporal_DFT) {
        write_DFT_files();
    }
    */

    foreach (blk; localFluidBlocks) {
        foreach (bc; blk.bc) { bc.finalize(); }
        blk.finalize();
    }
    GC.collect();
    GC.minimize();
    if (GlobalConfig.verbosity_level > 0  && GlobalConfig.is_master_task) {
        writeln("Step= ", SimState.step, " final-t= ", SimState.time);
    }
    if (GlobalConfig.is_master_task) {
        try {
            std.file.write(lmrCfg.progFile, "done\n");
        } catch (Exception e) {
            // do nothing
        }
    }
} // end finalize_simulation()

