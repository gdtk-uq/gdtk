/** simcore.d
 * Eilmer4 compressible-flow simulation code, core coordination functions.
 *
 * Author: Peter J. and Rowan G. 
 * First code: 2015-02-05
 */

module simcore;

import core.memory;
import std.math;
import std.stdio;
import std.file;
import std.conv;
import std.array;
import std.format;
import std.string;
import std.algorithm;
import std.typecons;
import std.datetime;
import std.parallelism;
import std.json;
import nm.complex;
import nm.number;

import util.lua;
import util.lua_service;
import lua_helper;
import fileutil;
import geom;
import gas;
import fvcore;
import globalconfig;
import globaldata;
import flowstate;
import fluidblock;
import sfluidblock;
import ufluidblock;
import ssolidblock;
import solidprops;
import solidfvinterface;
import solid_loose_coupling_update;
import bc;
import user_defined_source_terms;
import solid_udf_source_terms;
import shock_fitting_moving_grid;
import grid_motion;
import history;
import loads;
import conservedquantities;
import special_block_init;
version (opencl_gpu_chem) {
    import opencl_gpu_chem;
}
version (cuda_gpu_chem) {
    import cuda_gpu_chem;       
}
version(mpi_parallel) {
    import mpi;
    import mpi.util;
}

// State data for simulation.
// Needs to be seen by all of the coordination functions.
final class SimState {
    shared static double time;  // present simulation time, tracked by code
    shared static int step;
    shared static double dt_global;     // simulation time step determined by code
    shared static double dt_allow;      // allowable global time step determined by code
    // for STS
    shared static double dt_global_parab;
    shared static double dt_allow_parab; 
    //
    shared static double cfl_max;      // current max cfl determined by code
    shared static double dt_override = 0.0;  // A positive value will override a larger computed time step.
    shared static double target_time;  // simulate_in_time will work toward this value

    // We want to write sets of output files periodically.
    // The following periods set the cadence for output.
    shared static double t_plot;
    shared static double t_history;
    shared static double t_loads;
    // Once we write some data to files, we don't want to write another set of files
    // until we have done some more stepping.  The following flags help us remember
    // the state of the solution output.
    shared static bool output_just_written = true;
    shared static bool history_just_written = true;
    shared static bool loads_just_written = true;
    // We connect the sets of files to the simulation time at which they were written
    // with an index that gets incremented each time we write a set of files.
    shared static int current_tindx;
    shared static int current_loads_tindx;

    // For working out how long the simulation has been running.
    static SysTime wall_clock_start;
    static int maxWallClockSeconds;
} // end class SimState

// To avoid race conditions, there are a couple of locations where
// each block will put its result into the following arrays,
// then we will reduce across the arrays.
shared static double[] local_dt_allow;
shared static double[] local_dt_allow_parab;
shared static double[] local_cfl_max;
shared static int[] local_invalid_cell_count;

// The shared double[] flavour of GlobalConfig.userPad can give trouble,
// so we need a normal array for the MPI task to work with.
double[] userPad_copy;

//----------------------------------------------------------------------------

version(mpi_parallel) {
    void MPI_Sync_tasks()
    // This function is essentially an MPI_Barrier with check for liveness.
    // The MPI library does not have any timeout on MPI_Barrier so we first
    // use the nonblocking communication to test that our neighbour is still 
    // responding within a reasonable time.
    // If a neighbour does not respond within the timeout period, call a halt.
    {
        int my_rank = GlobalConfig.mpi_rank_for_local_task;
        //
        // Look to neighbour on left for incoming message.
        int left_rank = (my_rank == 0) ? GlobalConfig.mpi_size-1 : my_rank-1;
        int incoming_tag = left_rank;
        MPI_Request receive_request;
        int[1] incoming_buffer;
        MPI_Irecv(incoming_buffer.ptr, 1, MPI_INT, left_rank, incoming_tag,
                  MPI_COMM_WORLD, &receive_request);
        //
        // Send message to neighbour on right.
        MPI_Request send_request;
        int outgoing_tag = my_rank;
        int right_rank = (my_rank == GlobalConfig.mpi_size-1) ? 0 : my_rank+1;
        int[1] outgoing_buffer;
        outgoing_buffer[0] = my_rank;
        MPI_Isend(outgoing_buffer.ptr, 1, MPI_INT, right_rank, outgoing_tag,
                  MPI_COMM_WORLD, &send_request);
        //
        // Wait for message to be sent.
        long timeout_msecs = 10000; // Surely 10 seconds will be enough.
        SysTime startTime = Clock.currTime();
        MPI_Status send_status;
        int send_flag = 0;
        while (!send_flag) {
            int ierr = MPI_Test(&send_request, &send_flag, &send_status);
            long elapsedTime_msecs = (Clock.currTime() - startTime).total!"msecs"();
            if (elapsedTime_msecs > timeout_msecs) {
                // We do not expect our job to recover gracefully from this point.
                writeln("MPI_Sync_tasks time-out waiting to send. my_rank=%d, right_rank=%d",
                        my_rank, right_rank);
                MPI_Abort(MPI_COMM_WORLD, 3);
            }
        }
        //
        // Wait for expected message be received.
        startTime = Clock.currTime();
        MPI_Status receive_status;
        int receive_flag = 0;
        while (!receive_flag) {
            int ierr = MPI_Test(&receive_request, &receive_flag, &receive_status);
            long elapsedTime_msecs = (Clock.currTime() - startTime).total!"msecs"();
            if (elapsedTime_msecs > timeout_msecs) {
                // We do not expect our job to recover gracefully from this point.
                writeln("MPI_Sync_tasks time-out waiting for receive. my_rank=%d, left_rank=%d",
                        my_rank, left_rank);
                MPI_Abort(MPI_COMM_WORLD, 3);
            }
        }
        //
        // At this point, we know that this rank and its left-neighbour are alive,
        // so we wait for everyone to come to the same conclusion.
        MPI_Barrier(MPI_COMM_WORLD);
    } // end MPI_Sync_tasks()
}

//----------------------------------------------------------------------------

void init_simulation(int tindx, int nextLoadsIndx,
                     int maxCPUs, int threadsPerMPITask, int maxWallClock)
{
    if (GlobalConfig.verbosity_level > 0 && GlobalConfig.is_master_task) {
        writeln("Begin init_simulation()...");
    }
    SimState.maxWallClockSeconds = maxWallClock;
    SimState.wall_clock_start = Clock.currTime();
    read_config_file();  // most of the configuration is in here
    read_control_file(); // some of the configuration is in here
    //
    version(enable_fp_exceptions) {
        FloatingPointControl fpctrl;
        // Enable hardware exceptions for division by zero, overflow to infinity,
        // invalid operations, and uninitialized floating-point variables.
        // Copied from https://dlang.org/library/std/math/floating_point_control.html
        fpctrl.enableExceptions(FloatingPointControl.severeExceptions);
    }
    //
    if (GlobalConfig.grid_format == "rawbinary") { GlobalConfig.gridFileExt = "bin"; }
    if (GlobalConfig.flow_format == "rawbinary") { GlobalConfig.flowFileExt = "bin"; }
    setupIndicesForConservedQuantities(); 
    SimState.current_tindx = tindx;
    SimState.current_loads_tindx = nextLoadsIndx;
    if (SimState.current_loads_tindx == -1) {
        // This is used to indicate that we should look at the old -loads.times file
        // to find the index where the previous simulation stopped.
        string fname = "loads/" ~ GlobalConfig.base_file_name ~ "-loads.times";
        if (exists(fname)) {
            auto finalLine = readText(fname).splitLines()[$-1];
            if (finalLine[0] == '#') {
                // looks like we found a single comment line.
                SimState.current_loads_tindx = 0;
            } else {
                // assume we have a valid line to work with
                SimState.current_loads_tindx = to!int(finalLine.split[0]) + 1;
            }
        } else {
            SimState.current_loads_tindx = 0;
        }
    }
    auto job_name = GlobalConfig.base_file_name;
    if (GlobalConfig.nFluidBlocks == 0 && GlobalConfig.is_master_task) {
        throw new FlowSolverException("No FluidBlocks; no point in continuing to initialize simulation.");
    }
    version(mpi_parallel) {
        // Assign particular fluid blocks to this MPI task and keep a record
        // of the MPI rank for all blocks.
        int my_rank = GlobalConfig.mpi_rank_for_local_task;
        GlobalConfig.mpi_rank_for_block.length = GlobalConfig.nFluidBlocks;
        auto lines = readText("config/" ~ job_name ~ ".mpimap").splitLines();
        foreach (line; lines) {
            auto content = line.strip();
            if (content.startsWith("#")) continue; // Skip comment
            auto tokens = content.split();
            int blkid = to!int(tokens[0]);
            int taskid = to!int(tokens[1]);
            if (taskid >= GlobalConfig.mpi_size && GlobalConfig.is_master_task) {
                writefln("Number of MPI tasks (%d) is insufficient for "~
                         "taskid=%d that is associated with blockid=%d. Quitting.",
                         GlobalConfig.mpi_size, taskid, blkid);
                MPI_Abort(MPI_COMM_WORLD, 2);
            }
            GlobalConfig.mpi_rank_for_block[blkid] = taskid;
            if (taskid == my_rank) { localFluidBlocks ~= globalFluidBlocks[blkid]; }
        }
        MPI_Sync_tasks();
        if (localFluidBlocks.length == 0) {
            writefln("MPI-task with rank %d has no FluidBlocks. Quitting.", my_rank);
            MPI_Abort(MPI_COMM_WORLD, 2);
        }
    } else {
        // There is only one process and it deals with all blocks.
        foreach (blk; globalFluidBlocks) { localFluidBlocks ~= blk; }
    }
    foreach (blk; localFluidBlocks) { GlobalConfig.localBlockIds ~= blk.id; }
    //
    // Local blocks may be handled with thread-parallelism.
    auto nBlocksInThreadParallel = max(localFluidBlocks.length, GlobalConfig.nSolidBlocks);
    // There is no need to have more task threads than blocks local to the process.
    int extraThreadsInPool;
    version(mpi_parallel) {
        extraThreadsInPool = min(threadsPerMPITask-1, nBlocksInThreadParallel-1);
    } else {
        extraThreadsInPool = min(maxCPUs-1, nBlocksInThreadParallel-1);
    }
    defaultPoolThreads(extraThreadsInPool); // total = main thread + extra-threads-in-Pool
    if (GlobalConfig.verbosity_level > 0) {
        version(mpi_parallel) {
            writeln("MPI-task with rank ", my_rank, " running with ", extraThreadsInPool+1, " threads.");
            debug {
                foreach (blk; localFluidBlocks) { writeln("rank=", my_rank, " blk.id=", blk.id); }
            }
        } else {
            writeln("Single process running with ", extraThreadsInPool+1, " threads.");
            // Remember the +1 for the main thread.
        }
    }
    // At this point, note that we initialize the grid and flow arrays for blocks
    // that are in the current MPI-task or process, only.
    foreach (myblk; parallel(localFluidBlocks,1)) {
        if (GlobalConfig.grid_motion != GridMotion.none) {
            myblk.init_grid_and_flow_arrays(make_file_name!"grid"(job_name, myblk.id, SimState.current_tindx,
                                                                  GlobalConfig.gridFileExt)); 
        } else {
            // Assume there is only a single, static grid stored at tindx=0
            myblk.init_grid_and_flow_arrays(make_file_name!"grid"(job_name, myblk.id, 0, GlobalConfig.gridFileExt)); 
        }
        myblk.compute_primary_cell_geometric_data(0);
    }
    // Note that the global id is across all processes, not just the local collection of blocks.
    foreach (i, myblk; globalFluidBlocks) {
        myblk.globalCellIdStart = (i == 0) ? 0 : globalFluidBlocks[i-1].globalCellIdStart + globalFluidBlocks[i-1].ncells_expected;
    }
    shared double[] time_array;
    time_array.length = localFluidBlocks.length;
    foreach (i, myblk; parallel(localFluidBlocks,1)) {
        myblk.identify_reaction_zones(0);
        myblk.identify_turbulent_zones(0);
        myblk.identify_suppress_reconstruction_zones();
        time_array[i] = myblk.read_solution(make_file_name!"flow"(job_name, myblk.id, SimState.current_tindx,
                                                                  GlobalConfig.flowFileExt), false);
        if (myblk.myConfig.verbosity_level >= 2) { writefln("Cold start cells in block %d", myblk.id); }
        foreach (iface; myblk.faces) { iface.gvel.clear(); }
        foreach (cell; myblk.cells) {
            cell.encode_conserved(0, 0, myblk.omegaz);
            // Even though the following call appears redundant at this point,
            // fills in some gas properties such as Prandtl number that is
            // needed for both the cfd_check and the BaldwinLomax turbulence model.
            if (0 != cell.decode_conserved(0, 0, myblk.omegaz)) {
                throw new FlowSolverException("Bad cell decode_conserved while initializing.");
            }
        }
        myblk.set_cell_dt_chem(-1.0);
    }
    SimState.time = time_array[0]; // Pick one; they should all be the same.
    //
    version(mpi_parallel) { MPI_Sync_tasks(); }
    //
    // Now that the cells for all gas blocks have been initialized,
    // we can sift through the boundary condition effects and
    // set up the ghost-cell mapping for the appropriate boundaries.
    if (GlobalConfig.verbosity_level >= 2) { writeln("Prepare exchange of boundary information."); }
    // Serial loops because the cell-mapping function searches across
    // all blocks local to the process.
    // Also, there are several loops because the MPI communication,
    // if there is any, needs to be done in phases of posting of non-blocking reads,
    // followed by all of the sends and then waiting for all requests to be filled.
    //
    foreach (myblk; localFluidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellMappedCellCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping(); }
            }
        }
    }
    foreach (myblk; localFluidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellFullFaceCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping_phase0(); }
            }
        }
    }
    foreach (myblk; localFluidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellFullFaceCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping_phase1(); }
            }
        }
    }
    foreach (myblk; localFluidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellFullFaceCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping_phase2(); }
            }
        }
    }
    foreach (myblk; localFluidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy)gce;
                if (mygce1) { mygce1.exchange_geometry_phase0(); }
                auto mygce2 = cast(GhostCellFullFaceCopy)gce;
                if (mygce2) { mygce2.exchange_geometry_phase0(); }
            }
        }
    }
    foreach (myblk; localFluidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy)gce;
                if (mygce1) { mygce1.exchange_geometry_phase1(); }
                auto mygce2 = cast(GhostCellFullFaceCopy)gce;
                if (mygce2) { mygce2.exchange_geometry_phase1(); }
            }
        }
    }
    foreach (myblk; localFluidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy)gce;
                if (mygce1) { mygce1.exchange_geometry_phase2(); }
                auto mygce2 = cast(GhostCellFullFaceCopy)gce;
                if (mygce2) { mygce2.exchange_geometry_phase2(); }
            }
        }
    }
    //
    // Now that we know the ghost-cell locations, we can set up the least-squares subproblems for
    // 1. reconstruction prior to convective flux calculation for the unstructured-grid blocks
    // 2. calculation of flow gradients for the viscous fluxes with least-squares gradients.
    foreach (myblk; localFluidBlocks) { myblk.compute_least_squares_setup(0); }
    //
    version (gpu_chem) {
        initGPUChem();
    }
    //
    foreach (ref mySolidBlk; solidBlocks) {
        mySolidBlk.assembleArrays();
        mySolidBlk.bindFacesAndVerticesToCells();
        writeln("mySolidBlk= ", mySolidBlk);
        mySolidBlk.readGrid(make_file_name!"solid-grid"(job_name, mySolidBlk.id, 0, "gz")); // tindx==0 fixed grid
        mySolidBlk.readSolution(make_file_name!"solid"(job_name, mySolidBlk.id, tindx, "gz"));
        mySolidBlk.computePrimaryCellGeometricData();
        mySolidBlk.assignVtxLocationsForDerivCalc();
    }
    if (solidBlocks.length > 0) {
        initPropertiesAtSolidInterfaces(solidBlocks);
    }
    if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.loose) {
        initSolidLooseCouplingUpdate();
    }
    //
    // All cells are in place, so now we can initialise any history cell files.
    if (GlobalConfig.is_master_task) { ensure_directory_is_present(histDir); }
    version(mpi_parallel) { MPI_Sync_tasks(); }
    init_history_cell_files();
    //
    // create the loads directory, maybe
    if (GlobalConfig.write_loads && (SimState.current_loads_tindx == 0)) {
        if (GlobalConfig.is_master_task) { ensure_directory_is_present("loads"); }
        version(mpi_parallel) { MPI_Sync_tasks(); }
        init_loads_times_file();
    }
    version(mpi_parallel) { MPI_Sync_tasks(); }

    // For the shock fitting grid motion, we need to assign radial positions for all vertices
    if (GlobalConfig.grid_motion == GlobalConfig.grid_motion.shock_fitting) {
        foreach (myblk; localFluidBlocks) {
            version(mpi_parallel) {
                auto sblk = cast(SFluidBlock) myblk;
                assign_radial_dist_mpi(sblk);
            }
            else {
                if (myblk.bc[Face.west].type == "inflow_shock_fitting") {
                    auto sblk = cast(SFluidBlock) myblk;
                    assign_radial_dist(sblk);
                }
            }
        }
    }
    // Finally when both gas AND solid domains are setup..
    // Look for a solid-adjacent bc, if there is one,
    // then we can set up the cells and interfaces that
    // internal to the bc. They are only known after this point.
    if (GlobalConfig.apply_bcs_in_parallel) {
        foreach (myblk; parallel(localFluidBlocks,1)) {
            foreach (bc; myblk.bc) {
                foreach (bfe; bc.postDiffFluxAction) {
                    auto mybfe = cast(BFE_EnergyFluxFromAdjacentSolid)bfe;
                    if (mybfe) {
                        if (GlobalConfig.in_mpi_context) { throw new Error("[TODO] not available in MPI context."); }
                        auto adjSolidBC = to!BFE_EnergyFluxFromAdjacentSolid(mybfe);
                        adjSolidBC.initGasCellsAndIFaces();
                        adjSolidBC.initSolidCellsAndIFaces();
                    }
                }
            }
        }
    } else {
        foreach (myblk; localFluidBlocks) {
            foreach (bc; myblk.bc) {
                foreach (bfe; bc.postDiffFluxAction) {
                    auto mybfe = cast(BFE_EnergyFluxFromAdjacentSolid)bfe;
                    if (mybfe) {
                        if (GlobalConfig.in_mpi_context) { throw new Error("[TODO] not available in MPI context."); }
                        auto adjSolidBC = to!BFE_EnergyFluxFromAdjacentSolid(mybfe);
                        adjSolidBC.initGasCellsAndIFaces();
                        adjSolidBC.initSolidCellsAndIFaces();
                    }
                }
            }
        }
    }
    // For the MLP limiter (on unstructured grids only), we need access to the
    // gradients stored in the cloud of cells surrounding a vertex.
    if ((GlobalConfig.interpolation_order > 1) &&
        (GlobalConfig.unstructured_limiter == UnstructuredLimiter.mlp)) {
        foreach (myblk; localFluidBlocks) {
            auto ublock = cast(UFluidBlock) myblk;
            if (ublock) { ublock.build_cloud_of_cell_references_at_each_vertex(); }
        }
    }
    //
    // We can apply a special initialisation to the flow field, if requested.
    // This will take viscous boundary conditions and diffuse them into the
    // nearby domain.
    if (GlobalConfig.diffuseWallBCsOnInit) {
        writeln("Applying special initialisation to blocks: wall BCs being diffused into domain.");
        writefln("%d passes of the near-wall flow averaging operation will be performed.", GlobalConfig.nInitPasses);
        foreach (myblk; parallel(localFluidBlocks,1)) {
            diffuseWallBCsIntoBlock(myblk, GlobalConfig.nInitPasses, GlobalConfig.initTWall);
        }
    }
    //
    // We conditionally sort the local blocks, based on numbers of cells,
    // in an attempt to balance the load for shared-memory parallel runs.
    localFluidBlocksBySize.length = 0;
    if (GlobalConfig.block_marching || localFluidBlocks.length < 2) {
        // Keep the original block order, else we stuff up block-marching.
        // No point in sorting if there is only one local block.
        foreach (blk; localFluidBlocks) { localFluidBlocksBySize ~= blk; }
    } else {
        // We simply use the cell count as an estimate of load.
        Tuple!(int,int)[] blockLoads;
        blockLoads.length = localFluidBlocks.length;
        foreach (iblk, blk; localFluidBlocks) {
            // Here 'iblk' is equal to the local-to-process id of the block.
            blockLoads[iblk] = tuple(to!int(iblk), to!int(localFluidBlocks[iblk].cells.length)); 
        }
        sort!("a[1] > b[1]")(blockLoads);
        foreach (blkpair; blockLoads) {
            localFluidBlocksBySize ~= localFluidBlocks[blkpair[0]]; // [0] holds the block id
        }
    }
    //
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
    //
    // If we are using Lua supervisory script, fill in some more global
    // information for the interpreter.
    if (GlobalConfig.udf_supervisor_file.length > 0) {
        auto L = GlobalConfig.master_lua_State;
        lua_pushboolean(L, GlobalConfig.in_mpi_context);
        lua_setglobal(L, "in_mpi_context");
        lua_pushnumber(L, GlobalConfig.mpi_size);
        lua_setglobal(L, "mpi_size");
        lua_pushnumber(L, GlobalConfig.mpi_rank_for_local_task);
        lua_setglobal(L, "mpi_rank_for_local_task");
        lua_pushboolean(L, GlobalConfig.is_master_task);
        lua_setglobal(L, "is_master_task");
        push_array_to_Lua(L, GlobalConfig.localBlockIds, "localBlockIds");
        if (GlobalConfig.user_pad_length > 0) {
            push_array_to_Lua(L, GlobalConfig.userPad, "userPad");
        }
    }
    if (GlobalConfig.user_pad_length > 0) {
        // At this point, userPad has been initialized with values
        // from the job.config file, even if there are no supervisory functions.
        // All MPI tasks should see the same data and all interpreters
        // associated with blocks should get a copy of the data.
        broadcast_master_userPad();
        copy_userPad_into_block_interpreters();
    }
    //
    // Configure the run-time loads if required
    if (GlobalConfig.compute_run_time_loads) {
        string fileName = "config/" ~ GlobalConfig.base_file_name ~ ".config";
        string content = readText(fileName);
        JSONValue jsonData = parseJSON!string(content);
        initRunTimeLoads(jsonData["run_time_loads"]);
    }
    //
    // We want to get the corner-coordinates for each block copied
    // in case the Lua functions try to use them.
    synchronize_corner_coords_for_all_blocks();
    //
    // Keep our memory foot-print small.
    GC.collect();
    GC.minimize();
    //
    version(mpi_parallel) { MPI_Sync_tasks(); }
    if (GlobalConfig.verbosity_level > 0) {
        auto myStats = GC.stats();
        auto heapUsed = to!double(myStats.usedSize)/(2^^20);
        auto heapFree = to!double(myStats.freeSize)/(2^^20); 
        writefln("Heap memory used for task %d: %.2f  free: %.2f  total: %.1f MB",
                 GlobalConfig.mpi_rank_for_local_task, heapUsed, heapFree, heapUsed+heapFree);
        stdout.flush();
    }
    version(mpi_parallel) { MPI_Sync_tasks(); }
    if (GlobalConfig.verbosity_level > 0 && GlobalConfig.is_master_task) {
        // For reporting wall-clock time, convert to seconds with precision of milliseconds.
        double wall_clock_elapsed = to!double((Clock.currTime() - SimState.wall_clock_start).total!"msecs"())/1000.0;
        writefln("Done init_simulation() at wall-clock(WC)= %.1f sec", wall_clock_elapsed);
        stdout.flush();
    }
    return;
} // end init_simulation()

void write_solution_files()
{
    if (GlobalConfig.verbosity_level > 0 && GlobalConfig.is_master_task) {
        writeln("Write flow solution.");
        stdout.flush();
    }
    version(mpi_parallel) { MPI_Sync_tasks(); }
    SimState.current_tindx = SimState.current_tindx + 1;
    ensure_directory_is_present(make_path_name!"flow"(SimState.current_tindx));
    auto job_name = GlobalConfig.base_file_name;
    foreach (myblk; parallel(localFluidBlocksBySize,1)) {
        auto file_name = make_file_name!"flow"(job_name, myblk.id, SimState.current_tindx, GlobalConfig.flowFileExt);
        myblk.write_solution(file_name, SimState.time);
    }
    ensure_directory_is_present(make_path_name!"solid"(SimState.current_tindx));
    foreach (ref mySolidBlk; solidBlocks) {
        auto fileName = make_file_name!"solid"(job_name, mySolidBlk.id, SimState.current_tindx, "gz");
        mySolidBlk.writeSolution(fileName, SimState.time);
    }
    if (GlobalConfig.grid_motion != GridMotion.none) {
        ensure_directory_is_present(make_path_name!"grid"(SimState.current_tindx));
        if (GlobalConfig.verbosity_level > 0 && GlobalConfig.is_master_task) { writeln("Write grid"); }
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            blk.sync_vertices_to_underlying_grid(0);
            auto fileName = make_file_name!"grid"(job_name, blk.id, SimState.current_tindx, GlobalConfig.gridFileExt);
            blk.write_underlying_grid(fileName);
        }
    }
    // Update times file, connecting the tindx value to SimState.time.
    if (GlobalConfig.is_master_task) {
        auto writer = appender!string();
        formattedWrite(writer, "%04d %.18e %.18e\n", SimState.current_tindx, SimState.time, SimState.dt_global);
        append("config/" ~ GlobalConfig.base_file_name ~ ".times", writer.data);
    }
} // end write_solution_files()

void march_over_blocks()
{
    if (GlobalConfig.in_mpi_context) {
        throw new FlowSolverException("March over blocks is not available in MPI-parallel calculations.");
    }
    if (GlobalConfig.verbosity_level > 0) { writeln("March over blocks."); }
    // Organise the blocks into a regular array.
    int nib = GlobalConfig.nib;
    int njb = GlobalConfig.njb;
    int nkb = GlobalConfig.nkb;
    if (nib*njb*nkb != GlobalConfig.nFluidBlocks) {
        string errMsg = text("march_over_blocks(): inconsistent numbers of blocks\n",
                             "    nFluidBlocks=", GlobalConfig.nFluidBlocks,
                             " nib=", nib, " njb=", njb, " nkb=", nkb);
        throw new FlowSolverException(errMsg);
    }
    if (nkb != 1 && GlobalConfig.dimensions == 2) {
        string errMsg = text("march_over_blocks(): for 2D flow, expected nkb=1\n",
                             "    nkb=", nkb);
        throw new FlowSolverException(errMsg);
    }
    if (nib < 2) {
        string errMsg = text("march_over_blocks(): expected nib>=2\n",
                             "    nib=", nib);
        throw new FlowSolverException(errMsg);
    }
    FluidBlock[][][] gasBlockArray;
    gasBlockArray.length = nib;
    foreach (i; 0 .. nib) {
        gasBlockArray[i].length = njb;
        foreach (j; 0 .. njb) {
            gasBlockArray[i][j].length = nkb;
            foreach (k; 0 .. nkb) {
                int gid = k + nkb*(j + njb*i);
                gasBlockArray[i][j][k] = globalFluidBlocks[gid];
            }
        }
    }
    // Keep the first two slices active but deactivate the rest.
    foreach (i; 2 .. nib) {
        foreach (j; 0 .. njb) {
            foreach (k; 0 .. nkb) {
                gasBlockArray[i][j][k].active = false;
            }
        }
    }
    if (nib > 2) {
        // For the slice of blocks that are upstream of the inactive blocks,
        // set the downstream boundary condition to simple outflow,
        // saving the original boundary conditions so that they can be restored.
        foreach (j; 0 .. njb) {
            foreach (k; 0 .. nkb) {
                gasBlockArray[1][j][k].bc[Face.east].pushExtrapolateCopyAction();
            }
        }
    }
    //
    // Share the max_time for the overall calculation into smaller portions,
    // one for each block-slice.
    double time_slice = GlobalConfig.max_time / (nib - 1);
    // At most, we want to write out a set of solution files, only at the end
    // of integrating in time for each pair of block-slices.
    // To be sure that we don't write out lots of solutions,
    // we force dt_plot to be large enough.
    GlobalConfig.dt_plot = max(GlobalConfig.dt_plot, time_slice+1.0);
    // Let's start integrating for the first pair of block slices.
    if (integrate_in_time(SimState.time+time_slice) != 0) {
        throw new FlowSolverException("Integration failed for first pair of block slices.");
    }
    // Now, move along one block in i-direction at a time and do the rest.
    foreach (i; 2 .. nib) {
        // Deactivate the old-upstream slice, activate the new slice of blocks.
        // Restore boundary connections between currently-active upstream
        // and downstream block slices and
        // (so long as we are not at the last slice of blocks)
        // set downstream boundary condition for newly active slice to a simple outflow.
        foreach (j; 0 .. njb) {
            foreach (k; 0 .. nkb) {
                gasBlockArray[i-2][j][k].active = false;
                gasBlockArray[i-1][j][k].bc[Face.east].restoreOriginalActions();
                gasBlockArray[i][j][k].active = true;
                if (i < (nib-1)) { gasBlockArray[i][j][k].bc[Face.east].pushExtrapolateCopyAction(); }
            }
        }
        if (GlobalConfig.propagate_inflow_data) {
            exchange_ghost_cell_boundary_data(SimState.time, 0, 0);
            foreach (j; 0 .. njb) {
                foreach (k; 0 .. nkb) {
                    auto blk = gasBlockArray[i][j][k]; // our newly active block
                    // Get upstream flow data into ghost cells
                    blk.applyPreReconAction(SimState.time, 0, 0);
                    // and propagate it across the domain.
                    blk.propagate_inflow_data_west_to_east();
                }
            }
        }
        if (GlobalConfig.verbosity_level > 0) { writeln("march over blocks i=", i); }
        if (integrate_in_time(SimState.time+time_slice) != 0) {
            string msg = format("Integration failed for i=%d pair of block slices.", i);
            throw new FlowSolverException(msg);
        }
        if (GlobalConfig.save_intermediate_results) { write_solution_files(); }
    }
} // end march_over_blocks()

int integrate_in_time(double target_time_as_requested)
{
    version(enable_fp_exceptions) {
        FloatingPointControl fpctrl;
        // Enable hardware exceptions for division by zero, overflow to infinity,
        // invalid operations, and uninitialized floating-point variables.
        // Copied from https://dlang.org/library/std/math/floating_point_control.html
        fpctrl.enableExceptions(FloatingPointControl.severeExceptions);
    }
    //
    number mass_balance = to!number(0.0);
    number L2_residual = to!number(0.0);
    ConservedQuantities Linf_residuals = new ConservedQuantities(GlobalConfig.gmodel_master.n_species,
                                                                 GlobalConfig.gmodel_master.n_modes);
    version(mpi_parallel) { MPI_Sync_tasks(); }
    if (GlobalConfig.verbosity_level > 0 && GlobalConfig.is_master_task) {
        writeln("Integrate in time.");
        stdout.flush();
    }
    SimState.target_time = (GlobalConfig.block_marching) ? target_time_as_requested : GlobalConfig.max_time;
    // The next time for output...
    SimState.t_plot = SimState.time + GlobalConfig.dt_plot;
    SimState.t_history = SimState.time + GlobalConfig.dt_history;
    SimState.t_loads = SimState.time + GlobalConfig.dt_loads;
    // Overall iteration count.
    SimState.step = 0;
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
    local_dt_allow.length = localFluidBlocks.length; // prepare array for use later
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
    while ( !finished_time_stepping ) {
        try {
            // 0.0 Run-time configuration may change, a halt may be called, etc.
            check_run_time_configuration(target_time_as_requested);
            if (GlobalConfig.grid_motion != GridMotion.none) { synchronize_corner_coords_for_all_blocks(); }
            //
            // 1.0 Maintain a stable time step size, and other maintenance, as required.
            // The user may also have some start-of-time-step maintenance to do
            // via their Lua script file.  Let them have first go.
            if (GlobalConfig.udf_supervisor_file.length > 0) { call_UDF_at_timestep_start(); }            
            if (!GlobalConfig.fixed_time_step) { determine_time_step_size(); }
            if (GlobalConfig.divergence_cleaning) { update_ch_for_divergence_cleaning(); }
            // If using k-omega, we need to set mu_t and k_t BEFORE we call convective_update
            // because the convective update is where the interface values of mu_t and k_t are set.
            // only needs to be done on initial step, subsequent steps take care of setting these values 
            if ((SimState.step == 0) && (GlobalConfig.turbulence_model == TurbulenceModel.k_omega)) {
                k_omega_set_mu_and_k();
            }
            //
            // 2.0 Attempt a time step.
            // 2.1 Chemistry 1/2 step (if appropriate). 
            if(GlobalConfig.strangSplitting == StrangSplittingMode.half_R_full_T_half_R && GlobalConfig.with_local_time_stepping)
                assert(0, "Oops, strangSplitting.half_R_full_T_half_R and LTS aren't currently compatible");
            if (GlobalConfig.reacting && 
                (GlobalConfig.strangSplitting == StrangSplittingMode.half_R_full_T_half_R) &&
                (SimState.time > GlobalConfig.reaction_time_delay)) {
                chemistry_step(0.5*SimState.dt_global);
            }
            // 2.2 Update the convective terms.
            if (GlobalConfig.grid_motion == GridMotion.none) {
                if(GlobalConfig.with_super_time_stepping) { sts_gasdynamic_explicit_increment_with_fixed_grid(); }
		else { gasdynamic_explicit_increment_with_fixed_grid(); }
            } else {
                // Moving Grid - perform gas update for moving grid
                // [TODO] PJ 2018-01-20 Up to here with thinking about MPI parallel.
                set_grid_velocities();
                gasdynamic_explicit_increment_with_moving_grid();
                recalculate_all_geometry();
            }
            // 2.3 Solid domain update (if loosely coupled)
            // If tight coupling, then this has already been performed
            // in the gasdynamic_explicit_increment().
            if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.loose) {
                // Call Nigel's update function here.
                solid_domains_backward_euler_update(SimState.time, SimState.dt_global);
            }
            // 2.4 Chemistry step or 1/2 step (if appropriate). 
            if ( GlobalConfig.reacting && (SimState.time > GlobalConfig.reaction_time_delay)) {
                double mydt = (GlobalConfig.strangSplitting == StrangSplittingMode.full_T_full_R) ?
                    SimState.dt_global : 0.5*SimState.dt_global;
                chemistry_step(mydt);
            }
            //
            // 3.0 Update the time record and (occasionally) print status.
            SimState.step = SimState.step + 1;
            SimState.output_just_written = false;
            SimState.history_just_written = false;
            SimState.loads_just_written = false;
            if ((SimState.step % GlobalConfig.print_count) == 0) {
                // Print the current time-stepping status.
                auto writer = appender!string();
                formattedWrite(writer, "Step=%7d t=%10.3e dt=%10.3e cfl=%.2f ",
                               SimState.step, SimState.time, SimState.dt_global, SimState.cfl_max);
                // For reporting wall-clock time, convert to seconds with precision of milliseconds.
                double wall_clock_elapsed = to!double((Clock.currTime()-SimState.wall_clock_start).total!"msecs"())/1000.0;
                double wall_clock_per_step = wall_clock_elapsed / SimState.step;
                double WCtFT = (GlobalConfig.max_time - SimState.time) / SimState.dt_global * wall_clock_per_step;
                double WCtMS = (GlobalConfig.max_step - SimState.step) * wall_clock_per_step;
                formattedWrite(writer, "WC=%.1f WCtFT=%.1f WCtMS=%.1f", 
                               wall_clock_elapsed, WCtFT, WCtMS);
                if (GlobalConfig.verbosity_level > 0 && GlobalConfig.is_master_task) {
                    writeln(writer.data);
                    stdout.flush();
                }
                version(mpi_parallel) { MPI_Sync_tasks(); }
                if (GlobalConfig.report_residuals) {
                    // We also compute the residual information and write to screen
                    auto wallClock2 = 1.0e-3*(Clock.currTime() - SimState.wall_clock_start).total!"msecs"();
                    compute_mass_balance(mass_balance);
		    compute_L2_residual(L2_residual);
                    compute_Linf_residuals(Linf_residuals);
                    version(mpi_parallel) {
                        // Reduce residual values across MPI tasks.
                        double my_local_value = Linf_residuals.mass;
                        MPI_Allreduce(MPI_IN_PLACE, &my_local_value, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                        Linf_residuals.mass = my_local_value;
                        my_local_value = Linf_residuals.momentum.x;
                        MPI_Allreduce(MPI_IN_PLACE, &my_local_value, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                        Linf_residuals.momentum.refx = my_local_value;
                        my_local_value = Linf_residuals.momentum.y;
                        MPI_Allreduce(MPI_IN_PLACE, &my_local_value, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                        Linf_residuals.momentum.refy = my_local_value;
                        my_local_value = Linf_residuals.momentum.z;
                        MPI_Allreduce(MPI_IN_PLACE, &my_local_value, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                        Linf_residuals.momentum.refz = my_local_value;
                        my_local_value = Linf_residuals.total_energy;
                        MPI_Allreduce(MPI_IN_PLACE, &my_local_value, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                        Linf_residuals.total_energy = my_local_value;
                        my_local_value = mass_balance;
                        MPI_Allreduce(MPI_IN_PLACE, &my_local_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                        mass_balance = my_local_value;
			my_local_value = L2_residual;
                        MPI_Allreduce(MPI_IN_PLACE, &my_local_value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                        L2_residual = my_local_value;
		    }
                    if (GlobalConfig.is_master_task) {
                        auto writer2 = appender!string();
                        formattedWrite(writer2, "RESIDUALS: step= %7d WC= %.8f ",
                                       SimState.step, wallClock2);
                        formattedWrite(writer2, "MASS: %10.6e X-MOM: %10.6e Y-MOM: %10.6e Z-MOM: %10.6e ENERGY: %10.6e ",
                                       Linf_residuals.mass, Linf_residuals.momentum.x,
                                       Linf_residuals.momentum.y, Linf_residuals.momentum.z, Linf_residuals.total_energy);
                        formattedWrite(writer2, "L2: %10.6e ",
                                       fabs(L2_residual));
                        formattedWrite(writer2, "MASS_BALANCE: %10.6e",
                                       fabs(mass_balance));
                        writeln(writer2.data);
                        stdout.flush();
                    }
                } // end if report_residuals
                version(mpi_parallel) { MPI_Sync_tasks(); }
            } // end if (step...
            //
            // 4.0 (Occasionally) Write out an intermediate solution
            if ((SimState.time >= SimState.t_plot) && !SimState.output_just_written) {
                write_solution_files();
                if (GlobalConfig.udf_supervisor_file.length > 0) { call_UDF_at_write_to_file(); }
                SimState.output_just_written = true;
                SimState.t_plot = SimState.t_plot + GlobalConfig.dt_plot;
                GC.collect();
                GC.minimize();
            }
            //
            // 4.1 (Occasionally) Write out the cell history data and loads on boundary groups data
            if ((SimState.time >= SimState.t_history) && !SimState.history_just_written) {
                write_history_cells_to_files(SimState.time);
                SimState.history_just_written = true;
                SimState.t_history = SimState.t_history + GlobalConfig.dt_history;
                GC.collect();
                GC.minimize();
            }
            if (GlobalConfig.write_loads &&
                (SimState.time >= SimState.t_loads) && !SimState.loads_just_written) {
                write_boundary_loads_to_file(SimState.time, SimState.current_loads_tindx);
                update_loads_times_file(SimState.time, SimState.current_loads_tindx);
                SimState.loads_just_written = true;
                SimState.current_loads_tindx = SimState.current_loads_tindx + 1;
                SimState.t_loads = SimState.t_loads + GlobalConfig.dt_loads;
                GC.collect();
                GC.minimize();
            }
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
                            if (blk.active) { blk.viscous_flux(); }
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
        if(finished_time_stepping && GlobalConfig.verbosity_level >= 1 && GlobalConfig.is_master_task) {
            // Make an announcement about why we are finishing time-stepping.
            write("Integration stopped: ");
            if (caughtException) {
                writeln("An exception was caught during the time step.");
            }
            if (SimState.time >= SimState.target_time) {
                writefln("Reached target simulation time of %g seconds.", SimState.target_time);
            }
            if (SimState.step >= GlobalConfig.max_step) {
                writefln("Reached maximum number of steps with step=%d.", SimState.step);
            }
            if (GlobalConfig.halt_now == 1) { writeln("Halt set in control file."); }
            if (SimState.maxWallClockSeconds > 0 && (wall_clock_elapsed > SimState.maxWallClockSeconds)) {
                writefln("Reached maximum wall-clock time with elapsed time %s.", to!string(wall_clock_elapsed));
            }
            stdout.flush();
        }
    } // end while !finished_time_stepping

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
} // end integrate_in_time()

void check_run_time_configuration(double target_time_as_requested)
{
    // Alter configuration setting if necessary.
    if (GlobalConfig.control_count > 0 && (SimState.step % GlobalConfig.control_count) == 0) {
        read_control_file(); // Reparse the time-step control parameters occasionally.
        SimState.target_time = (GlobalConfig.block_marching) ? target_time_as_requested : GlobalConfig.max_time;
    }
    if (GlobalConfig.viscous && GlobalConfig.viscous_factor < 1.0 &&
        SimState.time > GlobalConfig.viscous_delay) {
        // We want to increment the viscous_factor that scales the viscous effects terms.
        double viscous_factor = GlobalConfig.viscous_factor;
        viscous_factor += GlobalConfig.viscous_factor_increment;
        viscous_factor = min(viscous_factor, 1.0);
        // Make sure that everyone is up-to-date.
        foreach (myblk; localFluidBlocksBySize) {
            myblk.myConfig.viscous_factor = viscous_factor; 
        }
        GlobalConfig.viscous_factor = viscous_factor;
    }
    // We might need to activate or deactivate the IgnitionZones depending on
    // what simulation time we are up to.
    if (SimState.time >= GlobalConfig.ignition_time_start && SimState.time <= GlobalConfig.ignition_time_stop) {
        foreach (blk; localFluidBlocksBySize) { blk.myConfig.ignition_zone_active = true; }
        GlobalConfig.ignition_zone_active = true;
    } else {
        foreach (blk; localFluidBlocksBySize) { blk.myConfig.ignition_zone_active = false; }
        GlobalConfig.ignition_zone_active = false;
    }
} // end check_run_time_configuration()

void synchronize_corner_coords_for_all_blocks()
{
    // Presently, the corner coordinates are only meaningful for structured-grid blocks.
    version(mpi_parallel) {
        // In MPI context, we can only see a subset of the block data.
        foreach (blk; globalFluidBlocks) {
            auto sblk = cast(SFluidBlock) blk;
            if (!sblk) { continue; }
            if (canFind(GlobalConfig.localBlockIds, sblk.id)) {
                // We can see this inside this block to get valid coordinate values.
                sblk.copy_current_corner_coords();
            } else {
                // Cannot see this block so fill in invalid coordinate values.
                sblk.set_current_corner_coords_to_infinity();
            }
        }
        // Now, propagate the valid coordinates across all tasks.
        foreach (blk; globalFluidBlocks) {
            auto sblk = cast(SFluidBlock) blk;
            if (sblk) {
                MPI_Allreduce(MPI_IN_PLACE, sblk.corner_coords.ptr, sblk.corner_coords.length,
                              MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            }
        }
    } else {
        // In shared-memory, we can see all blocks.
        foreach (blk; globalFluidBlocks) {
            auto sblk = cast(SFluidBlock) blk;
            if (sblk) { sblk.copy_current_corner_coords(); }
        }
    }
} // end synchronize_corner_coords_for_all_blocks()

void call_UDF_at_timestep_start()
{
    auto L = GlobalConfig.master_lua_State;
    lua_getglobal(L, "atTimestepStart");
    if (lua_isnil(L, -1)) {
        // There is no suitable Lua function.
        lua_pop(L, 1); // discard the nil item
    } else {
        if (GlobalConfig.user_pad_length > 0) {
            push_array_to_Lua(L, GlobalConfig.userPad, "userPad");
        }
        // Proceed to call the user's function.
        lua_pushnumber(L, SimState.time);
        lua_pushnumber(L, SimState.step);
        lua_pushnumber(L, SimState.dt_global);
        int number_args = 3;
        int number_results = 0;
        if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
            string errMsg = "ERROR: while running user-defined function atTimestepStart()\n";
            errMsg ~= to!string(lua_tostring(L, -1));
            throw new FlowSolverException(errMsg);
        }
        lua_getglobal(L, "dt_override");
        if (lua_isnumber(L, -1)) {
            SimState.dt_override = to!double(lua_tonumber(L, -1));
        } else {
            SimState.dt_override = 0.0;
        }
        lua_pop(L, 1); // dispose dt_override item
        //
        if (GlobalConfig.user_pad_length > 0) {
            get_array_from_Lua(L, GlobalConfig.userPad, "userPad");
        }
    }
    lua_settop(L, 0); // clear stack
    if (GlobalConfig.user_pad_length > 0) {
        broadcast_master_userPad();
        copy_userPad_into_block_interpreters();
    }
} // end call_UDF_at_timestep_start()

void broadcast_master_userPad()
{
    version(mpi_parallel) {
        // The userPad data in master MPI task is broadcast to all other MPI tasks.
        int nelem = to!int(GlobalConfig.userPad.length);
        assert(nelem == GlobalConfig.user_pad_length, "Oops, wrong lengths");
        // We allocate the array once.
        if (userPad_copy.length < nelem) { userPad_copy.length = nelem; }
        foreach (i, elem; GlobalConfig.userPad) { userPad_copy[i] = elem; }
        MPI_Bcast(userPad_copy.ptr, nelem, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (GlobalConfig.mpi_rank_for_local_task > 0) {
            foreach (i, elem; userPad_copy) { GlobalConfig.userPad[i] = elem; }
        }
    }
} // end broadcast_master_userPad()

void copy_userPad_into_block_interpreters()
{
    // Within the one task, broadcast userPad to the Lua interpreters
    // associated with the blocks and boundary-conditions.
    foreach (blk; localFluidBlocks) {
        push_array_to_Lua(blk.myL, GlobalConfig.userPad, "userPad");
        foreach (bc; blk.bc) {
            if (bc.myL) { push_array_to_Lua(bc.myL, GlobalConfig.userPad, "userPad"); }
        }
    }
} // end copy_userPad_into_block_interpreters()

void call_UDF_at_timestep_end()
{
    auto L = GlobalConfig.master_lua_State;
    lua_getglobal(L, "atTimestepEnd");
    if (lua_isnil(L, -1)) {
        // There is no suitable Lua function.
        lua_pop(L, 1); // discard the nil item
    } else {
        if (GlobalConfig.user_pad_length > 0) {
            push_array_to_Lua(L, GlobalConfig.userPad, "userPad");
        }
        // Proceed to call the user's function.
        lua_pushnumber(L, SimState.time);
        lua_pushnumber(L, SimState.step);
        lua_pushnumber(L, SimState.dt_global);
        int number_args = 3;
        int number_results = 0;
        if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
            string errMsg = "ERROR: while running user-defined function atTimestepEnd()\n";
            errMsg ~= to!string(lua_tostring(L, -1));
            throw new FlowSolverException(errMsg);
        }
        if (GlobalConfig.user_pad_length > 0) {
            get_array_from_Lua(L, GlobalConfig.userPad, "userPad");
        }
    }
    lua_settop(L, 0); // clear stack
    if (GlobalConfig.user_pad_length > 0) {
        broadcast_master_userPad();
        copy_userPad_into_block_interpreters();
    }
} // end call_UDF_at_timestep_end()

void call_UDF_at_write_to_file()
{
    // The user may also have some writing of data to do via their Lua script file.
    // This function is called just after writing the flow solution to file.
    auto L = GlobalConfig.master_lua_State;
    lua_getglobal(L, "atWriteToFile");
    if (lua_isnil(L, -1)) {
        // There is no suitable Lua function.
        lua_pop(L, 1); // discard the nil item
    } else {
        if (GlobalConfig.user_pad_length > 0) {
            push_array_to_Lua(L, GlobalConfig.userPad, "userPad");
        }
        //
        // Proceed to call the user's function.
        lua_pushnumber(L, SimState.time);
        lua_pushnumber(L, SimState.step);
        lua_pushnumber(L, SimState.dt_global);
        int number_args = 3;
        int number_results = 0;
        if ( lua_pcall(L, number_args, number_results, 0) != 0 ) {
            string errMsg = "ERROR: while running user-defined function atWriteToFile()\n";
            errMsg ~= to!string(lua_tostring(L, -1));
            throw new FlowSolverException(errMsg);
        }
        if (GlobalConfig.user_pad_length > 0) {
            get_array_from_Lua(L, GlobalConfig.userPad, "userPad");
        }
    }
    lua_settop(L, 0); // clear stack
    if (GlobalConfig.user_pad_length > 0) {
        broadcast_master_userPad();
        copy_userPad_into_block_interpreters();
    }
} // end call_UDF_at_write_to_file()

void determine_time_step_size()
{
    // Set the size of the time step to be the minimum allowed for any active block.
    // We will check it occasionally, if we have not elected to keep fixed time steps. 
    bool do_dt_check_now =
        ((SimState.step % GlobalConfig.cfl_count) == 0) ||
        (SimState.dt_override > 0.0);
    version(mpi_parallel) {
        // If one task is doing a time-step check, all tasks have to.
        int myFlag = to!int(do_dt_check_now);
        MPI_Allreduce(MPI_IN_PLACE, &myFlag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        do_dt_check_now = to!bool(myFlag);
    }
    if (do_dt_check_now) {
        // Adjust the time step...
        // First, check what each block thinks should be the allowable step size.
        // Also, if we have done some steps, check the CFL limits.
        foreach (i, myblk; parallel(localFluidBlocksBySize,1)) {
            // Note 'i' is not necessarily the block id but
            // that is not important here, just need a unique spot to poke into local_dt_allow.
            if (myblk.active) {
                local_dt_allow[i] = myblk.determine_time_step_size(SimState.dt_global,
                                                                   (SimState.step > 0))[0];
                local_cfl_max[i] = myblk.determine_time_step_size(SimState.dt_global,
                                                                   (SimState.step > 0))[1];
                local_dt_allow_parab[i] = myblk.determine_time_step_size(SimState.dt_global,
                                                                   (SimState.step > 0))[2];
            }
        }
        // Second, reduce this estimate across all local blocks.
	if (GlobalConfig.with_super_time_stepping) {
	    if (GlobalConfig.super_step_hyperbolic) SimState.dt_allow = 0.0; // to be sure it is replaced.
	    else SimState.dt_allow = double.max; // to be sure it is replaced.
	    SimState.cfl_max = 0.0; // to be sure it is replaced.
	} else {
	    SimState.dt_allow = double.max; //0.0; // to be sure it is replaced.
	    SimState.cfl_max = 0.0; // to be sure it is replaced.
	}
	SimState.dt_allow_parab = double.max;
        foreach (i, myblk; localFluidBlocks) { // serial loop
            if (myblk.active) { 
		if (GlobalConfig.with_super_time_stepping) {
		    if (GlobalConfig.super_step_hyperbolic) SimState.dt_allow = max(SimState.dt_allow, local_dt_allow[i]);
		    else  SimState.dt_allow = min(SimState.dt_allow, local_dt_allow[i]);
		    SimState.dt_allow_parab = min(SimState.dt_allow_parab, local_dt_allow_parab[i]);
		}
		else SimState.dt_allow = min(SimState.dt_allow, local_dt_allow[i]); 
	    } 
            if (myblk.active) { SimState.cfl_max = max(SimState.cfl_max, local_cfl_max[i]); } 
        }
        version(mpi_parallel) {
            double my_dt_allow = SimState.dt_allow;
            double my_cfl_max = SimState.cfl_max;
            MPI_Allreduce(MPI_IN_PLACE, &my_dt_allow, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &my_cfl_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            SimState.dt_allow = my_dt_allow;
            SimState.cfl_max = my_cfl_max;
        }
	if (GlobalConfig.with_super_time_stepping) { 
            if (SimState.step == 0) {
                // When starting out, we may override the computed value.
                // This might be handy for situations where the computed estimate
                // is likely to be not small enough for numerical stability.
                SimState.dt_allow = fmin(GlobalConfig.dt_init, SimState.dt_allow);
		SimState.dt_allow_parab = fmin(GlobalConfig.dt_init, SimState.dt_allow_parab);
            }
            // Now, change the actual time step, as needed.
            if (SimState.dt_allow <= SimState.dt_global) {
                // If we need to reduce the time step, do it immediately.
                //SimState.dt_global = SimState.dt_allow;
	    } else {
                // Make the transitions to larger time steps gentle.
		SimState.dt_global = min(SimState.dt_global*1.5, SimState.dt_allow);
		// The user may supply, explicitly, a maximum time-step size.
                SimState.dt_global = min(SimState.dt_global, GlobalConfig.dt_max);
	    }
            if (SimState.dt_allow_parab <= SimState.dt_global_parab) {
                // If we need to reduce the time step, do it immediately.
                //SimState.dt_global_parab = SimState.dt_allow_parab;
            } else {
                // Make the transitions to larger time steps gentle.
                SimState.dt_global_parab = min(SimState.dt_global_parab*1.5, SimState.dt_allow_parab);
                // The user may supply, explicitly, a maximum time-step size.
                SimState.dt_global_parab = min(SimState.dt_global_parab, GlobalConfig.dt_max);
            }
	    //SimState.dt_global = SimState.dt_allow;
	    //SimState.dt_global_parab = SimState.dt_allow_parab;
	} else if (GlobalConfig.with_local_time_stepping) { SimState.dt_global = SimState.dt_allow; }
	else { // do some global time-stepping checks
            if (SimState.step == 0) {
                // When starting out, we may override the computed value.
                // This might be handy for situations where the computed estimate
                // is likely to be not small enough for numerical stability.
                SimState.dt_allow = fmin(GlobalConfig.dt_init, SimState.dt_allow);
            }
            if (SimState.dt_override > 0.0) {
                // The user-defined supervisory function atTimestepStart may have set
                // dt_override because it knows something about the simulation conditions
                // that is not handled well by our generic CFL check.
                SimState.dt_allow = fmin(SimState.dt_override, SimState.dt_allow);
            }
            // Now, change the actual time step, as needed.
            if (SimState.dt_allow <= SimState.dt_global) {
                // If we need to reduce the time step, do it immediately.
                SimState.dt_global = SimState.dt_allow;
            } else {
                // Make the transitions to larger time steps gentle.
                SimState.dt_global = min(SimState.dt_global*1.5, SimState.dt_allow);
                // The user may supply, explicitly, a maximum time-step size.
                SimState.dt_global = min(SimState.dt_global, GlobalConfig.dt_max);
            }
        }
    } // end if do_cfl_check_now 
} // end determine_time_step_size()

void update_ch_for_divergence_cleaning()
{
    bool first = true;
    foreach (blk; localFluidBlocksBySize) {
        if (!blk.active) continue;
        if (first) {
            GlobalConfig.c_h = blk.update_c_h(SimState.dt_global);
            first = false;
        } else {
            GlobalConfig.c_h = fmin(blk.update_c_h(SimState.dt_global), GlobalConfig.c_h);
        }
    }
    version(mpi_parallel) {
        double my_c_h = GlobalConfig.c_h;
        MPI_Allreduce(MPI_IN_PLACE, &my_c_h, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        GlobalConfig.c_h = my_c_h;
    }
    // Now that we have a globally-reduced value, propagate that new value
    // into the block-local config structure.
    foreach (blk; localFluidBlocksBySize) {
        if (!blk.active) continue;
        blk.myConfig.c_h = GlobalConfig.c_h;
    }
} // end update_ch_for_divergence_cleaning()

void k_omega_set_mu_and_k()
{
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        if (blk.active) {
	    blk.flow_property_spatial_derivatives(0); 
            blk.estimate_turbulence_viscosity();
        }
    }
} // end k_omega_set_mu_and_k()

void chemistry_step(double dt)
{
    version (gpu_chem) {
        if (GlobalConfig.with_local_time_stepping)
            assert(0, "Oops, GPU accelerated chemistry and LTS aren't currently compatible.");
        GlobalConfig.gpuChem.thermochemical_increment(dt);
    } else {
        // without GPU accelerator
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) {
                double local_dt = dt;
                if (GlobalConfig.with_local_time_stepping) foreach (cell; blk.cells) { cell.thermochemical_increment(cell.dt_local); }
                else foreach (cell; blk.cells) { cell.thermochemical_increment(local_dt); }
            }
        }
    }
} // end chemistry_half_step()

void set_grid_velocities()
{
    final switch(GlobalConfig.grid_motion){
        case GridMotion.none:
            throw new Error("Should not be setting grid velocities in with GridMotion.none");
        case GridMotion.user_defined:
            // Rely on user to set vertex velocities.
            // Note that velocities remain unchanged if the user does nothing.
            assign_vertex_velocities_via_udf(SimState.time, SimState.dt_global);
            break;
        case GridMotion.shock_fitting:
            /++if (GlobalConfig.in_mpi_context) {
                throw new Error("oops, shock_fitting not compatible with MPI.");
                // [TODO] 2018-01-20 PJ should do something to lift this restriction.
            }++/
            // apply boundary conditions here because ...
            // shockfitting algorithm requires ghost cells to be up to date.
            if (SimState.time > GlobalConfig.shock_fitting_delay) {
                exchange_ghost_cell_boundary_data(SimState.time, 0, 0);
                foreach (blk; localFluidBlocksBySize) {
                    if (blk.active) { blk.applyPreReconAction(SimState.time, 0, 0); }
                }
                foreach (blk; localFluidBlocksBySize) {
                    SFluidBlock sblk = cast(SFluidBlock) blk;
                    
                    version(mpi_parallel) {
                        // Define a new communication group- need to be real explicit with MPI
                        assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                        MPI_Comm shock_fitting_comm;
                        int communicator_tag = 0, local_rank;
                        MPI_Comm_rank(MPI_COMM_WORLD, &local_rank);
                        foreach (indx; sblk.inflow_partners) {
                            communicator_tag += indx;
                        }
                        MPI_Comm_split(MPI_COMM_WORLD, communicator_tag, local_rank, &shock_fitting_comm);
                        // Define the variables that will be local to some processes
                        int[2] is_master_block; // Bit passed around to determine master rank
                        int[2] master_rank; // This contains the location of the shock fitted block
                        int shock_fitting_rank; // Rank of local process
                        double[] unpacked_vertex_velocities; // Array that will contain the vertex velocities broken into their spatial elements
                        Vector3[] inflow_vertex_velocities, packed_vertex_velocities; // Array that will contain the vertex velocities in vector form
                        MPI_Comm_rank(shock_fitting_comm, &shock_fitting_rank); // Assign the local rank
                        is_master_block[1] = shock_fitting_rank;

                        // Calculate the vertex velocities on the inflow boundary
                        if (blk.active && blk.bc[Face.west].type == "inflow_shock_fitting") {
                            is_master_block[0] = 1;     // Assign self as the master block
                            inflow_vertex_velocities = shock_fitting_vertex_velocities(sblk);   // Do calculation
                            assert(inflow_vertex_velocities.length == (sblk.jmax - sblk.jmin + 2), "the vertex velocity array is the wrong size");
                            unpacked_vertex_velocities = unpack_vertex_velocities(inflow_vertex_velocities); // Unpack the vertex velocities to a form that can be read by MPI
                        }

                        MPI_Allreduce(is_master_block.ptr, master_rank.ptr, 2, MPI_2INT, MPI_MAXLOC, shock_fitting_comm);     // Tell everyone who their master is
                        int ne = to!int(unpacked_vertex_velocities.length);     // Number of orders
                        MPI_Bcast(&ne, 1, MPI_INT, master_rank[1], shock_fitting_comm);     // Let all the other blocks know how much space to set aside
                        unpacked_vertex_velocities.length = ne;                             // Set aside that space
                        MPI_Bcast(unpacked_vertex_velocities.ptr, ne, MPI_DOUBLE, master_rank[1], shock_fitting_comm);  // Master gives his orders
                        packed_vertex_velocities = pack_vertex_velocities(unpacked_vertex_velocities);      // Put orders into usable format
                        assign_slave_velocities(sblk, packed_vertex_velocities);    // Assign inner vertex velocities
                        MPI_Comm_free(&shock_fitting_comm);
                    }
                    else {
                        if (blk.active && blk.bc[Face.west].type == "inflow_shock_fitting" && SimState.time > GlobalConfig.shock_fitting_delay) {
                            assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                            Vector3[] inflow_vertex_velocities = shock_fitting_vertex_velocities(sblk);
                            assert(inflow_vertex_velocities.length == (sblk.jmax - sblk.jmin + 2), "the vertex velocity array is the wrong size");

                            foreach (indx; sblk.inflow_partners) {
                                assign_slave_velocities(cast(SFluidBlock) globalFluidBlocks[indx], inflow_vertex_velocities);
                            }
                        }
                    }
                }
            }              
    }
} // end set_grid_velocities()

void recalculate_all_geometry()
{
    // Moving Grid - Recalculate all geometry, note that in the gas dynamic
    // update gtl level 2 is copied to gtl level 0 for the next step thus
    // we actually do want to calculate geometry at gtl 0 here.
    foreach (blk; localFluidBlocksBySize) {
        if (blk.active) {
            blk.compute_primary_cell_geometric_data(0);
            blk.compute_least_squares_setup(0);
        }
    }
} // end recalculate_all_geometry()

//---------------------------------------------------------------------------

void exchange_ghost_cell_boundary_data(double t, int gtl, int ftl)
// We have hoisted the exchange of ghost-cell data out of the GhostCellEffect class
// that used to live only inside the boundary condition attached to a block.
// The motivation for allowing this leakage of abstraction is that the MPI
// exchange of messages requires a coordination of actions that spans blocks.
// Sigh...  2017-01-24 PJ
// p.s. The data for that coordination is still buried in the FullFaceCopy class.
// No need to have all its guts hanging out.
{
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_geometry_phase0(); }
                auto mygce2 = cast(GhostCellFullFaceCopy) gce;
                if (mygce2) { mygce2.exchange_geometry_phase0(); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_geometry_phase1(); }
                auto mygce2 = cast(GhostCellFullFaceCopy) gce;
                if (mygce2) { mygce2.exchange_geometry_phase1(); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_geometry_phase2(); }
                auto mygce2 = cast(GhostCellFullFaceCopy) gce;
                if (mygce2) { mygce2.exchange_geometry_phase2(); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_flowstate_phase0(t, gtl, ftl); }
                auto mygce2 = cast(GhostCellFullFaceCopy) gce;
                if (mygce2) { mygce2.exchange_flowstate_phase0(t, gtl, ftl); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_flowstate_phase1(t, gtl, ftl); }
                auto mygce2 = cast(GhostCellFullFaceCopy) gce;
                if (mygce2) { mygce2.exchange_flowstate_phase1(t, gtl, ftl); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach(bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce1 = cast(GhostCellMappedCellCopy) gce;
                if (mygce1) { mygce1.exchange_flowstate_phase2(t, gtl, ftl); }
                auto mygce2 = cast(GhostCellFullFaceCopy) gce;
                if (mygce2) { mygce2.exchange_flowstate_phase2(t, gtl, ftl); }
            }
        }
    }
} // end exchange_ghost_cell_boundary_data()

//----------------------------------------------------------------------------
void sts_gasdynamic_explicit_increment_with_fixed_grid()
{
    shared double t0 = SimState.time;
    shared bool with_local_time_stepping = GlobalConfig.with_local_time_stepping;
    shared bool with_k_omega = (GlobalConfig.turbulence_model == TurbulenceModel.k_omega) &&
        !GlobalConfig.separate_update_for_k_omega_source;
    shared bool allow_high_order_interpolation = (SimState.time >= GlobalConfig.interpolation_delay);
    
    // compute number of sub-cycle stages (S)
    double alpha = (SimState.dt_global)/(SimState.dt_global_parab);
    double s = 0.5*(-1.0+sqrt(1+8.0*alpha));
    s = floor(s);
    if (fmod(s, 2) == 0.0 && s != 1.0) {
        s = s - 1.0;
    }
    
    // store as an int for looping later
    int S = to!int(s);

    // if S is less than 1 then just set S = 1 (Euler step)
    if (S <= 1) {
	S = 1;
	s = 1;
	SimState.dt_global = SimState.dt_global_parab;
    }
    else {
	// since we round S down to the nearest odd integer we should alter the time-step to be consistent
	SimState.dt_global = SimState.dt_global_parab * (S*S+S)/(2.0);
    }

    // for temporal MMS uncomment the following lines
    //S = 4;
    //s = 4;
    //SimState.dt_global = GlobalConfig.dt_init;
    
    //writeln("stages: ", S, ", dt_hyperbolic = ", SimState.dt_global, ", dt_parabolic = ", SimState.dt_global_parab, ", alpha = ", alpha);
    // --------------------------------------------------
    // j = 1
    // --------------------------------------------------
    int j = 1;

    // First-stage of gas-dynamic update.
    shared int ftl = 0; // time-level within the overall convective-update
    shared int gtl = 0; // grid time-level remains at zero for the non-moving grid

    // Preparation for the inviscid gas-dynamic flow update.
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
	if (blk.active) {
	    blk.clear_fluxes_of_conserved_quantities();
	    foreach (cell; blk.cells) {
		cell.clear_source_vector();
		cell.data_is_bad = false;
	    }
	}
    }
    exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
    if (GlobalConfig.apply_bcs_in_parallel) {
	foreach (blk; parallel(localFluidBlocksBySize,1)) {
	    if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
	}
    } else {
	foreach (blk; localFluidBlocksBySize) {
	    if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
	}
    }
    // And we'll do a first pass on solid domain bc's too in case we need up-to-date info.
    foreach (sblk; parallel(solidBlocks, 1)) {
	if (sblk.active) { sblk.applyPreSpatialDerivAction(SimState.time, ftl); }
    }
    // We've put this detector step here because it needs the ghost-cell data
    // to be current, as it should be just after a call to apply_convective_bc().
    if (GlobalConfig.flux_calculator == FluxCalculator.adaptive_hanel_ausmdv ||
        GlobalConfig.flux_calculator == FluxCalculator.adaptive_hlle_roe ||
	GlobalConfig.flux_calculator == FluxCalculator.adaptive_efm_ausmdv) {
	foreach (blk; parallel(localFluidBlocksBySize,1)) {
	    if (blk.active) { blk.detect_shock_points(); }
	}
    }
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
	if (blk.active) { blk.convective_flux_phase0(allow_high_order_interpolation, gtl); }
    }
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
	if (blk.active) { blk.convective_flux_phase1(allow_high_order_interpolation, gtl); }
    }
    if (GlobalConfig.apply_bcs_in_parallel) {
	foreach (blk; parallel(localFluidBlocksBySize,1)) {
	    if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
	}
    } else {
	foreach (blk; localFluidBlocksBySize) {
	    if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
	}
    }
    if (GlobalConfig.viscous && !GlobalConfig.separate_update_for_viscous_terms) {
	if (GlobalConfig.apply_bcs_in_parallel) {
	    foreach (blk; parallel(localFluidBlocksBySize,1)) {
		if (blk.active) {
		    blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
		    blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
		}
	    }
	} else {
	    foreach (blk; localFluidBlocksBySize) {
		if (blk.active) {
		    blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
		    blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
		}
	    }
	}
	foreach (blk; parallel(localFluidBlocksBySize,1)) {
	    if (blk.active) {
		blk.flow_property_spatial_derivatives(gtl); 
		blk.estimate_turbulence_viscosity();
	    }
	}
	// we exchange boundary data at this point to ensure the
	// ghost cells along block-block boundaries have the most
	// recent mu_t and k_t values.
	exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
	foreach (blk; parallel(localFluidBlocksBySize,1)) {
	    if (blk.active) {         
		blk.viscous_flux();
	    }
	}
	if (GlobalConfig.apply_bcs_in_parallel) {
	    foreach (blk; parallel(localFluidBlocksBySize,1)) {
		if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
	    }
	} else {
	    foreach (blk; localFluidBlocksBySize) {
		if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
	    }
	}
    } // end if viscous
    foreach (i, blk; parallel(localFluidBlocksBySize,1)) {
	if (!blk.active) continue;
	int local_ftl = ftl;
	int local_gtl = gtl;
	bool local_with_k_omega = with_k_omega;
	bool local_with_local_time_stepping = with_local_time_stepping;
	double local_dt_global = SimState.dt_global;
	double local_sim_time = SimState.time;
	foreach (cell; blk.cells) {
	    cell.add_inviscid_source_vector(local_gtl, blk.omegaz);
	    if (blk.myConfig.viscous && !blk.myConfig.separate_update_for_viscous_terms) {
		cell.add_viscous_source_vector(local_with_k_omega);
	    }
	    if (blk.myConfig.udf_source_terms) { // [TODO] may want to apply serially
		size_t i_cell = cell.id;
		size_t j_cell = 0;
		size_t k_cell = 0;
		if (blk.grid_type == Grid_t.structured_grid) {
		    auto sblk = cast(SFluidBlock) blk;
		    assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
		    auto ijk_indices = sblk.cell_id_to_ijk_indices(cell.id);
		    i_cell = ijk_indices[0];
		    j_cell = ijk_indices[1];
		    k_cell = ijk_indices[2];
		}
		addUDFSourceTermsToCell(blk.myL, cell, local_gtl, 
					local_sim_time, blk.myConfig.gmodel,
					blk.id, i_cell, j_cell, k_cell);
	    }
	    cell.time_derivatives(local_gtl, local_ftl, local_with_k_omega);
	}
	if (blk.myConfig.residual_smoothing) { blk.residual_smoothing_dUdt(local_ftl); }
	bool force_euler = false;
	foreach (cell; blk.cells) {
	    cell.rkl1_stage_1_update_for_flow_on_fixed_grid(local_dt_global, local_with_k_omega, j, S, false);
	    cell.decode_conserved(local_gtl, local_ftl+1, blk.omegaz);
	} // end foreach cell
	local_invalid_cell_count[i] = blk.count_invalid_cells(local_gtl, local_ftl+1);
    } // end foreach blk
    //
    int flagTooManyBadCells = 0;
    foreach (i, blk; localFluidBlocksBySize) { // serial loop
	if (local_invalid_cell_count[i] > GlobalConfig.max_invalid_cells) {
	    flagTooManyBadCells = 1;
	    writefln("Following first-stage gasdynamic update: %d bad cells in block[%d].",
		     local_invalid_cell_count[i], i);
	}
    }
    version(mpi_parallel) {
	MPI_Allreduce(MPI_IN_PLACE, &flagTooManyBadCells, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
    if (flagTooManyBadCells > 0) {
	throw new FlowSolverException("Too many bad cells; go home.");
    }
    //
    if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight) {
	// Next do solid domain update IMMEDIATELY after at same flow time level
	if (GlobalConfig.in_mpi_context && solidBlocks.length > 0) {
	    throw new Error("oops, should not be doing coupling with solid domains in MPI.");
	    // [TODO] 2018-01-20 PJ should do something to lift this restriction.
	}
	foreach (sblk; parallel(solidBlocks, 1)) {
	    if (!sblk.active) continue;
	    sblk.clearSources();
	    sblk.computeSpatialDerivatives(ftl);
	    sblk.computeFluxes();
	}
	if (GlobalConfig.apply_bcs_in_parallel) {
	    foreach (sblk; parallel(solidBlocks, 1)) {
		if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
	    }
	} else {
	    foreach (sblk; solidBlocks) {
		    if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
	    }
	}
	// We need to synchronise before updating
	foreach (sblk; parallel(solidBlocks, 1)) {
	    foreach (scell; sblk.activeCells) {
		if (GlobalConfig.udfSolidSourceTerms) {
		    addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time);
		}
		scell.timeDerivatives(ftl, GlobalConfig.dimensions);
		scell.stage1Update(SimState.dt_global);
		scell.T = updateTemperature(scell.sp, scell.e[ftl+1]);
	    } // end foreach scell
	} // end foreach sblk
    } // end if tight solid domain coupling.
    
    // --------------------------------------------------
    // 2 <= j <= S
    // --------------------------------------------------
    foreach (J; 2..S+1) {
	j = J;
	ftl = 1;
	
	// Preparation for the inviscid gas-dynamic flow update.
	foreach (blk; parallel(localFluidBlocksBySize,1)) {
	    if (blk.active) {
		blk.clear_fluxes_of_conserved_quantities();
		foreach (cell; blk.cells) {
		    cell.clear_source_vector();
		    cell.data_is_bad = false;
		}
	    }
	}
	exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
	if (GlobalConfig.apply_bcs_in_parallel) {
	    foreach (blk; parallel(localFluidBlocksBySize,1)) {
		if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
	    }
	} else {
	    foreach (blk; localFluidBlocksBySize) {
		if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
	    }
	}
	// And we'll do a first pass on solid domain bc's too in case we need up-to-date info.
	foreach (sblk; parallel(solidBlocks, 1)) {
	    if (sblk.active) { sblk.applyPreSpatialDerivAction(SimState.time, ftl); }
	}
	// We've put this detector step here because it needs the ghost-cell data
	// to be current, as it should be just after a call to apply_convective_bc().
	if (GlobalConfig.flux_calculator == FluxCalculator.adaptive_hanel_ausmdv ||
	    GlobalConfig.flux_calculator == FluxCalculator.adaptive_hlle_roe ||
	    GlobalConfig.flux_calculator == FluxCalculator.adaptive_efm_ausmdv) {
	    foreach (blk; parallel(localFluidBlocksBySize,1)) {
		if (blk.active) { blk.detect_shock_points(); }
	    }
	}
	foreach (blk; parallel(localFluidBlocksBySize,1)) {
	    if (blk.active) { blk.convective_flux_phase0(allow_high_order_interpolation, gtl); }
	}
	foreach (blk; parallel(localFluidBlocksBySize,1)) {
	    if (blk.active) { blk.convective_flux_phase1(allow_high_order_interpolation, gtl); }
	}
	if (GlobalConfig.apply_bcs_in_parallel) {
	    foreach (blk; parallel(localFluidBlocksBySize,1)) {
		if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
	    }
	} else {
	    foreach (blk; localFluidBlocksBySize) {
		if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
	    }
	}
	if (GlobalConfig.viscous && !GlobalConfig.separate_update_for_viscous_terms) {
	    if (GlobalConfig.apply_bcs_in_parallel) {
		foreach (blk; parallel(localFluidBlocksBySize,1)) {
		    if (blk.active) {
			blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
			blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
		    }
		}
	    } else {
		foreach (blk; localFluidBlocksBySize) {
		    if (blk.active) {
			blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
			blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
		    }
		}
	    }
	    foreach (blk; parallel(localFluidBlocksBySize,1)) {
		if (blk.active) {
		    blk.flow_property_spatial_derivatives(gtl); 
		    blk.estimate_turbulence_viscosity();
		}
	    }
	    // we exchange boundary data at this point to ensure the
	    // ghost cells along block-block boundaries have the most
	    // recent mu_t and k_t values.
	    exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
	    foreach (blk; parallel(localFluidBlocksBySize,1)) {
		if (blk.active) {         
		    blk.viscous_flux();
		}
	    }
	    if (GlobalConfig.apply_bcs_in_parallel) {
		foreach (blk; parallel(localFluidBlocksBySize,1)) {
		    if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
		}
	    } else {
		foreach (blk; localFluidBlocksBySize) {
		    if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
		}
	    }
	} // end if viscous
	foreach (i, blk; parallel(localFluidBlocksBySize,1)) {
	    if (!blk.active) continue;
	    int local_ftl = ftl;
	    int local_gtl = gtl;
	    bool local_with_k_omega = with_k_omega;
	    bool local_with_local_time_stepping = with_local_time_stepping;
	    double local_dt_global = SimState.dt_global;
	    double local_sim_time = SimState.time;
	    foreach (cell; blk.cells) {
		cell.add_inviscid_source_vector(local_gtl, blk.omegaz);
		if (blk.myConfig.viscous && !blk.myConfig.separate_update_for_viscous_terms) {
		    cell.add_viscous_source_vector(local_with_k_omega);
		}
		if (blk.myConfig.udf_source_terms) { // [TODO] may want to apply serially
		    size_t i_cell = cell.id;
		    size_t j_cell = 0;
		    size_t k_cell = 0;
		    if (blk.grid_type == Grid_t.structured_grid) {
			auto sblk = cast(SFluidBlock) blk;
			assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
			auto ijk_indices = sblk.cell_id_to_ijk_indices(cell.id);
			i_cell = ijk_indices[0];
			j_cell = ijk_indices[1];
			k_cell = ijk_indices[2];
		    }
		    addUDFSourceTermsToCell(blk.myL, cell, local_gtl, 
					local_sim_time, blk.myConfig.gmodel,
					    blk.id, i_cell, j_cell, k_cell);
		}
		cell.time_derivatives(local_gtl, local_ftl, local_with_k_omega);
	    }
	    if (blk.myConfig.residual_smoothing) { blk.residual_smoothing_dUdt(local_ftl); }
	    bool force_euler = false;
	    foreach (cell; blk.cells) {
		cell.rkl1_stage_j_update_for_flow_on_fixed_grid(local_dt_global, local_with_k_omega, j, S, false);
		cell.decode_conserved(local_gtl, local_ftl+1, blk.omegaz);
	    } // end foreach cell
	    local_invalid_cell_count[i] = blk.count_invalid_cells(local_gtl, local_ftl+1);
	} // end foreach blk
	//
	flagTooManyBadCells = 0;
	foreach (i, blk; localFluidBlocksBySize) { // serial loop
	    if (local_invalid_cell_count[i] > GlobalConfig.max_invalid_cells) {
		flagTooManyBadCells = 1;
		writefln("Following first-stage gasdynamic update: %d bad cells in block[%d].",
			 local_invalid_cell_count[i], i);
	    }
	}
	version(mpi_parallel) {
	    MPI_Allreduce(MPI_IN_PLACE, &flagTooManyBadCells, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
	}
	if (flagTooManyBadCells > 0) {
	    throw new FlowSolverException("Too many bad cells; go home.");
	}
	//
	if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight) {
	    // Next do solid domain update IMMEDIATELY after at same flow time level
	    if (GlobalConfig.in_mpi_context && solidBlocks.length > 0) {
		throw new Error("oops, should not be doing coupling with solid domains in MPI.");
		// [TODO] 2018-01-20 PJ should do something to lift this restriction.
	    }
	    foreach (sblk; parallel(solidBlocks, 1)) {
		if (!sblk.active) continue;
		sblk.clearSources();
		sblk.computeSpatialDerivatives(ftl);
		sblk.computeFluxes();
	    }
	    if (GlobalConfig.apply_bcs_in_parallel) {
		foreach (sblk; parallel(solidBlocks, 1)) {
		    if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
		}
	    } else {
		foreach (sblk; solidBlocks) {
		    if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
		}
	    }
	    // We need to synchronise before updating
	    foreach (sblk; parallel(solidBlocks, 1)) {
		foreach (scell; sblk.activeCells) {
		    if (GlobalConfig.udfSolidSourceTerms) {
			addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time);
		    }
		    scell.timeDerivatives(ftl, GlobalConfig.dimensions);
		    scell.stage1Update(SimState.dt_global);
		    scell.T = updateTemperature(scell.sp, scell.e[ftl+1]);
		} // end foreach scell
	    } // end foreach sblk
	} // end if tight solid domain coupling.

	// shuffle time-levels for next iteration (U1 goes to U0 & U2 goes to U1)
	foreach (blk; parallel(localFluidBlocksBySize,1)) {
	    if (blk.active) {
		foreach (cell; blk.cells) { 
		    swap(cell.U[1], cell.U[2]);
		    swap(cell.U[0], cell.U[2]); 
		}
	    }
	} // end foreach blk
	
    } // end foreach (J; 2..S+1)
    
    // Get the end conserved data into U[0] for next step.
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        if (blk.active) {
	    size_t end_indx = 1; // time-level holds current solution
	    foreach (cell; blk.cells) { swap(cell.U[0], cell.U[end_indx]); }
        }
    } // end foreach blk
    //
    foreach (sblk; solidBlocks) {
        if (sblk.active) {
            size_t end_indx = 1; // time-level holds current solution
            foreach (scell; sblk.activeCells) { scell.e[0] = scell.e[end_indx]; }
        }
    } // end foreach sblk
    
    // Finally, update the globally known simulation time for the whole step.
    SimState.time = t0 + SimState.dt_global;
} // end gasdynamic_explicit_increment_with_fixed_grid()

void gasdynamic_explicit_increment_with_fixed_grid()
{
    shared double t0 = SimState.time;
    shared bool with_local_time_stepping = GlobalConfig.with_local_time_stepping;
    shared bool with_k_omega = (GlobalConfig.turbulence_model == TurbulenceModel.k_omega) &&
        !GlobalConfig.separate_update_for_k_omega_source;
    shared bool allow_high_order_interpolation = (SimState.time >= GlobalConfig.interpolation_delay);
    // Set the time-step coefficients for the stages of the update scheme.
    shared double c2 = 1.0; // default for predictor-corrector update
    shared double c3 = 1.0; // default for predictor-corrector update
    final switch ( GlobalConfig.gasdynamic_update_scheme ) {
    case GasdynamicUpdate.euler:
    case GasdynamicUpdate.pc: c2 = 1.0; c3 = 1.0; break;
    case GasdynamicUpdate.midpoint: c2 = 0.5; c3 = 1.0; break;
    case GasdynamicUpdate.classic_rk3: c2 = 0.5; c3 = 1.0; break;
    case GasdynamicUpdate.tvd_rk3: c2 = 1.0; c3 = 0.5; break;
    case GasdynamicUpdate.denman_rk3: c2 = 1.0; c3 = 0.5; break;
    case GasdynamicUpdate.moving_grid_1_stage:
    case GasdynamicUpdate.moving_grid_2_stage: assert(false, "invalid option");
    }
    // Preparation for the predictor-stage of inviscid gas-dynamic flow update.
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        if (blk.active) {
            blk.clear_fluxes_of_conserved_quantities();
            foreach (cell; blk.cells) {
                cell.clear_source_vector();
                cell.data_is_bad = false;
            }
        }
    }
    // First-stage of gas-dynamic update.
    shared int ftl = 0; // time-level within the overall convective-update
    shared int gtl = 0; // grid time-level remains at zero for the non-moving grid
    exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
    if (GlobalConfig.apply_bcs_in_parallel) {
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
        }
    } else {
        foreach (blk; localFluidBlocksBySize) {
            if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
        }
    }
    // And we'll do a first pass on solid domain bc's too in case we need up-to-date info.
    foreach (sblk; parallel(solidBlocks, 1)) {
        if (sblk.active) { sblk.applyPreSpatialDerivAction(SimState.time, ftl); }
    }
    // We've put this detector step here because it needs the ghost-cell data
    // to be current, as it should be just after a call to apply_convective_bc().
    if (GlobalConfig.flux_calculator == FluxCalculator.adaptive_hanel_ausmdv ||
        GlobalConfig.flux_calculator == FluxCalculator.adaptive_hlle_roe ||
        GlobalConfig.flux_calculator == FluxCalculator.adaptive_efm_ausmdv) {
	foreach (blk; parallel(localFluidBlocksBySize,1)) {
	    if (blk.active) { blk.detect_shock_points(); }
	}
    }
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        if (blk.active) { blk.convective_flux_phase0(allow_high_order_interpolation, gtl); }
    }
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        if (blk.active) { blk.convective_flux_phase1(allow_high_order_interpolation, gtl); }
    }
    if (GlobalConfig.apply_bcs_in_parallel) {
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
        }
    } else {
        foreach (blk; localFluidBlocksBySize) {
            if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
        }
    }
    if (GlobalConfig.viscous && !GlobalConfig.separate_update_for_viscous_terms) {
        if (GlobalConfig.apply_bcs_in_parallel) {
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) {
                    blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
                    blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
                }
            }
        } else {
            foreach (blk; localFluidBlocksBySize) {
                if (blk.active) {
                    blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
                    blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
                }
            }
        }
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) {
                blk.flow_property_spatial_derivatives(gtl); 
                blk.estimate_turbulence_viscosity();
	    }
	}
	// we exchange boundary data at this point to ensure the
	// ghost cells along block-block boundaries have the most
	// recent mu_t and k_t values.
	exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
	foreach (blk; parallel(localFluidBlocksBySize,1)) {
	    if (blk.active) {         
		blk.viscous_flux();
            }
        }
        if (GlobalConfig.apply_bcs_in_parallel) {
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
            }
        } else {
            foreach (blk; localFluidBlocksBySize) {
                if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
            }
        }
    } // end if viscous
    foreach (i, blk; parallel(localFluidBlocksBySize,1)) {
        if (!blk.active) continue;
        int local_ftl = ftl;
        int local_gtl = gtl;
        bool local_with_k_omega = with_k_omega;
        bool local_with_local_time_stepping = with_local_time_stepping;
        double local_dt_global = SimState.dt_global;
        double local_sim_time = SimState.time;
        foreach (cell; blk.cells) {
            cell.add_inviscid_source_vector(local_gtl, blk.omegaz);
            if (blk.myConfig.viscous && !blk.myConfig.separate_update_for_viscous_terms) {
                cell.add_viscous_source_vector(local_with_k_omega);
            }
            if (blk.myConfig.udf_source_terms) { // [TODO] may want to apply serially
                size_t i_cell = cell.id;
                size_t j_cell = 0;
                size_t k_cell = 0;
                if (blk.grid_type == Grid_t.structured_grid) {
                    auto sblk = cast(SFluidBlock) blk;
                    assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                    auto ijk_indices = sblk.cell_id_to_ijk_indices(cell.id);
                    i_cell = ijk_indices[0];
                    j_cell = ijk_indices[1];
                    k_cell = ijk_indices[2];
                }
                addUDFSourceTermsToCell(blk.myL, cell, local_gtl, 
                                        local_sim_time, blk.myConfig.gmodel,
                                        blk.id, i_cell, j_cell, k_cell);
            }
            cell.time_derivatives(local_gtl, local_ftl, local_with_k_omega);
        }
        if (blk.myConfig.residual_smoothing) { blk.residual_smoothing_dUdt(local_ftl); }
        bool force_euler = false;
        foreach (cell; blk.cells) {
            cell.stage_1_update_for_flow_on_fixed_grid(local_dt_global, force_euler,
                                                       local_with_k_omega, local_with_local_time_stepping);
            cell.decode_conserved(local_gtl, local_ftl+1, blk.omegaz);
        } // end foreach cell
        local_invalid_cell_count[i] = blk.count_invalid_cells(local_gtl, local_ftl+1);
    } // end foreach blk
    //
    int flagTooManyBadCells = 0;
    foreach (i, blk; localFluidBlocksBySize) { // serial loop
        if (local_invalid_cell_count[i] > GlobalConfig.max_invalid_cells) {
            flagTooManyBadCells = 1;
            writefln("Following first-stage gasdynamic update: %d bad cells in block[%d].",
                     local_invalid_cell_count[i], i);
        }
    }
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &flagTooManyBadCells, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
    if (flagTooManyBadCells > 0) {
        throw new FlowSolverException("Too many bad cells; go home.");
    }
    //
    if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight) {
        // Next do solid domain update IMMEDIATELY after at same flow time level
        if (GlobalConfig.in_mpi_context && solidBlocks.length > 0) {
            throw new Error("oops, should not be doing coupling with solid domains in MPI.");
            // [TODO] 2018-01-20 PJ should do something to lift this restriction.
        }
        foreach (sblk; parallel(solidBlocks, 1)) {
            if (!sblk.active) continue;
            sblk.clearSources();
            sblk.computeSpatialDerivatives(ftl);
            sblk.computeFluxes();
        }
        if (GlobalConfig.apply_bcs_in_parallel) {
            foreach (sblk; parallel(solidBlocks, 1)) {
                if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
            }
        } else {
            foreach (sblk; solidBlocks) {
                if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
            }
        }
        // We need to synchronise before updating
        foreach (sblk; parallel(solidBlocks, 1)) {
            foreach (scell; sblk.activeCells) {
                if (GlobalConfig.udfSolidSourceTerms) {
                    addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time);
                }
                scell.timeDerivatives(ftl, GlobalConfig.dimensions);
                scell.stage1Update(SimState.dt_global);
                scell.T = updateTemperature(scell.sp, scell.e[ftl+1]);
            } // end foreach scell
        } // end foreach sblk
    } // end if tight solid domain coupling.
    //
    if (number_of_stages_for_update_scheme(GlobalConfig.gasdynamic_update_scheme) >= 2) {
        // Preparation for second-stage of gas-dynamic update.
        SimState.time = t0 + c2 * SimState.dt_global;
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) {
                blk.clear_fluxes_of_conserved_quantities();
                foreach (cell; blk.cells) cell.clear_source_vector();
            }
        }
        // Second stage of gas-dynamic update.
        ftl = 1;
        // We are relying on exchanging boundary data as a pre-reconstruction activity.
        exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
        if (GlobalConfig.apply_bcs_in_parallel) {
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
            }
        } else {
            foreach (blk; localFluidBlocksBySize) {
                if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
            }
        }
        // Let's set up solid domain bc's also before changing any flow properties.
        if (GlobalConfig.apply_bcs_in_parallel) {
            foreach (sblk; parallel(solidBlocks, 1)) {
                if (sblk.active) { sblk.applyPreSpatialDerivAction(SimState.time, ftl); }
            }
            foreach (sblk; parallel(solidBlocks, 1)) {
                if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
            }
        } else {
            foreach (sblk; solidBlocks) {
                if (sblk.active) { sblk.applyPreSpatialDerivAction(SimState.time, ftl); }
            }
            foreach (sblk; solidBlocks) {
                if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
            }
        }
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) { blk.convective_flux_phase0(allow_high_order_interpolation, gtl); }
        }
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) { blk.convective_flux_phase1(allow_high_order_interpolation, gtl); }
        }
        if (GlobalConfig.apply_bcs_in_parallel) {
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
            }
        } else {
            foreach (blk; localFluidBlocksBySize) {
                if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
            }
        }
        if (GlobalConfig.viscous && !GlobalConfig.separate_update_for_viscous_terms) {
            if (GlobalConfig.apply_bcs_in_parallel) {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) {
                        blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
                        blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
                    }
                }
            } else {
                foreach (blk; localFluidBlocksBySize) {
                    if (blk.active) {
                        blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
                        blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
                    }
                }
            }
	    foreach (blk; parallel(localFluidBlocksBySize,1)) {
		if (blk.active) {
		    blk.flow_property_spatial_derivatives(gtl); 
		    blk.estimate_turbulence_viscosity();
		}
	    }
	    // we exchange boundary data at this point to ensure the
	    // ghost cells along block-block boundaries have the most
	    // recent mu_t and k_t values.
	    exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
	    foreach (blk; parallel(localFluidBlocksBySize,1)) {
		if (blk.active) {         
		    blk.viscous_flux();
		}
	    }
            if (GlobalConfig.apply_bcs_in_parallel) {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
                }
            } else {
                foreach (blk; localFluidBlocksBySize) {
                    if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
                }
            }
        } // end if viscous
        foreach (i, blk; parallel(localFluidBlocksBySize,1)) {
            if (!blk.active) continue;
            int local_ftl = ftl;
            int local_gtl = gtl;
            bool local_with_k_omega = with_k_omega;
            bool local_with_local_time_stepping = with_local_time_stepping;
            double local_dt_global = SimState.dt_global;
            double local_sim_time = SimState.time;
            foreach (cell; blk.cells) {
                cell.add_inviscid_source_vector(local_gtl, blk.omegaz);
                if (blk.myConfig.viscous && !blk.myConfig.separate_update_for_viscous_terms) {
                    cell.add_viscous_source_vector(local_with_k_omega);
                }
                if (blk.myConfig.udf_source_terms) {
                    size_t i_cell = cell.id;
                    size_t j_cell = 0;
                    size_t k_cell = 0;
                    if (blk.grid_type == Grid_t.structured_grid) {
                        auto sblk = cast(SFluidBlock) blk;
                        assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                        auto ijk_indices = sblk.cell_id_to_ijk_indices(cell.id);
                        i_cell = ijk_indices[0];
                        j_cell = ijk_indices[1];
                        k_cell = ijk_indices[2];
                    }
                    addUDFSourceTermsToCell(blk.myL, cell, local_gtl,
                                            local_sim_time, blk.myConfig.gmodel,
                                            blk.id, i_cell, j_cell, k_cell);
                }
                cell.time_derivatives(local_gtl, local_ftl, local_with_k_omega);
            }
            if (blk.myConfig.residual_smoothing) { blk.residual_smoothing_dUdt(local_ftl); }
            foreach (cell; blk.cells) {
                cell.stage_2_update_for_flow_on_fixed_grid(local_dt_global, local_with_k_omega, local_with_local_time_stepping);
                cell.decode_conserved(local_gtl, local_ftl+1, blk.omegaz);
            } // end foreach cell
        local_invalid_cell_count[i] = blk.count_invalid_cells(local_gtl, local_ftl+1);
        } // end foreach blk
        //
        flagTooManyBadCells = 0;
        foreach (i, blk; localFluidBlocksBySize) { // serial loop
            if (local_invalid_cell_count[i] > GlobalConfig.max_invalid_cells) {
                flagTooManyBadCells = 1;
                writefln("Following second-stage gasdynamic update: %d bad cells in block[%d].",
                         local_invalid_cell_count[i], i);
            }
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &flagTooManyBadCells, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (flagTooManyBadCells > 0) {
            throw new FlowSolverException("Too many bad cells; go home.");
        }
        //
        if ( GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight ) {
            // Do solid domain update IMMEDIATELY after at same flow time level
            foreach (sblk; parallel(solidBlocks, 1)) {
                if (!sblk.active) continue;
                sblk.clearSources();
                sblk.computeSpatialDerivatives(ftl);
                sblk.computeFluxes();
            }
            if (GlobalConfig.apply_bcs_in_parallel) {
                foreach (sblk; parallel(solidBlocks, 1)) {
                    if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
                }
            } else {
                foreach (sblk; solidBlocks) {
                    if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
                }
            }
            // We need to synchronise before updating
            foreach (sblk; parallel(solidBlocks, 1)) {
                foreach (scell; sblk.activeCells) {
                    if (GlobalConfig.udfSolidSourceTerms) {
                        addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time);
                    }
                    scell.timeDerivatives(ftl, GlobalConfig.dimensions);
                    scell.stage2Update(SimState.dt_global);
                    scell.T = updateTemperature(scell.sp, scell.e[ftl+1]);
                } // end foreach cell
            } // end foreach blk
        } // end if tight solid domain coupling.
    } // end if number_of_stages_for_update_scheme >= 2 
    //
    if (number_of_stages_for_update_scheme(GlobalConfig.gasdynamic_update_scheme) >= 3) {
        // Preparation for third stage of gasdynamic update.
        SimState.time = t0 + c3 * SimState.dt_global;
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) {
                blk.clear_fluxes_of_conserved_quantities();
                foreach (cell; blk.cells) cell.clear_source_vector();
            }
        }
        // Third stage of gas-dynamic update.
        ftl = 2;
        // We are relying on exchanging boundary data as a pre-reconstruction activity.
        exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
        if (GlobalConfig.apply_bcs_in_parallel) {
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
            }
        } else {
            foreach (blk; localFluidBlocksBySize) {
                if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
            }
        }
        // Let's set up solid domain bc's also before changing any flow properties.
        if (GlobalConfig.apply_bcs_in_parallel) {
            foreach (sblk; parallel(solidBlocks, 1)) {
                if (sblk.active) { sblk.applyPreSpatialDerivAction(SimState.time, ftl); }
            }
            foreach (sblk; parallel(solidBlocks, 1)) {
                if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
            }
        } else {
            foreach (sblk; solidBlocks) {
                if (sblk.active) { sblk.applyPreSpatialDerivAction(SimState.time, ftl); }
            }
            foreach (sblk; solidBlocks) {
                if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
            }
        }
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) { blk.convective_flux_phase0(allow_high_order_interpolation, gtl); }
        }
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) { blk.convective_flux_phase1(allow_high_order_interpolation, gtl); }
        }
        if (GlobalConfig.apply_bcs_in_parallel) {
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
            }
        } else {
            foreach (blk; localFluidBlocksBySize) {
                if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
            }
        }
        if (GlobalConfig.viscous && !GlobalConfig.separate_update_for_viscous_terms) {
            if (GlobalConfig.apply_bcs_in_parallel) {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) {
                        blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
                        blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
                    }
                }
            } else {
                foreach (blk; localFluidBlocksBySize) {
                    if (blk.active) {
                        blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
                        blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
                    }
                }
            }
	    foreach (blk; parallel(localFluidBlocksBySize,1)) {
		if (blk.active) {
		    blk.flow_property_spatial_derivatives(gtl); 
		    blk.estimate_turbulence_viscosity();
		}
	    }
	    // we exchange boundary data at this point to ensure the
	    // ghost cells along block-block boundaries have the most
	    // recent mu_t and k_t values.
	    exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
	    foreach (blk; parallel(localFluidBlocksBySize,1)) {
		if (blk.active) {         
		    blk.viscous_flux();
		}
	    }
            if (GlobalConfig.apply_bcs_in_parallel) {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
                }
            } else {
                foreach (blk; localFluidBlocksBySize) {
                    if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
                }
            }
        } // end if viscous
        foreach (i, blk; parallel(localFluidBlocksBySize,1)) {
            if (!blk.active) continue;
            int local_ftl = ftl;
            int local_gtl = gtl;
            bool local_with_k_omega = with_k_omega;
            bool local_with_local_time_stepping = with_local_time_stepping;
            double local_dt_global = SimState.dt_global;
            double local_sim_time = SimState.time;
            foreach (cell; blk.cells) {
                cell.add_inviscid_source_vector(local_gtl, blk.omegaz);
                if (blk.myConfig.viscous && !blk.myConfig.separate_update_for_viscous_terms) {
                    cell.add_viscous_source_vector(local_with_k_omega);
                }
                if (blk.myConfig.udf_source_terms) {
                    size_t i_cell = cell.id;
                    size_t j_cell = 0;
                    size_t k_cell = 0;
                    if (blk.grid_type == Grid_t.structured_grid) {
                        auto sblk = cast(SFluidBlock) blk;
                        assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                        auto ijk_indices = sblk.cell_id_to_ijk_indices(cell.id);
                        i_cell = ijk_indices[0];
                        j_cell = ijk_indices[1];
                        k_cell = ijk_indices[2];
                    }
                    addUDFSourceTermsToCell(blk.myL, cell, local_gtl,
                                            local_sim_time, blk.myConfig.gmodel,
                                            blk.id, i_cell, j_cell, k_cell);
                }
                cell.time_derivatives(local_gtl, local_ftl, local_with_k_omega);
            }
            if (blk.myConfig.residual_smoothing) { blk.residual_smoothing_dUdt(local_ftl); }
            foreach (cell; blk.cells) {
                cell.stage_3_update_for_flow_on_fixed_grid(local_dt_global, local_with_k_omega, local_with_local_time_stepping);
                cell.decode_conserved(local_gtl, local_ftl+1, blk.omegaz);
            } // end foreach cell
        local_invalid_cell_count[i] = blk.count_invalid_cells(local_gtl, local_ftl+1);
        } // end foreach blk
        //
        flagTooManyBadCells = 0;
        foreach (i, blk; localFluidBlocksBySize) { // serial loop
            if (local_invalid_cell_count[i] > GlobalConfig.max_invalid_cells) {
                flagTooManyBadCells = 1;
                writefln("Following third-stage gasdynamic update: %d bad cells in block[%d].",
                         local_invalid_cell_count[i], i);
            }
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &flagTooManyBadCells, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (flagTooManyBadCells > 0) {
            throw new FlowSolverException("Too many bad cells; go home.");
        }
        //
        if ( GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.tight ) {
            // Do solid domain update IMMEDIATELY after at same flow time level
            foreach (sblk; parallel(solidBlocks, 1)) {
                if (!sblk.active) continue;
                sblk.clearSources();
                sblk.computeSpatialDerivatives(ftl);
                sblk.computeFluxes();
            }
            if (GlobalConfig.apply_bcs_in_parallel) {
                foreach (sblk; parallel(solidBlocks, 1)) {
                    if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
                }
            } else {
                foreach (sblk; solidBlocks) {
                    if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
                }
            }
            // We need to synchronise before updating
            foreach (sblk; parallel(solidBlocks, 1)) {
                foreach (scell; sblk.activeCells) {
                    if (GlobalConfig.udfSolidSourceTerms) {
                        addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time);
                    }
                    scell.timeDerivatives(ftl, GlobalConfig.dimensions);
                    scell.stage3Update(SimState.dt_global);
                    scell.T = updateTemperature(scell.sp, scell.e[ftl+1]);
                } // end foreach cell
            } // end foreach blk
        } // end if tight solid domain coupling.
    } // end if number_of_stages_for_update_scheme >= 3
    //
    // Get the end conserved data into U[0] for next step.
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        if (blk.active) {
            size_t end_indx = final_index_for_update_scheme(GlobalConfig.gasdynamic_update_scheme);
            foreach (cell; blk.cells) { swap(cell.U[0], cell.U[end_indx]); }
        }
    } // end foreach blk
    //
    foreach (sblk; solidBlocks) {
        if (sblk.active) {
            size_t end_indx = final_index_for_update_scheme(GlobalConfig.gasdynamic_update_scheme);
            foreach (scell; sblk.activeCells) { scell.e[0] = scell.e[end_indx]; }
        }
    } // end foreach sblk
    //
    // Finally, update the globally know simulation time for the whole step.
    SimState.time = t0 + SimState.dt_global;
} // end gasdynamic_explicit_increment_with_fixed_grid()

void gasdynamic_explicit_increment_with_moving_grid()
{
    // For moving grid simulations we move the grid on the first predictor step and then
    // leave it fixed in this position for the corrector steps.
    shared double t0 = SimState.time;
    shared bool with_local_time_stepping = GlobalConfig.with_local_time_stepping;
    shared bool with_k_omega = (GlobalConfig.turbulence_model == TurbulenceModel.k_omega) &&
        !GlobalConfig.separate_update_for_k_omega_source;
    shared bool allow_high_order_interpolation = (SimState.time >= GlobalConfig.interpolation_delay);
    // Set the time-step coefficients for the stages of the update scheme.
    shared double c2 = 1.0; // same for 1-stage or 2-stage update
    shared double c3 = 1.0; // ditto
    
    // Preparation for the predictor-stage of inviscid gas-dynamic flow update.
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        if (blk.active) {
            blk.clear_fluxes_of_conserved_quantities();
            foreach (cell; blk.cells) {
                cell.clear_source_vector();
                cell.data_is_bad = false;
            }
        }
    }
    // First-stage of gas-dynamic update.
    shared int ftl = 0; // time-level within the overall convective-update
    shared int gtl = 0; // grid time-level
    // Moving Grid - predict new vertex positions for moving grid              
    foreach (blk; localFluidBlocksBySize) {
        if (!blk.active) continue;
        auto sblk = cast(SFluidBlock) blk;
        assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
        // move vertices
        predict_vertex_positions(sblk, SimState.dt_global, gtl);
        // recalculate cell geometry with new vertex positions @ gtl = 1
        blk.compute_primary_cell_geometric_data(gtl+1);
        blk.compute_least_squares_setup(gtl+1);
        // determine interface velocities using GCL for gtl = 1
        set_gcl_interface_properties(sblk, gtl+1, SimState.dt_global);
    }
    gtl = 1; // update gtl now that grid has moved
    exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
    if (GlobalConfig.apply_bcs_in_parallel) {
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
        }
    } else {
        foreach (blk; localFluidBlocksBySize) {
            if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
        }
    }
    // And we'll do a first pass on solid domain bc's too in case we need up-to-date info.
    foreach (sblk; solidBlocks) {
        if (sblk.active) { sblk.applyPreSpatialDerivAction(SimState.time, ftl); }
    }
    foreach (sblk; solidBlocks) {
        if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
    }
    // We've put this detector step here because it needs the ghost-cell data
    // to be current, as it should be just after a call to apply_convective_bc().
    if (GlobalConfig.flux_calculator == FluxCalculator.adaptive_hanel_ausmdv ||
        GlobalConfig.flux_calculator == FluxCalculator.adaptive_hlle_roe ||
        GlobalConfig.flux_calculator == FluxCalculator.adaptive_efm_ausmdv) {
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) { blk.detect_shock_points(); }
        }
    }
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        if (blk.active) { blk.convective_flux_phase0(allow_high_order_interpolation, gtl); }
    }
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        if (blk.active) { blk.convective_flux_phase1(allow_high_order_interpolation, gtl); }
    }
    if (GlobalConfig.apply_bcs_in_parallel) {
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
        }
    } else {
        foreach (blk; localFluidBlocksBySize) {
            if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
        }
    }
    if (GlobalConfig.viscous && !GlobalConfig.separate_update_for_viscous_terms) {
        if (GlobalConfig.apply_bcs_in_parallel) {
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) {
                    blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
                    blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
                }
            }
        } else {
            foreach (blk; localFluidBlocksBySize) {
                if (blk.active) {
                    blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
                    blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
                }
            }
        }
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) {
                blk.flow_property_spatial_derivatives(gtl); 
                blk.estimate_turbulence_viscosity();
	    }
	}
	// we exchange boundary data at this point to ensure the
	// ghost cells along block-block boundaries have the most
	// recent mu_t and k_t values.
	exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
	foreach (blk; parallel(localFluidBlocksBySize,1)) {
	    if (blk.active) {         
		blk.viscous_flux();
            }
        }
        if (GlobalConfig.apply_bcs_in_parallel) {
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
            }
        } else {
            foreach (blk; localFluidBlocksBySize) {
                if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
            }
        }
    } // end if viscous
    foreach (i, blk; parallel(localFluidBlocksBySize,1)) {
        if (!blk.active) continue;
        int local_ftl = ftl;
        int local_gtl = gtl;
        bool local_with_k_omega = with_k_omega;
        bool local_with_local_time_stepping = with_local_time_stepping;
        double local_dt_global = SimState.dt_global;
        double local_sim_time = SimState.time;
        foreach (cell; blk.cells) {
            cell.add_inviscid_source_vector(local_gtl, blk.omegaz);
            if (blk.myConfig.viscous && !blk.myConfig.separate_update_for_viscous_terms) {
                cell.add_viscous_source_vector(local_with_k_omega);
            }
            if (blk.myConfig.udf_source_terms) { // [TODO] may want to apply serially
                size_t i_cell = cell.id;
                size_t j_cell = 0;
                size_t k_cell = 0;
                if (blk.grid_type == Grid_t.structured_grid) {
                    auto sblk = cast(SFluidBlock) blk;
                    assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                    auto ijk_indices = sblk.cell_id_to_ijk_indices(cell.id);
                    i_cell = ijk_indices[0];
                    j_cell = ijk_indices[1];
                    k_cell = ijk_indices[2];
                }
                addUDFSourceTermsToCell(blk.myL, cell, local_gtl, 
                                        local_sim_time, blk.myConfig.gmodel,
                                        blk.id, i_cell, j_cell, k_cell);
            }
            cell.time_derivatives(local_gtl, local_ftl, local_with_k_omega);
            bool force_euler = false;
            cell.stage_1_update_for_flow_on_moving_grid(local_dt_global, local_with_k_omega, local_with_local_time_stepping);
            cell.decode_conserved(local_gtl, local_ftl+1, blk.omegaz);
        } // end foreach cell
        local_invalid_cell_count[i] = blk.count_invalid_cells(local_gtl, local_ftl+1);
    } // end foreach blk
    //
    int flagTooManyBadCells = 0;
    foreach (i, blk; localFluidBlocksBySize) { // serial loop
        if (local_invalid_cell_count[i] > GlobalConfig.max_invalid_cells) {
            flagTooManyBadCells = 1;
            writefln("Following first-stage gasdynamic update: %d bad cells in block[%d].",
                     local_invalid_cell_count[i], i);
        }
    }
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE, &flagTooManyBadCells, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    }
    if (flagTooManyBadCells > 0) {
        throw new FlowSolverException("Too many bad cells; go home.");
    }
    //
    // Next do solid domain update IMMEDIATELY after at same flow time level
    foreach (sblk; solidBlocks) {
        if (!sblk.active) continue;
        sblk.clearSources();
        sblk.computeSpatialDerivatives(ftl);
        sblk.applyPostFluxAction(SimState.time, ftl);
        sblk.computeFluxes();
        sblk.applyPostFluxAction(SimState.time, ftl);
        foreach (scell; sblk.activeCells) {
            if (GlobalConfig.udfSolidSourceTerms) {
                addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time);
            }
            scell.timeDerivatives(ftl, GlobalConfig.dimensions);
            scell.stage1Update(SimState.dt_global);
            scell.T = updateTemperature(scell.sp, scell.e[ftl+1]);
        } // end foreach scell
    } // end foreach sblk
    /////
    if (number_of_stages_for_update_scheme(GlobalConfig.gasdynamic_update_scheme) == 2) {
        // Preparation for second-stage of gas-dynamic update.
        SimState.time = t0 + c2 * SimState.dt_global;
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) {
                blk.clear_fluxes_of_conserved_quantities();
                foreach (cell; blk.cells) { cell.clear_source_vector(); }
            }
        }
        // Second stage of gas-dynamic update.
        // Moving Grid - update geometry to gtl 2
        foreach (blk; localFluidBlocksBySize) {
            if (blk.active) {
                auto sblk = cast(SFluidBlock) blk;
                assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                // move vertices - this is a formality since pos[2] = pos[1]
                predict_vertex_positions(sblk, SimState.dt_global, gtl);
                // recalculate cell geometry with new vertex positions
                blk.compute_primary_cell_geometric_data(gtl+1);
                blk.compute_least_squares_setup(gtl+1);
                // grid remains at pos[gtl=1], thus let's use old interface velocities
                // thus no need to set_gcl_interface_properties(blk, 2, dt_global);
            }
        }
        ftl = 1;
        gtl = 2;
        // We are relying on exchanging boundary data as a pre-reconstruction activity.
        exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
        if (GlobalConfig.apply_bcs_in_parallel) {
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
            }
        } else {
            foreach (blk; localFluidBlocksBySize) {
                if (blk.active) { blk.applyPreReconAction(SimState.time, gtl, ftl); }
            }
        }
        // Let's set up solid domain bc's also before changing any flow properties.
        foreach (sblk; solidBlocks) {
            if (sblk.active) { sblk.applyPreSpatialDerivAction(SimState.time, ftl); }
        }
        foreach (sblk; solidBlocks) {
            if (sblk.active) { sblk.applyPostFluxAction(SimState.time, ftl); }
        }
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) { blk.convective_flux_phase0(allow_high_order_interpolation, 0); }
            // FIX-ME PJ 2018-07-25 Should this be gtl rather than 0?
        }
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) { blk.convective_flux_phase1(allow_high_order_interpolation, 0); }
            // FIX-ME PJ 2018-07-25 Should this be gtl rather than 0?
        }
        if (GlobalConfig.apply_bcs_in_parallel) {
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
            }
        } else {
            foreach (blk; localFluidBlocksBySize) {
                if (blk.active) { blk.applyPostConvFluxAction(SimState.time, gtl, ftl); }
            }
        }
        if (GlobalConfig.viscous && !GlobalConfig.separate_update_for_viscous_terms) {
            if (GlobalConfig.apply_bcs_in_parallel) {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) {
                        blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
                        blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
                    }
                }
            } else {
                foreach (blk; localFluidBlocksBySize) {
                    if (blk.active) {
                        blk.applyPreSpatialDerivActionAtBndryFaces(SimState.time, gtl, ftl);
                        blk.applyPreSpatialDerivActionAtBndryCells(SimState.time, gtl, ftl);
                    }
                }
            }
	    foreach (blk; parallel(localFluidBlocksBySize,1)) {
		if (blk.active) {
		    blk.flow_property_spatial_derivatives(gtl); 
		    blk.estimate_turbulence_viscosity();
		}
	    }
	    // we exchange boundary data at this point to ensure the
	    // ghost cells along block-block boundaries have the most
	    // recent mu_t and k_t values.
	    exchange_ghost_cell_boundary_data(SimState.time, gtl, ftl);
	    foreach (blk; parallel(localFluidBlocksBySize,1)) {
		if (blk.active) {         
		    blk.viscous_flux();
		}
	    }
            if (GlobalConfig.apply_bcs_in_parallel) {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
                }
            } else {
                foreach (blk; localFluidBlocksBySize) {
                    if (blk.active) { blk.applyPostDiffFluxAction(SimState.time, gtl, ftl); }
                }
            }
        } // end if viscous
        foreach (i, blk; parallel(localFluidBlocksBySize,1)) {
            if (!blk.active) continue;
            int local_ftl = ftl;
            int local_gtl = gtl;
            bool local_with_k_omega = with_k_omega;
            bool local_with_local_time_stepping = with_local_time_stepping;
            double local_dt_global = SimState.dt_global;
            double local_sim_time = SimState.time;
            foreach (cell; blk.cells) {
                cell.add_inviscid_source_vector(local_gtl, blk.omegaz);
                if (blk.myConfig.viscous && !blk.myConfig.separate_update_for_viscous_terms) {
                    cell.add_viscous_source_vector(local_with_k_omega);
                }
                if (blk.myConfig.udf_source_terms) {
                    size_t i_cell = cell.id;
                    size_t j_cell = 0;
                    size_t k_cell = 0;
                    if (blk.grid_type == Grid_t.structured_grid) {
                        auto sblk = cast(SFluidBlock) blk;
                        assert(sblk !is null, "Oops, this should be an SFluidBlock object.");
                        auto ijk_indices = sblk.cell_id_to_ijk_indices(cell.id);
                        i_cell = ijk_indices[0];
                        j_cell = ijk_indices[1];
                        k_cell = ijk_indices[2];
                    }
                    addUDFSourceTermsToCell(blk.myL, cell, local_gtl,
                                            local_sim_time, blk.myConfig.gmodel,
                                            blk.id, i_cell, j_cell, k_cell);
                }
                cell.time_derivatives(local_gtl, local_ftl, local_with_k_omega);
                cell.stage_2_update_for_flow_on_moving_grid(local_dt_global, local_with_k_omega, local_with_local_time_stepping);
                cell.decode_conserved(local_gtl, local_ftl+1, blk.omegaz);
            } // end foreach cell
            local_invalid_cell_count[i] = blk.count_invalid_cells(local_gtl, local_ftl+1);
        } // end foreach blk
        //
        flagTooManyBadCells = 0;
        foreach (i, blk; localFluidBlocksBySize) { // serial loop
            if (local_invalid_cell_count[i] > GlobalConfig.max_invalid_cells) {
                flagTooManyBadCells = 1;
                writefln("Following second-stage gasdynamic update: %d bad cells in block[%d].",
                         local_invalid_cell_count[i], i);
            }
        }
        version(mpi_parallel) {
            MPI_Allreduce(MPI_IN_PLACE, &flagTooManyBadCells, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        }
        if (flagTooManyBadCells > 0) {
            throw new FlowSolverException("Too many bad cells; go home.");
        }
        //
        // Do solid domain update IMMEDIATELY after at same flow time level
        foreach (sblk; solidBlocks) {
            if (!sblk.active) continue;
            sblk.clearSources();
            sblk.computeSpatialDerivatives(ftl);
            sblk.applyPostFluxAction(SimState.time, ftl);
            sblk.computeFluxes();
            sblk.applyPostFluxAction(SimState.time, ftl);
            foreach (scell; sblk.activeCells) {
                if (GlobalConfig.udfSolidSourceTerms) {
                    addUDFSourceTermsToSolidCell(sblk.myL, scell, SimState.time);
                }
                scell.timeDerivatives(ftl, GlobalConfig.dimensions);
                scell.stage2Update(SimState.dt_global);
                scell.T = updateTemperature(scell.sp, scell.e[ftl+1]);
            } // end foreach cell
        } // end foreach blk
    } // end if number_of_stages_for_update_scheme >= 2 
    //
    // Get the end conserved data into U[0] for next step.
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        if (blk.active) {
            size_t end_indx = final_index_for_update_scheme(GlobalConfig.gasdynamic_update_scheme);
            foreach (cell; blk.cells) { swap(cell.U[0], cell.U[end_indx]); }
        }
    }
    foreach (sblk; solidBlocks) {
        if (sblk.active) {
            size_t end_indx = final_index_for_update_scheme(GlobalConfig.gasdynamic_update_scheme);
            foreach (scell; sblk.activeCells) { scell.e[0] = scell.e[end_indx]; } 
        }
    }
    // update the latest grid level to the new step grid level 0
    foreach (blk; localFluidBlocksBySize) {
        if (blk.active) {
            foreach (cell; blk.cells) { cell.copy_grid_level_to_level(gtl, 0); }
        }
    }
    // Finally, update the globally known simulation time for the whole step.
    SimState.time = t0 + SimState.dt_global;
} // end gasdynamic_explicit_increment_with_moving_grid()

void compute_Linf_residuals(ConservedQuantities Linf_residuals)
{
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        blk.compute_Linf_residuals();
    }
    Linf_residuals.copy_values_from(localFluidBlocks[0].Linf_residuals);
    foreach (blk; localFluidBlocksBySize) {
        Linf_residuals.mass = fmax(Linf_residuals.mass, fabs(blk.Linf_residuals.mass));
        Linf_residuals.momentum.set(fmax(Linf_residuals.momentum.x, fabs(blk.Linf_residuals.momentum.x)),
                                    fmax(Linf_residuals.momentum.y, fabs(blk.Linf_residuals.momentum.y)),
                                    fmax(Linf_residuals.momentum.z, fabs(blk.Linf_residuals.momentum.z)));
        Linf_residuals.total_energy = fmax(Linf_residuals.total_energy, fabs(blk.Linf_residuals.total_energy));

    }
} // end compute_Linf_residuals()

void compute_L2_residual(ref number L2_residual)
{
    L2_residual = to!number(0.0);
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        blk.compute_L2_residual();
    }
    foreach (blk; localFluidBlocksBySize) {
        L2_residual += blk.L2_residual;
    }
} // end compute_Linf_residuals()

void compute_mass_balance(ref number mass_balance)
{
    mass_balance = to!number(0.0);
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        blk.compute_mass_balance();
    }
    foreach (blk; localFluidBlocksBySize) {
        mass_balance += blk.mass_balance;
    }
} // end compute_Linf_residuals()


void finalize_simulation()
{
    if (GlobalConfig.verbosity_level > 0 && GlobalConfig.is_master_task) {
        writeln("Finalize the simulation.");
    }
    if (!SimState.output_just_written) {
        write_solution_files();
        if (GlobalConfig.udf_supervisor_file.length > 0) { call_UDF_at_write_to_file(); }
    }
    if (!SimState.history_just_written) { write_history_cells_to_files(SimState.time); }
    GC.collect();
    GC.minimize();
    if (GlobalConfig.verbosity_level > 0  && GlobalConfig.is_master_task) {
        writeln("Step= ", SimState.step, " final-t= ", SimState.time);
    }
} // end finalize_simulation()
