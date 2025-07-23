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
import json_helper;
import fileutil;
import geom;
import geom.misc.kdtree;
import gas;
import globalconfig;
import globaldata;
import flowstate;
import fluidblock;
import sfluidblock;
import ufluidblock;
import fluidblockio_old;
import fluidblockio_new;
import fluidblockio;
import ssolidblock;
import solidprops;
import solidfvinterface;
import solid_full_face_copy;
import solid_gas_full_face_copy;
import bc.ghost_cell_effect.gas_solid_full_face_copy;
import bc;
import user_defined_source_terms;
import solid_udf_source_terms;
import grid_motion;
import grid_motion_udf;
import grid_motion_shock_fitting;
import history;
import loads;
import conservedquantities;
import special_block_init;
import field;
import stats;
version (opencl_gpu_chem) {
    import opencl_gpu_chem;
}
version (cuda_gpu_chem) {
    import cuda_gpu_chem;
}
version(mpi_parallel) {
    import mpi;
}
import simcore_gasdynamic_step;
import simcore_solid_step;
import simcore_exchange;
import simcore_io;
import celldata;
import util.time_utils;
version (catalyst) {
    import CatalystAdaptor;
}

// The shared double[] flavour of GlobalConfig.userPad can give trouble,
// so we need a normal array for the MPI task to work with.
double[] userPad_copy;

//----------------------------------------------------------------------------

int init_simulation(int tindx, int nextLoadsIndx,
                    int maxCPUs, int threadsPerMPITask, string maxWallClock)
// Returns with fail_flag == 0 for a successful initialization.
{
    if (GlobalConfig.verbosity_level > 0 && GlobalConfig.is_master_task) {
        writeln("Begin init_simulation()...");
    }
    if (GlobalConfig.is_master_task) {
        string progressFile = "config/"~GlobalConfig.base_file_name~"-progress.txt";
        string residualsFile = "config/"~GlobalConfig.base_file_name~"-residuals.txt";
        try {
            std.file.write(progressFile, "0\n");
            // Label columns of residuals for reference by gnuplot.
            if (tindx==0)
            std.file.write(residualsFile, "step time wall-clock mass x-mom y-mom"~
                           " z-mom energy L2 mass-balance\n");
        } catch (Exception e) {
            // do nothing
        }
    }
    //
    SimState.maxWallClockSeconds = timeStringToSeconds(maxWallClock);
    SimState.wall_clock_start = Clock.currTime();
    JSONValue config_jsonData = read_config_file();  // most of the configuration is in here
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

    if (GlobalConfig.new_flow_format) {
        switch (GlobalConfig.flow_format) {
            case "eilmer4text" : GlobalConfig.flowFileExt = "zip"; break;
            case "eilmer4binary" : GlobalConfig.flowFileExt = "zip"; break;
            default : throw new Error("Unrecognised flow format");
        }
    } else {
        switch (GlobalConfig.flow_format) {
            case "gziptext" : GlobalConfig.flowFileExt = "gz"; break;
            case "rawbinary" : GlobalConfig.flowFileExt = "bin"; break;
            default : throw new Error("Unrecognised flow format");
        }
    }

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
        // Assign particular fluid (& solid) blocks to this MPI task and keep a record
        // of the MPI rank for all blocks.
        int my_rank = GlobalConfig.mpi_rank_for_local_task;
        GlobalConfig.mpi_rank_for_block.length = GlobalConfig.nFluidBlocks + GlobalConfig.nSolidBlocks ;
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
            if (taskid == my_rank) {
                auto fblk = cast(FluidBlock) globalBlocks[blkid];
                if (fblk) { localFluidBlocks ~= fblk; }
                auto sblk = cast(SSolidBlock) globalBlocks[blkid];
                if (sblk) { localSolidBlocks ~= sblk; }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        if (localFluidBlocks.length == 0) {
            writefln("MPI-task with rank %d has no FluidBlocks. Quitting.", my_rank);
            MPI_Abort(MPI_COMM_WORLD, 2);
        }
    } else {
        // There is only one process and it deals with all blocks.
        foreach (blk; globalBlocks) {
            auto myfblk = cast(FluidBlock) blk;
            if (myfblk) { localFluidBlocks ~= myfblk; }
            auto mysblk = cast(SSolidBlock) blk;
            if (mysblk) { localSolidBlocks ~= mysblk; }
        }
    }
    foreach (blk; localFluidBlocks) { GlobalConfig.localFluidBlockIds ~= blk.id; }
    foreach (blk; localSolidBlocks) { GlobalConfig.localSolidBlockIds ~= blk.id; }
    //
    // Local blocks may be handled with thread-parallelism.
    auto nBlocksInThreadParallel = localFluidBlocks.length + localSolidBlocks.length;
    // There is no need to have more task threads than blocks local to the process.
    int extraThreadsInPool;
    version(mpi_parallel) {
        extraThreadsInPool = min(threadsPerMPITask-1, nBlocksInThreadParallel-1);
    } else {
        extraThreadsInPool = min(maxCPUs-1, nBlocksInThreadParallel-1);
    }
    defaultPoolThreads(extraThreadsInPool); // total = main thread + extra-threads-in-Pool
    version(mpi_parallel) {
        if (GlobalConfig.verbosity_level > 0) {
            debug {
            writeln("MPI-task with rank ", my_rank, " running with ", extraThreadsInPool+1, " threads.");
            //foreach (blk; localFluidBlocks) { writeln("rank=", my_rank, " blk.id=", blk.id); }
            }
        }
    } else {
        if (GlobalConfig.verbosity_level > 0) {
            writeln("Single process running with ", extraThreadsInPool+1, " threads.");
            // Remember the +1 for the main thread.
        }
    }
    // Now that we have finished assigning blocks to MPI tasks
    // and have identified the FluidBlocks in the local context,
    // we may finish the configuration of the blocks.
    // Note that we initialize the gas model, workspace and
    // the grid and flow arrays for blocks that are in the current
    // MPI-task or process, only.
    foreach (myblk; localFluidBlocks) {
        myblk.init_boundary_conditions(config_jsonData["block_" ~ to!string(myblk.id)]);
        myblk.myConfig.init_gas_model_bits();
        myblk.init_workspace();
        myblk.init_lua_globals();
        foreach (bci; myblk.bc) { bci.post_bc_construction(); }
        if (GlobalConfig.user_pad_length > 0) {
            push_array_to_Lua(myblk.myL, GlobalConfig.userPad, "userPad");
            foreach (bci; myblk.bc) {
                if (bci.myL) { push_array_to_Lua(bci.myL, GlobalConfig.userPad, "userPad"); }
            }
        }
        if (GlobalConfig.udf_source_terms) {
            luaL_dofile(myblk.myL, GlobalConfig.udf_source_terms_file.toStringz);
        }
        // After fully constructing the blocks and its boundary conditions,
        // we can optionally print their representation for checking.
        if (GlobalConfig.verbosity_level > 1) {
            writeln("  Block[", myblk.id, "]: ", myblk);
        }
    }
    // When initializing the grid for each block, there is a possibility that
    // an exception will be thrown if there are negative or zero cell volumes.
    // Having messed up grids is a fairly common occurrence, so we should be ready.
    bool any_block_fail = 0;
    foreach (myblk; parallel(localFluidBlocks,1)) {
        try {
            if (GlobalConfig.grid_motion != GridMotion.none) {
                myblk.init_grid_and_flow_arrays(make_file_name!"grid"(job_name, myblk.id, SimState.current_tindx,
                                                                      GlobalConfig.gridFileExt));
            } else {
                // Assume there is only a single, static grid stored at tindx=0
                myblk.init_grid_and_flow_arrays(make_file_name!"grid"(job_name, myblk.id, 0, GlobalConfig.gridFileExt));
            }
            myblk.compute_primary_cell_geometric_data(0);
            myblk.add_IO();
        } catch (Exception e) {
            writefln("Block[%d] failed to initialize geometry, msg=%s", myblk.id, e.msg);
            any_block_fail = true;
        }
    }
    version(mpi_parallel) {
        int myFlag = to!int(any_block_fail);
        MPI_Allreduce(MPI_IN_PLACE, &myFlag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        any_block_fail = to!bool(myFlag);
    }
    if (any_block_fail) { return -1; }
    //
    // Note that the global id is across all processes, not just the local collection of blocks.
    foreach (i, blk; globalBlocks) {
        auto fluidblk = cast(FluidBlock) blk;
        if (fluidblk) {
            if ( i == 0 ) {
                fluidblk.globalCellIdStart = 0;
            } else {
                auto prev_fluidblk = cast(FluidBlock) globalBlocks[i-1];
                fluidblk.globalCellIdStart = prev_fluidblk.globalCellIdStart + prev_fluidblk.ncells_expected;
            }
        }
    }
    shared double[] time_array;
    time_array.length = localFluidBlocks.length;
    foreach (i, myblk; parallel(localFluidBlocks,1)) {
        myblk.identify_reaction_zones(0);
        myblk.identify_turbulent_zones(0);
        myblk.identify_suppress_reconstruction_zones();
        myblk.identify_suppress_viscous_stresses_zones();
        if (GlobalConfig.new_flow_format) {
            time_array[i] = myblk.read_solution(make_file_name("CellData", CellData.tag, job_name, myblk.id, SimState.current_tindx,
                                                            GlobalConfig.flowFileExt), false);
        } else {
            time_array[i] = myblk.read_solution(make_file_name!"flow"(job_name, myblk.id, SimState.current_tindx,
                                                            GlobalConfig.flowFileExt), false);
        }

        if (myblk.myConfig.verbosity_level >= 2) { writefln("Cold start cells in block %d", myblk.id); }
        // Note that, even for grid_motion==none simulations, we use the grid velocities for setting
        // the gas velocities at boundary faces.  These will need to be set to zero for correct viscous simulation.
        foreach (iface; myblk.faces) { iface.gvel.clear(); }
        foreach (cell; myblk.cells) {
            cell.encode_conserved(0, 0, myblk.omegaz);
            // Even though the following call appears redundant at this point,
            // fills in some gas properties such as Prandtl number that is
            // needed for both the cfl_check and the BaldwinLomax turbulence model.
            if (0 != cell.decode_conserved(0, 0, myblk.omegaz)) {
                writefln("Block[%d] Bad cell decode_conserved while initializing flow.", myblk.id);
                any_block_fail = true;
            }
        }
        myblk.set_cell_dt_chem(-1.0);
    }
    version(mpi_parallel) {
        myFlag = to!int(any_block_fail);
        MPI_Allreduce(MPI_IN_PLACE, &myFlag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        any_block_fail = to!bool(myFlag);
    }
    if (any_block_fail) { return -1; }
    //
    SimState.time = time_array[0]; // Pick one; they should all be the same.
    //
    version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }
    //
    // Now that the cells for all gas blocks have been initialized,
    // we can sift through the boundary condition effects and
    // set up the ghost-cell mapping for the appropriate boundaries.
    if (GlobalConfig.verbosity_level >= 2) { writeln("Prepare exchange of boundary information."); }
    //
    // First, check that full-face exchange boundaries are paired correctly.
    foreach (my_blk; localFluidBlocks) {
        foreach (j, bc; my_blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto my_gce = cast(GhostCellFullFaceCopy)gce;
                if (my_gce) {
                    // The local block thinks that it has an exchange boundary with another block,
                    // so we need to check the ghost-cell effects of the other block's face to see
                    // that it points back to the local block face.

                    size_t other_blk_id   = my_gce.neighbourBlock.id;
                    size_t other_blk_bndy = my_gce.neighbourFace;
                    bool ok = check_full_face_copy(my_blk.id, j, other_blk_id, other_blk_bndy, config_jsonData);
                    if (!ok) {
                        string msg = format("FullFaceCopy for local blk_id=%d face=%d", my_blk.id, j);
                        msg ~= format(" is not correctly paired with other block id=%d face=%d.",
                                      other_blk_id, other_blk_bndy);
                        writeln(msg);
                        any_block_fail = true;
                    }
                } // end if (my_copy_gce)
            } // end foreach (gce;
        } // end foreach (j, bc;
    } // end foreach (my_blk;
    version(mpi_parallel) {
        myFlag = to!int(any_block_fail);
        MPI_Allreduce(MPI_IN_PLACE, &myFlag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        any_block_fail = to!bool(myFlag);
    }
    if (any_block_fail) { return -1; }
    //
    // Serial loops follow because the cell-mapping function searches across
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
                if (mygce && (mygce.check_cell_mapping() != 0)) { any_block_fail = true; }
            }
        }
    }
    version(mpi_parallel) {
        myFlag = to!int(any_block_fail);
        MPI_Allreduce(MPI_IN_PLACE, &myFlag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        any_block_fail = to!bool(myFlag);
    }
    if (any_block_fail) { return -1; }
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
    // Earlier, compute_primary_cell_geometric_data filled in approximate values for the ghost
    // cells' volume, position, area etc. But for the internal ghost cells we really want the
    // exact numbers. These are copied in here from the real cells in their corresponding
    // neighbour blocks.
    exchange_ghost_cell_geometry_data();
    //
    // These fixed up ghost cells might have counterparts on the actual domain boundaries.
    // For viscous walls in particular, it is important that these ghost cells mirror their
    // interior counterparts precisely. So, here we do a second pass on all of the non-shared
    // boundaries, in case they need to be updated.
    foreach (myblk; localFluidBlocks) { myblk.update_nonshared_ghost_cell_positions(0); }
    //
    // Now that we know the ghost-cell locations, we can set up the least-squares subproblems for
    // 1. reconstruction prior to convective flux calculation for the unstructured-grid blocks
    // 2. calculation of flow gradients for the viscous fluxes with least-squares gradients.
    foreach (myblk; localFluidBlocks) { myblk.compute_least_squares_setup(0); }
    //
    // We can also setup the precomputed stencil data for the structured interpolation
    // In an unstructured block, the following call returns without doing anything
    foreach (myblk; localFluidBlocks) { myblk.precompute_stencil_data(0); }

    version (gpu_chem) {
        initGPUChem();
    }
    // solid blocks assigned prior to this
    foreach (ref mySolidBlk; localSolidBlocks) {
        mySolidBlk.assembleArrays();
        mySolidBlk.bindFacesAndVerticesToCells();
        writeln("mySolidBlk= ", mySolidBlk);
        mySolidBlk.readGrid(make_file_name!"solid-grid"(job_name, mySolidBlk.id, 0, "gz")); // tindx==0 fixed grid
        mySolidBlk.readSolution(make_file_name!"solid"(job_name, mySolidBlk.id, tindx, "gz"));
        mySolidBlk.computePrimaryCellGeometricData();
        mySolidBlk.assignCellLocationsForDerivCalc();
    }
    // setup communication across blocks
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(SolidGCE_SolidGhostCellFullFaceCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping_phase0(); }
            }
        }
    }
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(SolidGCE_SolidGhostCellFullFaceCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping_phase1(); }
            }
        }
    }
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(SolidGCE_SolidGhostCellFullFaceCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping_phase2(); }
            }
        }
    }
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(SolidGCE_SolidGhostCellFullFaceCopy)gce;
                if (mygce) { mygce.exchange_solid_data_phase0(); }
            }
        }
    }
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(SolidGCE_SolidGhostCellFullFaceCopy)gce;
                if (mygce) { mygce.exchange_solid_data_phase1(); }
            }
        }
    }
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(SolidGCE_SolidGhostCellFullFaceCopy)gce;
                if (mygce) { mygce.exchange_solid_data_phase2(); }
            }
        }
    }
    // Now that we know the ghost-cell locations, we can set up the least-squares subproblems for
    // calculation of temperature gradients for the solid solver with least-squares gradients.
    foreach (ref mySolidBlk; localSolidBlocks) {
        mySolidBlk.setupSpatialDerivativeCalc();
    }
    if (localSolidBlocks.length > 0) {
        initPropertiesAtSolidInterfaces(localSolidBlocks);
    }
    if (GlobalConfig.coupling_with_solid_domains == SolidDomainCoupling.loose) {
        // currently we use an explicit-explicit loose coupling approach
        // so there is no need to initialise anything at this stage.
        //initSolidLooseCouplingUpdate();
    }
    //
    // All cells are in place, so now we can initialise any history cell files.
    if (GlobalConfig.is_master_task) { ensure_directory_is_present(histDir); }
    version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }
    init_history_cell_files();
    //
    // Create the loads directory, maybe.
    if (GlobalConfig.write_loads && (SimState.current_loads_tindx == 0)) {
        if (GlobalConfig.is_master_task) {
            ensure_directory_is_present("loads");
            init_loads_times_file();
        }
    }
    version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }
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
                fba.read_rails_file(format("config/fba-%04d.rails", i));
                fba.read_velocity_weights(format("config/fba-%04d.weights", i));
                fba.read_shockfitting_inflow("config/"~GlobalConfig.base_file_name~".config", i);
            }
        }
    }
    // Finally when both gas AND solid domains are setup..
    // Look for a solid-adjacent bc, if there is one,
    // then we can set up the cells and interfaces that
    // internal to the bc. They are only known after this point.
    foreach (myblk; localFluidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellGasSolidFullFaceCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping_phase0(); }
            }
        }
    }
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(GhostCellSolidGasFullFaceCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping_phase0(); }
            }
        }
    }
    foreach (myblk; localFluidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellGasSolidFullFaceCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping_phase1(); }
            }
        }
    }
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(GhostCellSolidGasFullFaceCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping_phase1(); }
            }
        }
    }
    foreach (myblk; localFluidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellGasSolidFullFaceCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping_phase2(); }
            }
        }
    }
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(GhostCellSolidGasFullFaceCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping_phase2(); }
            }
        }
    }
    foreach (myblk; localFluidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellGasSolidFullFaceCopy)gce;
                if (mygce) { mygce.exchange_solidstate_phase0(); }
            }
        }
    }
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(GhostCellSolidGasFullFaceCopy)gce;
                if (mygce) { mygce.exchange_fluidstate_phase0(); }
            }
        }
    }
    foreach (myblk; localFluidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellGasSolidFullFaceCopy)gce;
                if (mygce) { mygce.exchange_fluidstate_phase1(); }
            }
        }
    }
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(GhostCellSolidGasFullFaceCopy)gce;
                if (mygce) { mygce.exchange_solidstate_phase1(); }
            }
        }
    }
    foreach (myblk; localFluidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellGasSolidFullFaceCopy)gce;
                if (mygce) { mygce.exchange_solidstate_phase2(); }
            }
        }
    }
    foreach (myblk; localSolidBlocks) {
        foreach (bc; myblk.bc) {
            foreach (gce; bc.preSpatialDerivActionAtBndryCells) {
                auto mygce = cast(GhostCellSolidGasFullFaceCopy)gce;
                if (mygce) { mygce.exchange_fluidstate_phase2(); }
            }
        }
    }
    // For the MLP limiter (on unstructured grids only), we need access to the
    // gradients stored in the cloud of cells surrounding a vertex.
    if ((GlobalConfig.interpolation_order > 1) &&
        (GlobalConfig.unstructured_limiter == UnstructuredLimiter.venkat_mlp)) {
        foreach (myblk; localFluidBlocks) {
            auto ublock = cast(UFluidBlock) myblk;
            if (ublock) { ublock.build_cloud_of_cell_references_at_each_vertex(); }
        }
    }
    //
    // We can apply a special initialisation to the flow field, if requested.
    // This will take viscous boundary conditions and diffuse them into the
    // nearby domain.
    if (GlobalConfig.diffuseWallBCsOnInit && (SimState.time == 0.0)) {
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
    SimState.catalyst_just_called = true;
    // When starting a new calculation,
    // set the global time step to the initial value.
    // If we have elected to run with a variable time step,
    // this value may get revised on the very first step.
    SimState.dt_global = GlobalConfig.dt_init;
    //
    // Fill in some more global information for the master interpreter.
    auto L = GlobalConfig.master_lua_State;
    lua_pushboolean(L, GlobalConfig.in_mpi_context);
    lua_setglobal(L, "in_mpi_context");
    lua_pushnumber(L, GlobalConfig.mpi_size);
    lua_setglobal(L, "mpi_size");
    lua_pushnumber(L, GlobalConfig.mpi_rank_for_local_task);
    lua_setglobal(L, "mpi_rank_for_local_task");
    lua_pushboolean(L, GlobalConfig.is_master_task);
    lua_setglobal(L, "is_master_task");
    push_array_to_Lua(L, GlobalConfig.localFluidBlockIds, "localFluidBlockIds");
    if (GlobalConfig.user_pad_length > 0) {
        // At this point, userPad has been initialized with values
        // from the job.config file, even if there are no supervisory functions.
        push_array_to_Lua(L, GlobalConfig.userPad, "userPad");
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
    if (GlobalConfig.turb_model.needs_dwall) compute_wall_distances();
    //

    if (GlobalConfig.solve_electric_field){
        eField = new ElectricField(localFluidBlocks, GlobalConfig.field_conductivity_model);
    }
    //
    if (GlobalConfig.dt_catalyst > 0.0) {
        version(catalyst) {
            do_catalyst_initialization();
            InitializeCatalystData(&catalyst_data, localFluidBlocks, GlobalConfig.dimensions);
        }
    }
    // if stats TODO
    if (GlobalConfig.save_flowstats){
        flowstats = new FlowStats(GlobalConfig.turb_model.nturb, GlobalConfig.n_species, 20);
        foreach(blk; localFluidBlocks) {
            flowstats.init_map_arrays_for_block(blk.celldata.positions, blk.celldata.cell_to_bin_map,
                                                blk.celldata.bin_weights);
        }
        flowstats.sum_weights_over_mpi(GlobalConfig.is_master_task);
    }
    // Keep our memory foot-print small.
    GC.collect();
    GC.minimize();
    //
    version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }

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
        writefln("simcore memory used: %.0f MB, unused: %.0f MB, total: %.0f MB (%.0f-%.0f MB per task)",
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

    return 0; // Successfully initialized simulation.
} // end init_simulation()


void march_over_blocks()
{
    // Notes:
    // (1) We assume that the full collection of blocks was assembled
    //     as a single structured block with indices 0<=i<nib, 0<=j<njb, 0<=k<nkb
    //     and that blocks 0<=i<nib for a particular j,k are on one MPI task.
    // (2) There may be more than on j,k slice on any MPI task.
    // (3) We need to be careful that the local MPI task does not fiddle with
    //     blocks assigned to another MPI task.
    if (GlobalConfig.verbosity_level > 0 && GlobalConfig.is_master_task) {
        writeln("March over blocks.");
    }
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
    bool[][][] gasBlockIsLocal;
    gasBlockArray.length = nib;
    gasBlockIsLocal.length = nib;
    foreach (i; 0 .. nib) {
        gasBlockArray[i].length = njb;
        gasBlockIsLocal[i].length = njb;
        foreach (j; 0 .. njb) {
            gasBlockArray[i][j].length = nkb;
            gasBlockIsLocal[i][j].length = nkb;
            foreach (k; 0 .. nkb) {
                int gid = k + nkb*(j + njb*i);
                auto fluidblk = cast(FluidBlock) globalBlocks[gid];
                if (fluidblk) {
                    gasBlockArray[i][j][k] = fluidblk;
                    gasBlockIsLocal[i][j][k] = canFind(GlobalConfig.localFluidBlockIds, fluidblk.id);
                } else {
                    string errMsg = text("march_over_blocks(): expected a FluidBlock with id=", gid);
                    throw new FlowSolverException(errMsg);
                }
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
                if (!gasBlockIsLocal[1][j][k]) { continue; } // skip over non-local blocks
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
                if (!gasBlockIsLocal[i][j][k]) { continue; } // skip over non-local blocks.
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
                    if (!gasBlockIsLocal[i][j][k]) { continue; } // skip over non-local blocks.
                    auto blk = gasBlockArray[i][j][k]; // our newly active block
                    // Get upstream flow data into ghost cells
                    blk.applyPreReconAction(SimState.time, 0, 0);
                    // and propagate it across the domain.
                    blk.propagate_inflow_data_west_to_east();
                }
            }
        }
        if (GlobalConfig.verbosity_level > 0 && GlobalConfig.is_master_task) {
            writeln("march over blocks i=", i);
        }
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
    auto Linf_residuals = new_ConservedQuantities(GlobalConfig.cqi.n);
    version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }
    if (GlobalConfig.verbosity_level > 0 && GlobalConfig.is_master_task) {
        writeln("Integrate in time.");
        stdout.flush();
    }
    SimState.target_time = (GlobalConfig.block_marching) ? target_time_as_requested : GlobalConfig.max_time;
    // The next time for output...
    SimState.t_plot = SimState.time + GlobalConfig.dt_plot;
    SimState.t_history = SimState.time + GlobalConfig.dt_history;
    SimState.t_loads = SimState.time + GlobalConfig.dt_loads;
    SimState.t_catalyst = SimState.time + GlobalConfig.dt_catalyst;
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
            check_run_time_configuration(target_time_as_requested);
            if (GlobalConfig.grid_motion != GridMotion.none) { synchronize_corner_coords_for_all_blocks(); }
            //
            // 1.0 Maintain a stable time step size, and other maintenance, as required.
            // The user may also have some start-of-time-step maintenance to do
            // via their Lua script file.  Let them have first go.
            if (GlobalConfig.udf_supervisor_file.length > 0) { call_UDF_at_timestep_start(); }
            determine_time_step_size();
            if (GlobalConfig.divergence_cleaning) { update_ch_for_divergence_cleaning(); }
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
            if (GlobalConfig.solve_electric_field){
                eField.solve_efield(localFluidBlocks);
                eField.compute_boundary_current(localFluidBlocks);
            }

            // 3.0 Update the time record and (occasionally) print status.
            SimState.step = SimState.step + 1;
            if (GlobalConfig.is_master_task) {
                string progressFile = "config/"~GlobalConfig.base_file_name~"-progress.txt";
                try {
                    std.file.write(progressFile, format("%d\n", SimState.step));
                } catch (Exception e) {
                    // do nothing
                }
            }
            //
            SimState.output_just_written = false;
            SimState.history_just_written = false;
            SimState.loads_just_written = false;
            SimState.catalyst_just_called = false;
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
                version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }
                if (GlobalConfig.report_residuals) {
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
                        string residualsFile = "config/"~GlobalConfig.base_file_name~"-residuals.txt";
                        string txt = format("%7d %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e\n",
                                            SimState.step, SimState.time, wall_clock_elapsed,
                                            Linf_residuals[cqi.mass].re,
                                            Linf_residuals[cqi.xMom].re,
                                            Linf_residuals[cqi.yMom].re,
                                            (cqi.threeD) ? Linf_residuals[cqi.zMom].re : 0.0,
                                            Linf_residuals[cqi.totEnergy].re,
                                            fabs(L2_residual.re), fabs(mass_balance.re));
                        std.file.append(residualsFile, txt);
                    }
                } // end if report_residuals
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

            // update any auxiliary values
            if (GlobalConfig.new_flow_format) {
                foreach (myblk; parallel(localFluidBlocksBySize,1)) {
                    myblk.update_aux(SimState.dt_global, SimState.time, SimState.step);
                }
            }

            //
            // 4.0 (Occasionally) Write out an intermediate solution
            if ( SimState.step == GlobalConfig.write_flow_solution_at_step ) {
                write_solution_files();
                GC.collect();
                GC.minimize();
            }
            if ((SimState.time >= SimState.t_plot) && !SimState.output_just_written) {
                write_solution_files();
                if (GlobalConfig.udf_supervisor_file.length > 0) { call_UDF_at_write_to_file(); }
                SimState.output_just_written = true;
                SimState.t_plot = SimState.t_plot + GlobalConfig.dt_plot;
                GC.collect();
                GC.minimize();
            }
            //
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
            //
            // 4.2 (Occasionally) Write out the cell history data and loads on boundary groups data
            if ((SimState.time >= SimState.t_history) && !SimState.history_just_written) {
                write_history_cells_to_files(SimState.time);
                SimState.history_just_written = true;
                SimState.t_history = SimState.t_history + GlobalConfig.dt_history;
                GC.collect();
                GC.minimize();
            }
            if (GlobalConfig.write_loads &&
                ( ((SimState.time >= SimState.t_loads) && !SimState.loads_just_written) ||
                  SimState.step == GlobalConfig.write_loads_at_step )) {
                if (GlobalConfig.is_master_task) {
                    init_current_loads_tindx_dir(SimState.current_loads_tindx);
                }
                version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }
                wait_for_current_tindx_dir(SimState.current_loads_tindx);
                write_boundary_loads_to_file(SimState.time, SimState.current_loads_tindx);
                if (GlobalConfig.is_master_task) {
                    update_loads_times_file(SimState.time, SimState.current_loads_tindx);
                }
                SimState.loads_just_written = true;
                SimState.current_loads_tindx = SimState.current_loads_tindx + 1;
                SimState.t_loads = SimState.t_loads + GlobalConfig.dt_loads;
                GC.collect();
                GC.minimize();
            }
            //
            // 4.3 Increment the DFT in each cell
            if ((SimState.step % GlobalConfig.DFT_step_interval == 0) && GlobalConfig.do_temporal_DFT) {
                foreach (blk; localFluidBlocks) {
                    blk.increment_DFT(SimState.step / GlobalConfig.DFT_step_interval - 1);
                }
            }
            // 4.4 Possibly execute a catalyst pipeline, if configured for it.
            version(catalyst){
                if ((GlobalConfig.dt_catalyst > 0.0) && ((SimState.time >= SimState.t_catalyst) && !SimState.catalyst_just_called)){
                    if (GlobalConfig.is_master_task) {
                        writefln("    Execute catalyst pipeline: (step %d, time %03.03e)", SimState.step, SimState.time);
                    }
                    UpdateCatalystFieldData(&catalyst_data, localFluidBlocks);
                    do_catalyt_execute(SimState.step, SimState.time, &catalyst_data);
                    SimState.catalyst_just_called = true;
                    SimState.t_catalyst = SimState.t_catalyst + GlobalConfig.dt_catalyst;
                }
            }
            // 4.5 FlowStats if they are asked for TODO
            if (GlobalConfig.save_flowstats) {
                foreach(blk; localFluidBlocks) {
                    flowstats.update(SimState.time, blk.celldata.cell_to_bin_map,
                                     blk.celldata.bin_weights, blk.celldata.flowstates);
                }

                flowstats.increment_time_index();
                if (flowstats.tidx==flowstats.NTIME) {
                    if (GlobalConfig.is_master_task) {
                        writef(" -- Flushing stats buffer to file...");
                    }
                    flowstats.dump_stats_to_file_binary(GlobalConfig.is_master_task, GlobalConfig.turb_model, GlobalConfig.gmodel_master);
                    flowstats.reset_buffers();
                    flowstats.reset_time_index();
                }
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
    //
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

//-------------------------------------------------------------------------

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
        foreach (blk; globalBlocks) {
            auto sblk = cast(SFluidBlock) blk;
            if (!sblk) { continue; }
            if (canFind(GlobalConfig.localFluidBlockIds, sblk.id)) {
                // We can see this inside this block to get valid coordinate values.
                sblk.copy_current_corner_coords();
            } else {
                // Cannot see this block so fill in invalid coordinate values.
                sblk.set_current_corner_coords_to_infinity();
            }
        }
        // Now, propagate the valid coordinates across all tasks.
        foreach (blk; globalBlocks) {
            auto sblk = cast(SFluidBlock) blk;
            if (sblk) {
                MPI_Allreduce(MPI_IN_PLACE, sblk.corner_coords.ptr, sblk.corner_coords.length,
                              MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            }
        }
    } else {
        // In shared-memory, we can see all blocks.
        foreach (blk; globalBlocks) {
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
            fill_array_from_Lua(L, GlobalConfig.userPad, "userPad");
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
            fill_array_from_Lua(L, GlobalConfig.userPad, "userPad");
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
            fill_array_from_Lua(L, GlobalConfig.userPad, "userPad");
        }
    }
    lua_settop(L, 0); // clear stack
    if (GlobalConfig.user_pad_length > 0) {
        broadcast_master_userPad();
        copy_userPad_into_block_interpreters();
    }
} // end call_UDF_at_write_to_file()

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

void set_mu_and_k()
{
    version(turbulence){
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            if (blk.active) {
                blk.flow_property_spatial_derivatives(0);
                blk.estimate_turbulence_viscosity();
            }
        }
    }
} // end set_mu_and_k()

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
                if (GlobalConfig.with_local_time_stepping) {
                    foreach (cell; blk.cells) { cell.thermochemical_increment(cell.dt_local); }
                } else {
                    foreach (cell; blk.cells) { cell.thermochemical_increment(dt); }
                }
            }
        }
    }
} // end chemistry_half_step()


void compute_Linf_residuals(ref ConservedQuantities Linf_residuals)
{
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        blk.compute_Linf_residuals();
    }
    Linf_residuals.copy_values_from(localFluidBlocks[0].Linf_residuals);
    auto cqi = GlobalConfig.cqi;
    foreach (blk; localFluidBlocksBySize) {
        foreach (i; 0 .. cqi.n) {
            Linf_residuals[i] = fmax(Linf_residuals[i], fabs(blk.Linf_residuals[i]));
        }
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
} // end compute_L2_residuals()

void compute_mass_balance(ref number mass_balance)
{
    mass_balance = to!number(0.0);
    foreach (blk; parallel(localFluidBlocksBySize,1)) {
        blk.compute_mass_balance();
    }
    foreach (blk; localFluidBlocksBySize) {
        mass_balance += blk.mass_balance;
    }
} // end compute_mass_balance()


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
    if (!SimState.loads_just_written) {
        if (GlobalConfig.is_master_task) {
            init_current_loads_tindx_dir(SimState.current_loads_tindx);
        }
        version(mpi_parallel) { MPI_Barrier(MPI_COMM_WORLD); }
        wait_for_current_tindx_dir(SimState.current_loads_tindx);
        write_boundary_loads_to_file(SimState.time, SimState.current_loads_tindx);
        if (GlobalConfig.is_master_task) {
            update_loads_times_file(SimState.time, SimState.current_loads_tindx);
        }
    }
    if (!GlobalConfig.new_flow_format && GlobalConfig.do_temporal_DFT) {
        write_DFT_files();
    }
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
        string progressFile = "config/"~GlobalConfig.base_file_name~"-progress.txt";
        try {
            std.file.write(progressFile, "done\n");
        } catch (Exception e) {
            // do nothing
        }
    }
    if (GlobalConfig.dt_catalyst > 0.0) {
        version(catalyst) {
            FinalizeCatalystData(&catalyst_data);
            do_catalyt_finalization();
        }
    }
} // end finalize_simulation()

void compute_wall_distances() {
    /*
        Compute the distance from each cell to the nearest viscous wall,
        for any turbulence models that need it.

        @author: Nick Gibbons
    */
    // First count many viscous wall faces are in our local blocks
    int nfaces = 0;
    foreach (blk; localFluidBlocksBySize) {
        foreach(bc; blk.bc) {
            if (bc.is_wall_with_viscous_effects) nfaces+=bc.faces.length;
        }
    }
    //
    // Now pack their centroid positions and normal vectors and lengths into a special buffer
    double[] facepos;
    size_t ii=0;
    facepos.length = nfaces*7;
    int this_rank = GlobalConfig.mpi_rank_for_local_task;
    //
    foreach(blk; localFluidBlocksBySize) {
        foreach(bc; blk.bc) {
            if (bc.is_wall_with_viscous_effects){
                foreach(face; bc.faces){
                    version(complex_numbers){
                    facepos[ii++] = face.pos.x.re;
                    facepos[ii++] = face.pos.y.re;
                    facepos[ii++] = face.pos.z.re;
                    facepos[ii++] = face.n.x.re;
                    facepos[ii++] = face.n.y.re;
                    facepos[ii++] = face.n.z.re;
                    facepos[ii++] = face.length.re;
                    } else {
                    facepos[ii++] = face.pos.x;
                    facepos[ii++] = face.pos.y;
                    facepos[ii++] = face.pos.z;
                    facepos[ii++] = face.n.x;
                    facepos[ii++] = face.n.y;
                    facepos[ii++] = face.n.z;
                    facepos[ii++] = face.length;
                    }
                }
            }
        }
    }
    // Now we need to accumulate faces from across the other MPI tasks
    version(mpi_parallel) {
        int my_rank = GlobalConfig.mpi_rank_for_local_task;
        int mpi_worldsize = GlobalConfig.mpi_size;
        //
        int globalsize = to!int(facepos.length);
        MPI_Allreduce(MPI_IN_PLACE,&globalsize,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        //
        double[] globalfacepos;
        globalfacepos.length = globalsize;
        //
        int ngathered = 0;
        double[] taskbuffer;
        int ntask;
        foreach (task; 0 .. mpi_worldsize){
            if (my_rank == task) ntask = to!int(facepos.length);
            MPI_Bcast(&ntask, 1, MPI_INT, task, MPI_COMM_WORLD);
            //
            taskbuffer.length = ntask;
            if (my_rank == task) foreach(i; 0 .. ntask) taskbuffer[i] = facepos[i];
            MPI_Bcast(taskbuffer.ptr, ntask, MPI_DOUBLE, task, MPI_COMM_WORLD);
            //
            foreach(i; 0 .. ntask) globalfacepos[ngathered+i] = taskbuffer[i];
            ngathered += ntask;
            taskbuffer.length=0;
        }
        // Now clean up by copying the globaldata back into facepos
        facepos.length = globalfacepos.length;
        foreach(i; 0 .. globalfacepos.length)  facepos[i] = globalfacepos[i];
        nfaces = to!int(facepos.length)/7;
    } // end version(mpi_parallel)
    //
    if (nfaces == 0) {
        throw new Exception("Turbulence model requires wall distance, but no walls found!");
    }
    // These centroids need to be assembled into a kdtree
    size_t totalfaces = nfaces;
    Node[] nodes;
    foreach(i; 0 .. nfaces) {
        Node node = {[facepos[7*i+0], facepos[7*i+1], facepos[7*i+2]]};
        node.cellid = i;
        nodes ~= node;
    }
    auto root = makeTree(nodes);
    //
    // Now loop over the nodes in each of our local blocks and set dwall
    foreach(blk; localFluidBlocksBySize) {
        foreach(i, cell; blk.cells){
            Node cellnode = {[cell.pos[0].x.re, cell.pos[0].y.re, cell.pos[0].z.re]};
            const(Node)* found = null;
            double bestDist = 0.0;
            size_t nVisited = 0;
            root.fast_nearest(cellnode, 0, found, bestDist, nVisited);
            //double dist = sqrt(bestDist); // This caused problems on distorted grid cells.

            // Unpack face info into a point and a normal, and compute the normal distance to the resulting plane
            double[3] n = [facepos[7*found.cellid+3],
                           facepos[7*found.cellid+4],
                           facepos[7*found.cellid+5]];
            double[3] px = [cell.pos[0].x.re - found.x[0],
                            cell.pos[0].y.re - found.x[1],
                            cell.pos[0].z.re - found.x[2]];
            double dnormal = fabs(px[0]*n[0] + px[1]*n[1] + px[2]*n[2]);
            // We also care about the parallel direction, since we could be far from a wall but still close to its
            // infinite plane. See derivation 28/10/21 by NNG for the new correction factor offwalldistance.
            double celllength = facepos[7*found.cellid+6];
            double dpx_squared = px[0]*px[0]  + px[1]*px[1]  + px[2]*px[2];
            double dparallel = sqrt(dpx_squared - dnormal*dnormal);
            double offwalldistance = fmax(dparallel - celllength/2, 0.0);
            double dwall = sqrt(dnormal*dnormal + offwalldistance*offwalldistance);
            cell.dwall = dwall;
            blk.celldata.wall_distances[i] = dwall;
        }
    }
} // end compute_wall_distances()

bool check_full_face_copy(size_t blk_id, size_t blk_bndy, size_t oblk_id, size_t oblk_bndy, JSONValue config_jsonData)
{
/*
    Previously, we checked the actual other block to make sure our full face copy had the
    correct parameters. With the move to initialisating BCs on local blocks only, we need
    to check the config json instead, which is effectively the same thing.

    @author: Nick Gibbons
*/

    string oblk_face_string = face_name[oblk_bndy];

    auto oblk_json      = config_jsonData["block_" ~ to!string(oblk_id)];
    auto oblk_bndy_json = oblk_json["boundary_" ~ oblk_face_string];
    auto oblk_bces      = oblk_bndy_json["pre_recon_action"].array;

    if (oblk_bces.length == 0) return false;
    auto oblk_bce = oblk_bces[0];

    string oblk_bce_type = getJSONstring(oblk_bce, "type", "");
    if (oblk_bce_type != "full_face_copy") return false;

    int oblk_oblk = getJSONint(oblk_bce, "other_block", -1);
    if (oblk_oblk != blk_id) return false;

    string oblk_oface = getJSONstring(oblk_bce, "other_face", "");
    if (face_index(oblk_oface) != blk_bndy) return false;

    return true;
}
