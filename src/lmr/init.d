/**
 * Functions to initialise simulation to various levels of readiness.
 *
 * Authors: PAJ and RJG
 * Date: 2023-05-07
 *
 * History:
 *   2023-05-07 Functionality has been pulled together from
 *              eilmer/init_simulation.d
 **/

module lmr.init;

import std.algorithm : min, sort, find;
import std.conv : to;
import std.file : rename, readText, dirEntries, SpanMode, write;
import std.format : format, formattedWrite;
import std.json;
import std.math: pow;
import std.parallelism : parallel, defaultPoolThreads;
import std.stdio : File, writeln, writefln, stdout;
import std.string;
import std.typecons : Tuple, tuple;

import util.json_helper;
import util.lua;
import util.lua_service;

import geom.elements.nomenclature;
import lmr.bc.ghost_cell_effect.gas_solid_full_face_copy;
import lmr.bc.boundary_vertex_full_face_copy;
import lmr.bc;
import lmr.blockio : BinaryBlockIO, GzipBlockIO;
import lmr.fileutil : ensure_directory_is_present;
import lmr.fluidblock : FluidBlock;
import lmr.fvcell : FVCell;
import lmr.fvcellio;
import lmr.globalconfig;
import lmr.globaldata;
import lmr.lmrconfig;
import lmr.lmrexceptions : LmrException;
import lmr.lmrerrors;
import lmr.loads : init_loads_metadata_file, initRunTimeLoads;
import lmr.lua_helper;
import lmr.sfluidblock : SFluidBlock;
import lmr.simcore;
import lmr.simcore_exchange : exchange_ghost_cell_geometry_data;
import lmr.solid.solid_full_face_copy : SolidGCE_SolidGhostCellFullFaceCopy;
import lmr.solid.solid_gas_full_face_copy;
import lmr.solid.solidfvinterface : initPropertiesAtSolidInterfaces;
import lmr.solid.ssolidblock : SSolidBlock;
import lmr.ufluidblock : UFluidBlock;
import lmr.efield.efield: ElectricField;

version(mpi_parallel) {
    import mpi;
}

/**
 * Read from config file and set some global parameters.
 *
 * At the end of this function, configuration is set for:
 *  + global parameters; and
 *  + configuration for each block
 *
 * Authors: RJG and PAJ
 * Date: 2023-05-7
 */
JSONValue initConfiguration()
{
    // Read in config file and set parameters
    auto cfgData = readJSONfile(lmrCfg.cfgFile);
    set_config_for_core(cfgData);
    set_config_for_blocks(cfgData);
    return cfgData;
}

/**
 * Read from the control file.
 *
 * We will need to do this at initialisation, and then repeatedly
 * throughout the simulation to control certain behaviours
 * of the time-marching mode. How often we re-read this is
 * set by GlobalConfig.control_count.
 *
 * Authors: RJG and PAJ
 * Date: 2024-02-07
 */
void readControl()
{
    alias cfg = GlobalConfig;
    if (cfg.verbosity_level > 2) writeln("readControl()");
    JSONValue jsonData = readJSONfile(lmrCfg.ctrlFile);
    mixin(update_double("dt_init", "dt_init"));
    mixin(update_double("dt_max", "dt_max"));
    mixin(update_double("cfl_scale_factor", "cfl_scale_factor"));
    mixin(update_bool("stringent_cfl", "stringent_cfl"));
    mixin(update_double("viscous_signal_factor", "viscous_signal_factor"));
    mixin(update_double("turbulent_signal_factor", "turbulent_signal_factor"));
    mixin(update_enum("residual_smoothing_type", "residual_smoothing_type", "residual_smoothing_type_from_name"));
    mixin(update_double("residual_smoothing_weight", "residual_smoothing_weight"));
    mixin(update_int("residual_smoothing_iterations", "residual_smoothing_iterations"));
    mixin(update_bool("fixed_time_step", "fixed_time_step"));
    mixin(update_int("print_count", "print_count"));
    mixin(update_int("cfl_count", "cfl_count"));
    mixin(update_double("max_time", "max_time"));
    mixin(update_int("max_step", "max_step"));
    // mixin(update_double("dt_plot", "dt_plot")); // 2024-08-07 moved to dt_plot_schedule in config file.
    mixin(update_double("dt_history", "dt_history"));
    mixin(update_double("dt_loads", "dt_loads"));
    mixin(update_int("write_loads_at_step", "write_loads_at_step"));
    mixin(update_int("write_flow_solution_at_step", "write_flow_solution_at_step"));
    mixin(update_int("snapshot_count", "snapshotCount"));
    mixin(update_int("number_total_snapshots", "nTotalSnapshots"));
    //
    mixin(update_int("halt_now", "halt_now"));
    //
    if (cfg.verbosity_level > 2) {
        writeln("  dt_init: ", cfg.dt_init);
        writeln("  dt_max: ", cfg.dt_max);
        writeln("  cfl_scale_factor: ", cfg.cfl_scale_factor);
        writeln("  stringent_cfl: ", cfg.stringent_cfl);
        writeln("  viscous_signal_factor: ", cfg.viscous_signal_factor);
        writeln("  turbulent_signal_factor: ", cfg.turbulent_signal_factor);
        writeln("  residual_smoothing_type: ", cfg.residual_smoothing_type);
        writeln("  residual_smoothing_weight: ", cfg.residual_smoothing_weight);
        writeln("  residual_smoothing_iterations: ", cfg.residual_smoothing_iterations);
        writeln("  fixed_time_step: ", cfg.fixed_time_step);
        writeln("  print_count: ", cfg.print_count);
        writeln("  cfl_count: ", cfg.cfl_count);
        writeln("  max_time: ", cfg.max_time);
        writeln("  max_step: ", cfg.max_step);
        // writeln("  dt_plot: ", cfg.dt_plot); // 2024-08-07 moved to dt_plot_schedule in config file.
        writeln("  dt_history: ", cfg.dt_history);
        writeln("  dt_loads: ", cfg.dt_loads);
        writeln("  write_loads_at_step: ", cfg.write_loads_at_step);
        writeln("  write_flow_solution_at_step: ", cfg.write_flow_solution_at_step);
        writeln("  snapshot_count: ", cfg.snapshotCount);
        writeln("  number_total_snapshots: ", cfg.nTotalSnapshots);
        writeln("  halt_now: ", cfg.halt_now);
    }
    // Propagate new values to the local copies of config.
    foreach (localConfig; dedicatedConfig) {
        localConfig.update_control_parameters();
    }
}

/**
 * Write 0th entry into progress file for time-marching simulations
 *
 * Creating the 0th entry *is* the initialisation.
 *
 * Authors: RJG and PJ
 * Date: 2025-04-02
 */
void initTimeMarchingProgressFile() {
    if (GlobalConfig.is_master_task) {
        try {
            write(lmrCfg.progFile, format("%d\n", SimState.step));
        }
        catch (Exception e) {
            lmrErrorExit(LmrError.inputOutput, "Couldn't write to file: " ~ lmrCfg.progFile);
        }
    }
}

/**
 * Assign (fluid/solid) blocks to the local(Fluid/Solid)Blocks container and set IDs.
 *
 * Each process only works on updating blocks that reside in local(Fluid)Blocks.
 * For MPI simulations, local(Fluid)Blocks contains different blocks on different ranks.
 * This function is responsible for doing that assignment of blocks at the
 * initialisation stage.
 *
 * For a shared-memory simulation, *all* blocks are assigned to the localFluidBlocks container.
 *
 * Additionally the block ids are set in local(Fluid)BlockIds
 *
 * Authors: RJG and PAJ
 * Date: 2023-05-07
 * History:
 *   2024-02-25 -- add in handling for solid blocks
 *
 */
void initLocalBlocks()
{
    alias cfg = GlobalConfig;
    version(mpi_parallel) {
        // Assign particular fluid (& solid) blocks to this MPI task and keep a record
        // of the MPI rank for all blocks.
        int my_rank = cfg.mpi_rank_for_local_task;
        cfg.mpi_rank_for_block.length = cfg.nFluidBlocks + cfg.nSolidBlocks ;
        auto lines = readText(lmrCfg.mpimapFile).splitLines();
        foreach (line; lines) {
            auto content = line.strip();
            if (content.startsWith("#")) continue; // Skip comment
            auto tokens = content.split();
            int blkid = to!int(tokens[0]);
            int taskid = to!int(tokens[1]);
            if (taskid >= cfg.mpi_size && cfg.is_master_task) {
                writefln("Number of MPI tasks (%d) is insufficient for "~
                         "taskid=%d that is associated with blockid=%d. Quitting.",
                         cfg.mpi_size, taskid, blkid);
                MPI_Abort(MPI_COMM_WORLD, 2);
            }
            cfg.mpi_rank_for_block[blkid] = taskid;
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
    }
    else {
        foreach (blk; globalBlocks) {
            auto fblk = cast(FluidBlock) blk;
            if (fblk) { localFluidBlocks ~= fblk; }
            auto mysblk = cast(SSolidBlock) blk;
            if (mysblk) { localSolidBlocks ~= mysblk; }
        }
    }
    // Set block IDs
    foreach (blk; localFluidBlocks) { cfg.localFluidBlockIds ~= blk.id; }
    foreach (blk; localSolidBlocks) { cfg.localSolidBlockIds ~= blk.id; }
} // end initLocalBlocks()

/**
 * Set the pool of threads on a per process basis.
 *
 * This function computes number of threads for each process and
 * sets up the thread pool with a call to defaultPoolThreads().
 *
 * Authors: RJG and PAJ
 * Date: 2023-05-07
 */
void initThreadPool(int maxCPUs, int threadsPerMPITask)
{
    // Local blocks may be handled with thread-parallelism.
    auto nBlocksInThreadParallel = localFluidBlocks.length + localSolidBlocks.length;
    // There is no need to have more task threads than blocks local to the process.
    int extraThreadsInPool;
    version(mpi_parallel) {
	extraThreadsInPool = min(threadsPerMPITask-1, nBlocksInThreadParallel-1);
    }
    else {
	extraThreadsInPool = min(maxCPUs-1, nBlocksInThreadParallel-1);
    }
    defaultPoolThreads(extraThreadsInPool); // total = main thread + extra-threads-in-Pool
    version(mpi_parallel) {
	if (GlobalConfig.verbosity_level > 0) {
	    debug {
		int my_rank = GlobalConfig.mpi_rank_for_local_task;
		writeln("MPI-task with rank ", my_rank, " running with ", extraThreadsInPool+1, " threads.");
	    }
	}
    }
    else {
	if (GlobalConfig.verbosity_level > 0) {
	    writeln("Single process running with ", extraThreadsInPool+1, " threads.");
	}
    }
} // end initThreadPool()

/**
 * Initialise basic components of all fluid blocks.
 *
 * This function instructs blocks to:
 *   + initialise gas model
 *   + initialise any workspace in the block
 *   + set globals for Lua state used by user-defined hooks
 *   + complete the construction of boundary conditions
 *     (now that more info is available to the block)
 *
 * Authors: RJG and PAJ
 * Date: 2023-05-07
 */
void initFluidBlocksBasic(JSONValue cfgData, bool withUserPad=false)
{
    foreach (myblk; localFluidBlocks) {
        myblk.myConfig.init_gas_model_bits();
        myblk.init_workspace();
        myblk.init_lua_globals();
        myblk.init_boundary_conditions(cfgData["block_" ~ to!string(myblk.id)]);
        foreach (bci; myblk.bc) { bci.post_bc_construction(); } // TODO: Possible remove this.
        if (withUserPad && GlobalConfig.user_pad_length > 0) {
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
}

/**
 * Perform the bulk of the memory allocation for all blocks.
 *
 * To do memory allocation, we'll need to read grids and do
 * some associated geometry calculations.
 *
 * This function instructs blocks to:
 *   + read grids
 *   + initialise array space for grid and flow storage
 *   + compute and store geometry data for primary cells
 *   + set up the space for IO per block
 *
 * Authors: RJG and PAJ
 * Date: 2023-05-07
 */
void initFluidBlocksMemoryAllocation()
{
    shared bool anyBlockFail = false;
    foreach (blk; parallel(localFluidBlocks,1)) {
        try {
            string gName;
            if (GlobalConfig.grid_motion != GridMotion.none) {
                gName = gridFilename(SimState.current_tindx, blk.id);
            }
            else {
                gName = gridFilename(lmrCfg.initialFieldDir, blk.id);
            }
            if (GlobalConfig.verbosity_level > 1)
                writeln("Calling init_grid_and_flow_arrays for grid: ", gName);
            blk.init_grid_and_flow_arrays(gName);
            blk.compute_primary_cell_geometric_data(0);
        }
        catch (Exception e) {
            writefln("Block[%d] failed to initialise geometry, msg=%s", blk.id, e.msg);
            anyBlockFail = true;
        }
    }
    version(mpi_parallel) {
        int myFlag = to!int(anyBlockFail);
        MPI_Allreduce(MPI_IN_PLACE, &myFlag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        anyBlockFail = to!bool(myFlag);
    }
    if (anyBlockFail) {
        throw new LmrException("Failed at initialisation stage during grid reading and geometry calculations.");
    }
}

/**
 * Set start ID on individual blocks so that global cell IDs can be computed when needed.
 *
 * Authors: RJG and PAJ
 * Date: 2023-05-07
 */
void initFluidBlocksGlobalCellIDStarts()
{
    // Note that the global id is across all processes, not just the local collection of blocks.
    foreach (i, blk; globalBlocks) {
        auto fluidblk = cast(FluidBlock) blk;
        if (fluidblk) {
            if (i == 0) {
                fluidblk.globalCellIdStart = 0;
            } else {
                auto prev_fluidblk = cast(FluidBlock) globalBlocks[i-1];
                fluidblk.globalCellIdStart = prev_fluidblk.globalCellIdStart + prev_fluidblk.ncells_expected;
            }
        }
    }
}


/**
 * Initialise special zones for all blocks.
 *
 * Authors: RJG and PAJ
 * Date: 2023-05-07
 */
void initFluidBlocksZones()
{
    foreach (blk; localFluidBlocks) {
        // 2023-07-04 PJ: Split this loop out and make it serial
        // because we are having trouble with parallel and GC.
        blk.identify_reaction_zones(0);
        blk.identify_turbulent_zones(0);
        blk.identify_suppress_reconstruction_zones();
        blk.identify_suppress_viscous_stresses_zones();
    }
}

/**
 * Initiliase flow fields for all blocks.
 *
 * This function reads from a snapshot to fill in the flow field
 * data as for all blocks to be used as an initial condition.
 *
 * Authors: RJG and PAJ
 * Date: 2023-05-07
 */
void initFluidBlocksFlowField(int snapshotStart)
{
    // Note that we are not doing the work here in parallel, for a shared-memory run.
    // [TODO] PJ 2024-03-05 Maybe we could.
    bool anyBlockFail = false;
    if (GlobalConfig.field_format == "rawbinary")
        fluidBlkIO = new BinaryBlockIO();
    else
        fluidBlkIO = new GzipBlockIO();

    fluidBlkIO.readMetadataFromFile(lmrCfg.fluidMetadataFile);

    // Note that we have already just computed the cell centroids based on the vertex
    // positions, so we won't read them in from the snapshot file upon start-up.
    string[3] varsToSkip = ["pos.x", "pos.y", "pos.z"];

    foreach (blk; localFluidBlocks) {
        FVCell[] cells;
        cells.length = blk.cells.length;
        foreach (i, ref c; cells) c = blk.cells[i];
        fluidBlkIO.readVariablesFromFile(fluidFilename(snapshotStart, blk.id), cells, varsToSkip);
        // Note that, even for grid_motion==none simulations, we use the grid velocities for setting
        // the gas velocities at boundary faces.  These will need to be set to zero for correct viscous simulation.
        foreach (iface; blk.faces) iface.gvel.clear();
        foreach (cell; blk.cells) {
            // Set rho_s now. It's not read in, and it's not part of encode_conserved
            version(multi_species_gas) {
                foreach (i; 0 .. cell.fs.gas.massf.length) cell.fs.gas.rho_s[i] = cell.fs.gas.massf[i] * cell.fs.gas.rho;
            }
            cell.encode_conserved(0, 0, blk.omegaz);
            // Even though the following call appears redundant at this point,
            // fills in some gas properties such as Prandtl number that is
            // needed for both the cfl_check and the BaldwinLomax turbulence model.
            if (0 != cell.decode_conserved(0, 0, blk.omegaz)) {
                writefln("Block[%d] Bad cell decode_conserved while initialising flow.", blk.id);
                anyBlockFail = true;
            }
        }
        blk.set_cell_dt_chem(-1.0);
    }
    version(mpi_parallel) {
        int myFlag = to!int(anyBlockFail);
        MPI_Allreduce(MPI_IN_PLACE, &myFlag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        anyBlockFail = to!bool(myFlag);
    }
    if (anyBlockFail) {
        throw new LmrException("Failed at initialisation stage during flow field initialisation.");
    }
    /*
    writeln("blk-0, cell-0");
    writeln(localFluidBlocks[0].cells[0]);
    writeln("blk-1, cell-0");
    writeln(localFluidBlocks[1].cells[0]);
    import core.stdc.stdlib : exit;
    exit(1);
    */
}

/**
 * Initialise the full-face exchanges at block connections.
 * Notes:
 *  - Edited by NNG to use config file instead of the globalBlocks
 *    boundary conditions + named change (2025-02-11)
 *
 *
 * Authors: RJG and PAJ
 * Date: 2023-05-07
 */
void initFullFaceDataExchange(JSONValue cfgData)
{
    bool anyBlockFail = false;
    foreach (blk; localFluidBlocks) {
        foreach (j, bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto my_gce = cast(GhostCellFullFaceCopy)gce;
                if (my_gce) {
                    // The local block thinks that it has an exchange boundary with another block,
                    // so we need to check the ghost-cell effects of the other block's face to see
                    // that it points back to the local block face.

                    size_t other_blk_id   = my_gce.neighbourBlock.id;
                    size_t other_blk_bndy = my_gce.neighbourFace;
                    bool ok = checkFullFaceCopyOkay(blk.id, j, other_blk_id, other_blk_bndy, cfgData);
                    if (!ok) {
                        string msg = format("FullFaceCopy for local blk_id=%d face=%d", blk.id, j);
                        msg ~= format(" is not correctly paired with other block id=%d face=%d.",
                                      other_blk_id, my_gce.neighbourFace);
                        writeln(msg);
                        anyBlockFail = true;
                    }
                }
            }
        }
    }
    version(mpi_parallel) {
        int myFlag = to!int(anyBlockFail);
        MPI_Allreduce(MPI_IN_PLACE, &myFlag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        anyBlockFail = to!bool(myFlag);
    }
    if (anyBlockFail) {
        throw new LmrException("Failed at initialisation stage during full-face boundary data exchange.");
    }
}

bool checkFullFaceCopyOkay(size_t blk_id, size_t blk_bndy, size_t oblk_id, size_t oblk_bndy, JSONValue config_jsonData)
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


/**
 * Initialise the mapped-cell exchanges at block connections.
 *
 * Authors: RJG and PAJ
 * Date: 2023-05-07
 */
void initMappedCellDataExchange()
{
    // Serial loops follow because the cell-mapping function searches across
    // all blocks local to the process.
    // Also, there are several loops because the MPI communication,
    // if there is any, needs to be done in phases of posting of non-blocking reads,
    // followed by all of the sends and then waiting for all requests to be filled.
    //
    bool anyBlockFail = false;
    foreach (blk; localFluidBlocks) {
        foreach (bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellMappedCellCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping(); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach (bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellFullFaceCopy)gce;
                if (mygce && (mygce.check_cell_mapping() != 0)) { anyBlockFail = true; }
            }
        }
    }
    version(mpi_parallel) {
        int myFlag = to!int(anyBlockFail);
        MPI_Allreduce(MPI_IN_PLACE, &myFlag, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        anyBlockFail = to!bool(myFlag);
    }
    if (anyBlockFail) {
        throw new LmrException("Failed at initialisation stage during locating mapped-cell boundaries.");
    }

    foreach (blk; localFluidBlocks) {
        foreach (bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellFullFaceCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping_phase0(); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach (bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellFullFaceCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping_phase1(); }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach (bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellFullFaceCopy)gce;
                if (mygce) { mygce.set_up_cell_mapping_phase2(); }
            }
        }
    }
}

void initVertexPositionExchange()
{
    // Instantiate the exchangers. The vertex exchanger is instantiated for each
    // GhostCellFullFaceCopy, and gets all the information it needs from the
    // GhostCelLFullFaceCopy.
    foreach (blk; localFluidBlocks) {
        foreach (bc; blk.bc) {
            foreach (gce; bc.preReconAction) {
                auto mygce = cast(GhostCellFullFaceCopy) gce;
                if (mygce) {
                    bc.vertex_exchange = new BoundaryVertexFullFaceCopy(mygce.blk.id,
                                                                        mygce.which_boundary, 
                                                                        mygce.neighbourBlock.id,
                                                                        mygce.neighbourFace);
                    bc.vertex_exchange.setup_vertex_mapping_phase0();
                }
            }
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach (bc; blk.bc) {
            if (bc.vertex_exchange) bc.vertex_exchange.setup_vertex_mapping_phase1();
        }
    }
    foreach (blk; localFluidBlocks) {
        foreach (bc; blk.bc) {
            if (bc.vertex_exchange) bc.vertex_exchange.setup_vertex_mapping_phase2();
        }
    }
}

void initShockFitting(JSONValue cfgData)
{
    // For a simulation with shock fitting, the files defining the rails for
    // vertex motion and the scaling of vertex velocities throughout the blocks
    // will have been written by prep.lua + output.lua.
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
            fba.read_shockfitting_inflow(cfgData, i);
        }
    }

    initVertexPositionExchange();
}

/**
 * Initialise the ghost cell geometry.
 *
 * This function can be called once a block knows
 * everything about its own cell geometry AND
 * who it is connected to. At that point, that
 * information can be exchanged with connecting neighbours.
 *
 * Authors: RJG and PAJ
 * Date: 2023-05-07
 */
void initGhostCellGeometry()
{
    exchange_ghost_cell_geometry_data();
}

/**
 * Initialise the stencils for least-squares reconstruction.
 *
 * Authors: KAD and PAJ
 * Date: 2023-05-07
 */
void initLeastSquaresStencils()
{
    foreach (blk; localFluidBlocks) blk.compute_least_squares_setup(0);
}

/**
 * Compute the structured interpolation data.
 *
 * Authors: NNG
 * Date: 2024-11-06
 */
void initStructuredStencilData()
{
    foreach (blk; localFluidBlocks) blk.precompute_stencil_data(0);
}
/**
 * Initialise the unstructured grid limiters (for unstructured blocks).
 *
 * Authors: KAD
 * Date: 2024-08-06
 */
void initUSGlimiters()
{
    auto sblock = cast(SFluidBlock) localFluidBlocks[0];
    if (sblock) return; // we don't need a reference length for structured grid simulations

    double dim = to!double(GlobalConfig.dimensions);
    double length = 0.0;

    // compute a reference length for the computational domain
    foreach (blk; localFluidBlocks) {
        foreach (cell; blk.cells) { length += cell.volume[0].re; }
    }
    version(mpi_parallel) {
        MPI_Allreduce(MPI_IN_PLACE,&length,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    }
    length = pow(length, 1.0/dim);

    if (GlobalConfig.is_master_task) { writeln("Computational Domain Reference Length (m): ", length); }

    // set the reference length in the gradients objects for later use in computing the unstructured grid limiter
    foreach (blk; localFluidBlocks) {
        foreach (cell; blk.cells) { cell.gradients.lref = length; }
    }
}

/**
 * Initialise the MPL limiter (for unstructured blocks).
 *
 * Authors: KAD and RJG
 * Date: 2023-05-07
 */
void initMLPlimiter()
{
    foreach (blk; localFluidBlocks) {
        auto ublock = cast(UFluidBlock) blk;
        if (ublock) { ublock.build_cloud_of_cell_references_at_each_vertex(); }
    }
}

/**
 * Initialise the master Lua state (for global user-defined actions).
 *
 * Authors: RJG and PAJ
 * Date: 2023-05-07
 */
void initMasterLuaState()
{
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
    lua_pushnumber(L, GlobalConfig.nBlocks);
    lua_setglobal(L, "nBlocks");
    lua_pushnumber(L, GlobalConfig.nFluidBlocks);
    lua_setglobal(L, "nFluidBlocks");
    lua_pushnumber(L, GlobalConfig.nSolidBlocks);
    lua_setglobal(L, "nSolidBlocks");
    lua_pushnumber(L, GlobalConfig.dimensions);
    lua_setglobal(L, "dimensions");
    // [TODO] think about user_pad -- does it have a use case in steady-state?
}

/**
 * Set up global view of all block corners.
 *
 * Authors: RJG and PAJ
 * Date: 2023-05-07
 */
void initCornerCoordinates()
{
    synchronize_corner_coords_for_all_blocks();
}

/**
 * Compute distances to wall.
 *
 * Authors: RJG and PAJ
 * Date: 2023-05-07
 */
void initWallDistances()
{
    compute_wall_distances();
}

/**
 * Initialise run-time loads computation.
 *
 * Authors: RJG
 * Date: 2024-08-12
 */
void initRunTimeLoads()
{
    string content = readText(lmrCfg.cfgFile);
    JSONValue jsonData = parseJSON!string(content);
    initRunTimeLoads(jsonData["run_time_loads"]);
}

/**
 * Initialise the electric field solver
 *
 * Authors: NNG
 * Date: 2025-09-04
 */
void initElectricField()
{
    // eField lives in globaldata.d
    eField = new ElectricField(localFluidBlocks, GlobalConfig.conductivity_model_name);
}


void orderBlocksBySize()
{
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

/**
 * Initialise solid blocks.
 *
 * Authors: RJG and KAD
 * Date: 2024-02-25
 */
void initSolidBlocks()
{
    solidBlkIO = (GlobalConfig.field_format == "rawbinary") ? new BinaryBlockIO() : new GzipBlockIO();
    solidBlkIO.readMetadataFromFile(lmrCfg.solidMetadataFile);
    foreach (ref solidBlk; localSolidBlocks) {
        solidBlk.assembleArrays();
        solidBlk.bindFacesAndVerticesToCells();
        solidBlk.bindCellsToFaces();
        auto gName = gridFilename(SimState.current_tindx, solidBlk.id);
        solidBlk.readGrid(gName);
        FVCell[] cells;
        cells.length = solidBlk.cells.length;
        foreach (i, ref c; cells) c = solidBlk.cells[i];
        solidBlkIO.readVariablesFromFile(solidFilename(SimState.current_tindx, solidBlk.id), cells);
        solidBlk.computePrimaryCellGeometricData();
        solidBlk.assignCellLocationsForDerivCalc();
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
}

/**
 * Initialise the boundary conditions related to fluid-solid interfaces.
 *
 * Authors: KAD and RJG
 * Date: 2024-03-04
 */
void initFluidSolidExchangeBoundaries()
{
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
}

/**
 * Initialise the loads writing infrastructure.
 *
 * Authors: RJG
 * Date: 2024-03-11
 */
void initLoadsFiles()
{
    ensure_directory_is_present(lmrCfg.loadsDir);
    init_loads_metadata_file();
}
