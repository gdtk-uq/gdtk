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

module init;

import std.algorithm : min, sort;
import std.conv : to;
import std.parallelism : parallel, defaultPoolThreads;
import std.file : rename, readText;
import std.stdio : File, writeln, writefln, stdout;
import std.format : format, formattedWrite;
import std.string;
import std.typecons : Tuple, tuple;

import util.lua;
import util.lua_service;
import lua_helper;

import json_helper;
import lmrexceptions : LmrException;
import lmrconfig;
import globalconfig;
import globaldata;
import simcore;
import simcore_exchange : exchange_ghost_cell_geometry_data;
import fluidblockio_new : read_zip_solution;
import bc;
import fluidblock : FluidBlock;
import sfluidblock : SFluidBlock;
import ufluidblock : UFluidBlock;
import blockio : blkIO, BinaryBlockIO, GzipBlockIO;
import fvcellio;

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
void initConfiguration()
{
    // Read in config file and set parameters
    auto cfgData = readJSONfile(lmrCfg.cfgFile);
    set_config_for_core(cfgData);
    set_config_for_blocks(cfgData);
}

/**
 * Assign fluid blocks to the localFluidBlocks container and set IDs.
 *
 * Each process only works on updating blocks that reside in localFluidBlocks.
 * For MPI simulations, localFluidBlocks contains different blocks on different ranks.
 * This function is resposnsible for doing that assignment of blocks at the
 * initialisation stage.
 *
 * For a shared-memory simulation, *all* blocks are assigned to the localFluidBlocks container.
 *
 * Additionally the block ids are set in loclFluidBlockIds
 *
 * Authors: RJG and PAJ
 * Date: 2023-05-07
 */
void initLocalFluidBlocks()
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
		/+ [TODO] Add in solid blocks.
                auto sblk = cast(SSolidBlock) globalBlocks[blkid];
                if (sblk) { localSolidBlocks ~= sblk; }
		+/
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
	    /+ [TODO] add in solid blocks
	     auto mysblk = cast(SSolidBlock) blk;
	     if (mysblk) { localSolidBlocks ~= mysblk; }
	     +/
	}
    }

    // Set block IDs
    foreach (blk; localFluidBlocks) cfg.localFluidBlockIds ~= blk.id;

}

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
    auto nBlocksInThreadParallel = localFluidBlocks.length; // [TODO] add solid blocks
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
}

/**
 * Initialise basic components of all fluid blocks.
 *
 * This function instructs blocks to:
 *   + initialise gas moodel
 *   + initialise any workspace in the block
 *   + set globals for Lua state used by user-defined hooks
 *   + complete the construction of boundary conditions
 *     (now that more info is available to the block)
 *
 * Authors: RJG and PAJ
 * Date: 2023-05-07
 */
void initFluidBlocksBasic()
{
    foreach (myblk; localFluidBlocks) {
        myblk.myConfig.init_gas_model_bits();
        myblk.init_workspace();
        myblk.init_lua_globals();
        foreach (bci; myblk.bc) { bci.post_bc_construction(); }
        // NOTE: Removed userPad in NK solver.
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
    bool anyBlockFail = false;
    foreach (blk; parallel(localFluidBlocks,1)) {
        try {
            string gName = gridFilename(lmrCfg.initialFieldDir, blk.id);
            debug { writeln("Calling init_grid_and_flow_arrays for grid: ", gName); }
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
    foreach (blk; parallel(localFluidBlocks,1)) {
        blk.identify_reaction_zones(0);
        blk.identify_turbulent_zones(0);
        blk.identify_suppress_reconstruction_zones();
        blk.identify_suppress_viscous_stresses_zones();
    }
}

/**
 * Initiliase flow fields for all blocks for steady mode iterations.
 *
 * This function reads from a snapshot to fill in the flow field
 * data as for all blocks to be used as an initial condition.
 *
 * Authors: RJG and PAJ
 * Date: 2023-05-07
 */
void initFluidBlocksFlowFieldSteadyMode(int snapshotStart)
{
    bool anyBlockFail = false;
    if (GlobalConfig.flow_format == "rawbinary")
	blkIO = new BinaryBlockIO();
    else
	blkIO = new GzipBlockIO();

    blkIO.readMetadataFromFile(lmrCfg.flowMetadataFile);

    foreach (blk; parallel(localFluidBlocks,1)) {
        blkIO.readVariablesFromFile(flowFilename(snapshotStart, blk.id), blk.cells);
        foreach (iface; blk.faces) iface.gvel.clear();
        foreach (cell; blk.cells) {
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
}

/**
 * Initialise the full-face exchanges at block connections.
 *
 * Authors: RJG and PAJ
 * Date: 2023-05-07
 */
void initFullFaceDataExchange()
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
                    auto other_blk = my_gce.neighbourBlock;
                    bool ok = false;
                    auto other_blk_bc = other_blk.bc[my_gce.neighbourFace];
                    foreach (gce2; other_blk_bc.preReconAction) {
                        auto other_gce = cast(GhostCellFullFaceCopy)gce2;
                        if (other_gce &&
                            (other_gce.neighbourBlock.id == blk.id) &&
                            (other_gce.neighbourFace == j)) {
                            ok = true;
                        }
                    }
                    if (!ok) {
                        string msg = format("FullFaceCopy for local blk_id=%d face=%d", blk.id, j);
                        msg ~= format(" is not correctly paired with other block id=%d face=%d.",
                                      other_blk.id, my_gce.neighbourFace);
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
