// simcore_io.d
// 2021-04-15: extracted from simcore.d
//

module simcore_io;

import std.math;
import std.stdio;
import std.file;
import std.conv;
import std.array;
import std.format;
import std.string;
import std.algorithm;
import std.typecons;
import std.parallelism;
import nm.complex;
import nm.number;

import geom;
import geom.misc.kdtree;
import gas;
import fileutil;
import fvcore;
import globalconfig;
import globaldata;
import flowstate;
import fluidblock;
import sfluidblock;
import ufluidblock;
import fluidblockio;
import fluidblockio_old;
import fluidblockio_new;
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
version(mpi_parallel) {
    import mpi;
    import mpi.util;
}

// Keep a record of simulation time and dt for snapshots
Tuple!(double, "t", double, "dt")[] snapshotInfo;


void write_solution_files()
{

    bool legacy = is_legacy_format(GlobalConfig.flow_format);
    FluidBlockIO[] io_list = localFluidBlocks[0].block_io;

    SimState.current_tindx = SimState.current_tindx + 1;
    if (GlobalConfig.verbosity_level > 0 && GlobalConfig.is_master_task) {
        writeln("Write flow solution.");
        stdout.flush();
    }
    if (GlobalConfig.is_master_task) {
        if (legacy) {
            ensure_directory_is_present(make_path_name!"flow"(SimState.current_tindx));
        } else {
            foreach(io; io_list) {
                string path = "CellData/"~io.tag;
                if (io.do_save()) ensure_directory_is_present(make_path_name(path, SimState.current_tindx));
            }
        }
        
        ensure_directory_is_present(make_path_name!"solid"(SimState.current_tindx));
        if (GlobalConfig.grid_motion != GridMotion.none) {
            ensure_directory_is_present(make_path_name!"grid"(SimState.current_tindx));
        }
    }
    version(mpi_parallel) {
        version(mpi_timeouts) {
            MPI_Sync_tasks();
        } else {
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    if (legacy) {
        wait_for_directory_to_be_present(make_path_name!"flow"(SimState.current_tindx));
    } else {
        foreach (io; io_list) {
            string path = "CellData/"~io.tag;
            if (io.do_save()) wait_for_directory_to_be_present(make_path_name(path,SimState.current_tindx));
        }
    }
    auto job_name = GlobalConfig.base_file_name;
    foreach (myblk; parallel(localFluidBlocksBySize,1)) {
        if (legacy) {
            auto file_name = make_file_name!"flow"(job_name, myblk.id, SimState.current_tindx, GlobalConfig.flowFileExt);
            myblk.write_solution(file_name, SimState.time);
        } else {
            foreach(io; myblk.block_io) {
                auto file_name = make_file_name("CellData",io.tag,job_name, myblk.id, SimState.current_tindx, GlobalConfig.flowFileExt);
                if (io.do_save()) io.save_to_file(file_name, SimState.time); 
            }
        }
    }
    wait_for_directory_to_be_present(make_path_name!"solid"(SimState.current_tindx));
    foreach (ref mySolidBlk; localSolidBlocks) {
        auto fileName = make_file_name!"solid"(job_name, mySolidBlk.id, SimState.current_tindx, "gz");
        mySolidBlk.writeSolution(fileName, SimState.time);
    }
    if (GlobalConfig.grid_motion != GridMotion.none) {
        wait_for_directory_to_be_present(make_path_name!"grid"(SimState.current_tindx));
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

void write_snapshot_files()
{
    if (GlobalConfig.nTotalSnapshots == 0) { return; }
    if (snapshotInfo.length == 0) { snapshotInfo.length = GlobalConfig.nTotalSnapshots; }
    if (SimState.nWrittenSnapshots >= GlobalConfig.nTotalSnapshots) {
        // We need to shuffle the existing snapshots down one slot
        auto job_name = GlobalConfig.base_file_name;
        foreach (iSnap; 1 .. GlobalConfig.nTotalSnapshots) {
            foreach (blk; parallel(localFluidBlocksBySize,1)) {
                auto fromName = make_snapshot_file_name("flow", job_name, blk.id, iSnap, GlobalConfig.flowFileExt);
                auto toName = make_snapshot_file_name("flow", job_name, blk.id, iSnap-1, GlobalConfig.flowFileExt);
                rename(fromName, toName);
            }
            foreach (ref solidBlk; localSolidBlocks) {
                auto fromName = make_snapshot_file_name("solid", job_name, solidBlk.id, iSnap, "gz");
                auto toName = make_snapshot_file_name("solid", job_name, solidBlk.id, iSnap-1, "gz");
                rename(fromName, toName);
            }
            if (GlobalConfig.grid_motion != GridMotion.none) {
                foreach (blk; parallel(localFluidBlocksBySize,1)) {
                    auto fromName = make_snapshot_file_name("grid", job_name, blk.id, iSnap, GlobalConfig.gridFileExt);
                    auto toName = make_snapshot_file_name("grid", job_name, blk.id, iSnap-1, GlobalConfig.gridFileExt);
                    rename(fromName, toName);
                }
            }
            snapshotInfo[iSnap-1] = snapshotInfo[iSnap];
        }
    }

    int snapshotIdx = (SimState.nWrittenSnapshots < GlobalConfig.nTotalSnapshots) ? SimState.nWrittenSnapshots : GlobalConfig.nTotalSnapshots-1;
    snapshotInfo[snapshotIdx] = tuple!("t", "dt")(SimState.time, SimState.dt_global);

    if (GlobalConfig.is_master_task) {
        ensure_directory_is_present(make_snapshot_path_name("flow", snapshotIdx));
        ensure_directory_is_present(make_snapshot_path_name("solid", snapshotIdx));
        if (GlobalConfig.grid_motion != GridMotion.none) {
            ensure_directory_is_present(make_snapshot_path_name("grid", snapshotIdx));
        }
    }

    version(mpi_parallel) {
        version(mpi_timeouts) {
            MPI_Sync_tasks();
        } else {
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    wait_for_directory_to_be_present(make_snapshot_path_name("flow", snapshotIdx));
    auto job_name = GlobalConfig.base_file_name;
    foreach (myblk; parallel(localFluidBlocksBySize,1)) {
        auto file_name = make_snapshot_file_name("flow", job_name, myblk.id, snapshotIdx, GlobalConfig.flowFileExt);
        myblk.write_solution(file_name, SimState.time);
    }
    wait_for_directory_to_be_present(make_snapshot_path_name("solid", snapshotIdx));
    foreach (ref mySolidBlk; localSolidBlocks) {
        auto fileName = make_snapshot_file_name("solid", job_name, mySolidBlk.id, snapshotIdx, "gz");
        mySolidBlk.writeSolution(fileName, SimState.time);
    }
    if (GlobalConfig.grid_motion != GridMotion.none) {
        wait_for_directory_to_be_present(make_snapshot_path_name("grid", snapshotIdx));
        if (GlobalConfig.verbosity_level > 0 && GlobalConfig.is_master_task) { writeln("   Write grid"); }
        foreach (blk; parallel(localFluidBlocksBySize,1)) {
            blk.sync_vertices_to_underlying_grid(0);
            auto fileName = make_snapshot_file_name("grid", job_name, blk.id, snapshotIdx, GlobalConfig.gridFileExt);
            blk.write_underlying_grid(fileName);
        }
    }
    // Write out .snapshots file
    // NOTE: We re-write this completely every time a snapshot is saved.
    // This is not a performance hit because snapshots are infrequent,
    // and this file with time records is small.
    if (GlobalConfig.is_master_task) {
        auto fname = format("config/%s.snapshots", GlobalConfig.base_file_name);
        auto f = File(fname, "w");
        f.writeln("# snapshot-indx sim_time dt_global");
        foreach (idx, snap; snapshotInfo) {
            f.writefln("%04d %.18e %.18e", idx, snap.t, snap.dt);
        }
        f.close();
    }

    SimState.nWrittenSnapshots = SimState.nWrittenSnapshots + 1;
} // end write_snapshot_files()

void write_DFT_files()
{
    if (GlobalConfig.is_master_task) {
        ensure_directory_is_present("DFT");
    }

    version(mpi_parallel) {
        version(mpi_timeouts) {
            MPI_Sync_tasks();
        } else {
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    auto job_name = GlobalConfig.base_file_name;

    foreach (myblk; parallel(localFluidBlocksBySize, 1)) {
        auto file_name = make_DFT_file_name(job_name, myblk.id, GlobalConfig.flowFileExt);
        myblk.write_DFT(file_name);
    }
}

