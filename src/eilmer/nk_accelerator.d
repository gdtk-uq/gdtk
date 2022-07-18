/** steadystate_solver.d
 * Newton-Krylov updates for steady-state convergence.
 *
 * Author: Rowan G.
 * Date: 2016-10-09
 *
 * Note: This is the first attempt at 'production' code.
 * Some test implementations began on 2016-03-29.
 *
 * History:
 *   2016-10-17 : Implement right-preconditioned GMRES and
 *                add a preconditioner.
 *   2016-10-28 : Add right-preconditioned flexible GMRES so that
 *                we can use a variable preconditioning step.
 *   2016-11-02 : Add restarted GMRES method.
 *   2018-06-07 : Separated solver into two files,
 *                steadystate_solver.d & steadystate_core.d (Author: K. Damm).
 *   2020-04-11 : Renamed nk_accelerator.d and added MPI capability
 */

import core.thread;
import core.runtime;
import std.stdio;
import std.file;
import std.format;
import std.conv;
import std.parallelism;
import std.algorithm;
import std.getopt;
import std.string;
import std.array;
import std.math;
import std.datetime;

import nm.smla;
import nm.bbla;
import nm.complex;
import nm.number;

import steadystate_core;
import special_block_init;
import fluidblock;
import fluidblockio_old;
import sfluidblock;
import globaldata;
import globalconfig;
import simcore : init_simulation;
import fileutil;
import user_defined_source_terms;
import conservedquantities;
import postprocess : readTimesFile;
import loads;
//import shape_sensitivity_core : sss_preconditioner_initialisation, sss_preconditioner;
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
        writeln("Eilmer 4.0 compressible-flow simulation code -- using Newton-Krylov accelerator.");
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
        msg ~= "Usage: e4-nk-dist [OPTIONS]\n";
    } else {
        msg ~= "Usage: e4-nk-shared [OPTIONS]\n";
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

    // Additional memory allocation specific to steady-state solver
    allocate_global_fluid_workspace();
    foreach (blk; localFluidBlocks) {
        blk.allocate_GMRES_workspace();
    }
    allocate_global_solid_workspace();
    foreach (sblk; localSolidBlocks) {
        sblk.allocate_GMRES_workspace();
    }

    /* Check that items are implemented. */
    bool goodToProceed = true;
    /*
    if (GlobalConfig.gmodel_master.n_species > 1) {
        if (GlobalConfig.is_master_task) {
            writeln("Newton-Krylov accelerator not implemented for multiple-species calculations.");
            stdout.flush();
        }
        goodToProceed = false;
    }
    */
    if (!goodToProceed) {
        if (GlobalConfig.is_master_task) {
            writeln("One or more options are not yet available for the Newton-Krylov accelerator.");
            writeln("Bailing out!");
            stdout.flush();
        }
        exitFlag = 1;
        return exitFlag;
    }

    iterate_to_steady_state(snapshotStart, maxCPUs, threadsPerMPITask);

    /* Write residuals to file before exit. */
    if (GlobalConfig.is_master_task) {
        ensure_directory_is_present("residuals");
    }
    version(mpi_parallel) {
        MPI_Barrier(MPI_COMM_WORLD);
    }
    foreach (blk; localFluidBlocks) {
        auto fileName = format!"residuals/%s.residuals.b%04d.gz"(jobName, blk.id);
        blk.write_residuals(fileName);
    }

    if (GlobalConfig.is_master_task) {
        writeln("Done simulation.");
        stdout.flush();
    }
    exitFlag = 0;
    return exitFlag;
}
