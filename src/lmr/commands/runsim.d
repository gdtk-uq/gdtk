/**
 * Module for launching a steady-state simulation.
 *
 * Authors: RJG, KAD, NNG, PJ
 * Date: 2022-08-13
 * History:
 *   2024-02-12 -- renamed module runsteady --> runsim
 *                 this command will use "solver_mode"
 *                 to decide how the execution is delegated
 */

module runsim;

import core.runtime;
import core.stdc.stdlib : system;
import std.getopt;
import std.stdio : File, write, writeln, writefln;
import std.string;
import std.file : exists;
import core.stdc.stdlib : exit;
import std.parallelism : totalCPUs;
import std.math : FloatingPointControl;
import std.conv : to;
import dyaml;

import json_helper : readJSONfile;

import lmrconfig;
import globalconfig;
import command;
import newtonkrylovsolver : initNewtonKrylovSimulation, performNewtonKrylovUpdates;
import timemarching: initTimeMarchingSimulation, integrateInTime, finalizeSimulation_timemarching;

version(mpi_parallel) {
    import mpi;
}

enum NumberType {default_type, real_values, complex_values};

int determineNumberOfSnapshots()
{
    if (!exists(lmrCfg.restartFile)) return 1;

    auto f = File(lmrCfg.restartFile, "r");
    auto line = f.readln().strip();
    int count = 0;
    while (line.length > 0) {
        if (line[0] != '#') ++count;
        line = f.readln().strip();
    }
    f.close();
    return count + 1; // +1 because restarts file does not have a 0 entry
}

int determineNumberOfTimesEntries()
{
    if (!exists(lmrCfg.timesFile)) return 1;

    Node timesData = dyaml.Loader.fromFile(lmrCfg.timesFile).load();
    return to!int(timesData.length);
}

Command runCmd;
string cmdName = "run";

static this()
{
    // no main field, treated specially
    runCmd.description = "Run a simulation with Eilmer.";
    runCmd.shortDescription = runCmd.description;
    runCmd.helpMsg = format(
`lmr %s [options]

Run an Eilmer simulation.

When invoking this command, the shared memory model of execution is used.
This command assumes that a simulation has been pre-processed
and is available in the working directory.

For distributed memory (using MPI), use the stand-alone executable 'lmr-mpi-run'.
For example:

   $ mpirun -np 4 lmr-mpi-run

options ([+] can be repeated):

 -s, --snapshot-start
     Index of snapshot to use when starting iterations.
     examples:
       -s 1 : start from snapshot 1
       --snapshot-start=3 : start from snapshot 3
       -s -1 : start from final snapshot
       -s=0 : (special case) start from initial condition
     default: none

     NOTE: if the requested snapshot index is greater than
           number of snapshots available, then the iterations will
           begin from the final snapshot.

 --start-with-cfl <or>
 --cfl
     Override the starting CFL on the command line.

     --start-with-cfl=100 : start stepping with cfl of 100
     --cfl 3.5 : start stepping with cfl of 3.5
     default: no override

     NOTE: When not set, the starting CFL comes from input file
           or is computed for the case of a restart.

 --max-cpus=<int>
     Sets maximum number of CPUs for shared-memory parallelism.
     default: %d (on this machine)

     NOTE: in solver_mode=steady, this option has no effect.

 --threads-per-mpi-task=<int>
     Sets threads for MPI tasks when running in MPI mode.
     Leave the default value at 1 unless you know what you're doing
     and know about the distributed/shared memory parallel processing
     model used in Eilmer.
     default: 1

     NOTE: in solver_mode=steady, this option has no effect.

 --max-wall-clock=hh:mm:ss
     This the maximum simulation duration given in hours, minutes and seconds.
     default: 24:00:00
 <OR>
 --max-wall-clock=s
     A single integer value in seconds can also be given.
     A value of -1 disables the max-wall-clock stopping criterion.


 -v, --verbose [+]
     Increase verbosity during progression of the simulation.

`, cmdName, totalCPUs);
}

int delegateAndExecute(string[] args, NumberType numberType)
{
    // We just want to pull out solver mode at this point, so that we can direct the execution flow.
    // We choose to execute these few lines rather than invoking the overhead of reading the entire
    // config file into GlobalConfig.
    auto cfgJSON = readJSONfile(lmrCfg.cfgFile);
    string solverModeStr = cfgJSON["solver_mode"].str;
    auto solverMode = solverModeFromName(solverModeStr);

    string shellStr;

    final switch(solverMode) {
    case SolverMode.steady:
        // Our first choice is to run with complex values.
        final switch(numberType) {
        case NumberType.default_type:
        case NumberType.complex_values:
            shellStr = args[0] ~ "Z-" ~ args[1];
            break;
        case NumberType.real_values:
            shellStr = args[0] ~ "-" ~ args[1];
             break;
	    }
        break;
    case SolverMode.transient:
    case SolverMode.block_marching:
        // Our first choice is to run with real values.
        final switch(numberType) {
        case NumberType.default_type:
        case NumberType.real_values:
            shellStr = args[0] ~ "-" ~ args[1];
            break;
        case NumberType.complex_values:
            shellStr = args[0] ~ "Z-" ~ args[1];
            break;
        }
        break;
    }

    foreach (s; args[2 .. $]) {
        shellStr ~= " " ~ s;
    }

    return system(shellStr.toStringz);
}


/* Why the version(run_main)?
 * Eilmer has a small set of exectuables and each requires a "main" function.
 * The main "main" is the Eilmer command dispatcher and it lives in main.d.
 * Here we need the capability to build two flavours of "run" executable:
 * one for shared memory and one for distributed memory (MPI).
 * We require a "main" function exposed at compile time, but don't want
 * this to collide with the one in main.d.
 * So we wrap this in a version clause and only enable when building those
 * specific executables.
 *
 * RJG, 2024-02-12
 */

version(run_main)
{

int main(string[] args)
{

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
        writeln("Eilmer simulation code.");
        writeln("Revision-id: ", lmrCfg.revisionId);
        writeln("Revision-date: ", lmrCfg.revisionDate);
        writeln("Compiler-name: ", lmrCfg.compilerName);
        writeln("Parallel-flavour: PUT_PARALLEL_FLAVOUR_HERE");
        writeln("Number-type: PUT_NUMBER_TYPE_HERE");
        write("Build-flavour: ");
        version(flavour_debug) { writeln("debug"); }
        version(flavour_profile) { writeln("profile"); }
        version(flavour_fast) { writeln("fast"); }
        writeln("Build-date: ", lmrCfg.buildDate);
    }

    int verbosity = 0;
    int snapshotStart = 0;
    int numberSnapshots = 0;
    int maxCPUs = totalCPUs;
    int threadsPerMPITask = 1;
    double startCFL = -1.0;
    string maxWallClock = "24:00:00";

    getopt(args,
           config.bundling,
           "v|verbose+", &verbosity,
           "s|snapshot-start", &snapshotStart,
           "start-with-cfl|cfl", &startCFL,
           "max-cpus", &maxCPUs,
           "threads-per-mpi-task", &threadsPerMPITask,
           "max-wall-clock", &maxWallClock);

    GlobalConfig.verbosity_level = verbosity;

    if (verbosity > 0 && GlobalConfig.is_master_task) {
        writeln("lmr run: Begin simulation.");
        version(mpi_parallel) {
            writefln("lmr-mpi-run: number of MPI ranks= %d", size);
        }
    }

    // Enable hardware exceptions for division by zero, overflow to infinity,
    // invalid operations, and uninitialized floating-point variables.
    // Copied from https://dlang.org/library/std/math/hardware/floating_point_control.html
    // fpctrl.enableExceptions(FloatingPointControl.severeExceptions);
    //
    // Unfortunately invalidException objects to even writing a NaN
    // to the screen. Since this behaviour is undesirable we really
    // can't leave this one enabled, though the other two are fine.
    // (NNG, Oct 23)
    //
    // When compiling in fast mode (e.g. with -O1, -O2, etc.) the LDC2 compiler
    // optimizes the code to pre-calculate both branches of conditional statements.
    // This can cause problems when floating point exceptions are enabled, since
    // one of the branches, the invalid or 'forbidden' branch, may peform an invalid
    // floating point operation, such as division by zero. We have observed that
    // the routine in the complex number module which handles complex number division
    // exhibits this behaviour. This does not occur for the debug version of the code
    // since no optimizations are enabled, and thus only the 'safe' branch of the code
    // is executed. For this reason we only switch on floating point exceptions for
    // the debug flavour of the complex number version.
    // (KAD, Feb 24)
    version(complex_numbers) {
        debug {
            version(enable_fp_exceptions) {
                FloatingPointControl fpctrl;
                //fpctrl.enableExceptions(FloatingPointControl.invalidException);
                fpctrl.enableExceptions(FloatingPointControl.divByZeroException);
                fpctrl.enableExceptions(FloatingPointControl.overflowException);
            }
        }
    }
    else {
        version(enable_fp_exceptions) {
            FloatingPointControl fpctrl;
            //fpctrl.enableExceptions(FloatingPointControl.invalidException);
            fpctrl.enableExceptions(FloatingPointControl.divByZeroException);
            fpctrl.enableExceptions(FloatingPointControl.overflowException);
        }
    }

    // We just want to pull out solver mode at this point, so that we can direct the execution flow.
    // We choose to execute these few lines rather than invoking the overhead of reading the entire
    // config file into GlobalConfig.
    auto cfgJSON = readJSONfile(lmrCfg.cfgFile);
    string solverModeStr = cfgJSON["solver_mode"].str;
    auto solverMode = solverModeFromName(solverModeStr);

    // Figure out which snapshot to start from
    if (GlobalConfig.is_master_task) {
        if (solverMode == SolverMode.steady) {
            numberSnapshots = determineNumberOfSnapshots();
        }
        else {
            numberSnapshots = determineNumberOfTimesEntries();
        }
        if (snapshotStart == -1) {
            snapshotStart = numberSnapshots-1;
            if (verbosity > 1) {
                writeln("lmr run: snapshot requested is '-1' -- final snapshot");
                writefln("lmr run: starting from final snapshot, index= %02d", snapshotStart);
            }
        }
        if (snapshotStart >= numberSnapshots) {
            if (verbosity > 1) {
                writefln("lmr run: snapshot requested is %02d; this is greater than number of available snapshots", snapshotStart);
                writefln("lmr run: starting from final snapshot, index= %02d", numberSnapshots);
            }
            snapshotStart = numberSnapshots-1;
        }
    }

    version(mpi_parallel) {
        MPI_Bcast(&snapshotStart, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    final switch (solverMode) {
    case SolverMode.steady:
        if (verbosity > 0 && GlobalConfig.is_master_task) writeln("lmr run: Initialise Newton-Krylov simulation.");
        initNewtonKrylovSimulation(snapshotStart, maxCPUs, threadsPerMPITask, maxWallClock);

        if (verbosity > 0 && GlobalConfig.is_master_task) writeln("lmr run: Perform Newton steps.");
        performNewtonKrylovUpdates(snapshotStart, startCFL, maxCPUs, threadsPerMPITask);
        break;
    case SolverMode.transient:
        if (verbosity > 0 && GlobalConfig.is_master_task) writeln("lmr run: Initialise transient simulation.");
        initTimeMarchingSimulation(snapshotStart, maxCPUs, threadsPerMPITask, maxWallClock);

        if (verbosity > 0 && GlobalConfig.is_master_task) writeln("lmr run: Perform integration in time.");
        auto flag = integrateInTime(GlobalConfig.max_time);
        if (flag != 0 && GlobalConfig.is_master_task) writeln("Note that integrateInTime failed.");
        finalizeSimulation_timemarching();
        break;
    case SolverMode.block_marching:
        writeln("NOT IMPLEMENTED: block marching is not available right now.");
	    exit(1);
    }

    GlobalConfig.finalize();

    return 0;
}

} // end: version(run_main)
