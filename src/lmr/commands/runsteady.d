/**
 * Module for launching a steady-state simulation.
 *
 * Authors: RJG, KAD, NNG, PJ
 * Date: 2022-08-13
 */

module runsteady;

import core.runtime;
import std.getopt;
import std.stdio : writeln;

import globalconfig;
import command;
import newtonkrylovsolver : initNewtonKrylovSimulation, performNewtonKrylovUpdates;

version(mpi_parallel) {
    import mpi;
}

Command runSteadyCmd;

static this()
{
    runSteadyCmd.main = &command.callShellCommand;
    runSteadyCmd.description = "Run a steady-state simulation with Eilmer.";
    runSteadyCmd.shortDescription = runSteadyCmd.description;
    runSteadyCmd.helpMsg =
`lmr run-steady

Run a simulation in steady-state mode.

When invoking this command, the shared memory model of execution is used.
This command assumes that a simulation has been pre-processed
and is available in the working directory.

For distributed memory (using MPI), use the stand-alone executable 'lmr-mpi-run-steady'.
For example:

   $ mpirun -np 4 lmr-mpi-run-steady

`;
}

version(runsteady_main)
{

void main(string[] args)
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
        writeln("Eilmer simulation code: steady-state solver mode.");
        writeln("Revision-id: PUT_REVISION_STRING_HERE");
        writeln("Revision-date: PUT_REVISION_DATE_HERE");
        writeln("Compiler-name: PUT_COMPILER_NAME_HERE");
	writeln("Parallel-flavour: PUT_PARALLEL_FLAVOUR_HERE");
        writeln("Build-date: PUT_BUILD_DATE_HERE");

	writeln("debug.....");
	writeln("args= ", args);
	writeln("size= ", GlobalConfig.mpi_size);
    }
    
    int verbosity = 0;
    int snapshotStart = 0;
    int maxCPUs = 1;
    int threadsPerMPITask = 1;
    string maxWallClock = "24:00:00";

    getopt(args,
           config.bundling,
           "v|verbose+", &verbosity,
           "s|snapshot", &snapshotStart,
           "max-cpus", &maxCPUs,
           "threads-per-mpi-task", &threadsPerMPITask,
           "max-wall-clock", &maxWallClock);

    if (verbosity > 0) {
        writeln("lmr run-steady: Begin Newton-Krylov simulation.");
    }

    if (verbosity > 1) writeln("lmr run-steady: Initialise simulation.");
    initNewtonKrylovSimulation(snapshotStart, maxCPUs, threadsPerMPITask, maxWallClock);

    if (verbosity > 1) writeln("lmr run-steady: Perform Newton steps.");
    performNewtonKrylovUpdates(snapshotStart, maxCPUs, threadsPerMPITask);

    return;
}

}
