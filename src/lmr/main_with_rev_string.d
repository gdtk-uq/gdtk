/**
 * Module for launching a steady-state simulation.
 *
 * Authors: RJG, KAD, NNG, PJ
 * Date: 2022-08-13
 */

module runsteady;

import core.runtime;
import std.getopt;
import std.stdio : File, writeln, writefln;
import std.string;
import std.file : exists;

import lmrconfig;
import globalconfig;
import command;
import newtonkrylovsolver : initNewtonKrylovSimulation, performNewtonKrylovUpdates;

version(mpi_parallel) {
    import mpi;
}

int determineNumberOfSnapshots()
{
    if (!exists(lmrCfg.restartFile))
	return 0;
	
    auto f = File(lmrCfg.restartFile, "r");
    auto line = f.readln().strip();
    int count = 0;
    while (line.length > 0) {
	if (line[0] != '#')
	    ++count;
	line = f.readln().strip();
    }
    f.close();
    return count;
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

 --start-with-cfl|--cfl
     Override the starting CFL on the command line.

     --start-with-cfl=100 : start stepping with cfl of 100
     --cfl 3.5 : start stepping with cfl of 3.5


     NOTE: When not set, the starting CFL comes from input file
           or is computed for the case of a restart.

 -v, --verbose [+]
     Increase verbosity during progression of iterative solver.

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
        writeln("Revision-id: 29bd1e9c");
        writeln("Revision-date: Sun Jul 2 19:11:09 2023 +1000");
        writeln("Compiler-name: ldc2");
	writeln("Parallel-flavour: shared");
        writeln("Build-date: Mon 03 Jul 2023 13:00:19 AEST");
    }
    
    int verbosity = 0;
    int snapshotStart = 0;
    int numberSnapshots = 0;
    int maxCPUs = 1;
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
	writeln("lmr run-steady: Begin Newton-Krylov simulation.");
	version(mpi_parallel) {
	    writefln("lmr-mpi-run-steady: number of MPI ranks= %d", size);
	}
    }

    // Figure out which snapshot to start from
    if (GlobalConfig.is_master_task) {
	numberSnapshots = determineNumberOfSnapshots();
	if (snapshotStart == -1) {
	    snapshotStart = numberSnapshots;
	    if (verbosity > 1) {
		writeln("lmr run-steady: snapshot requested is '-1' -- final snapshot");
		writefln("lmr run-steady: starting from final snapshot, index= %02d", snapshotStart);
	    }
	}
	if (snapshotStart > numberSnapshots) {
	    if (verbosity > 1) {
		writefln("lmr run-steady: snapshot requested is %02d; this is greater than number of available snapshots", snapshotStart);
		writefln("lmr run-steady: starting from final snapshot, index= %02d", numberSnapshots);
	    }
	    snapshotStart = numberSnapshots;
	}
    }
    version(mpi_parallel) {
	MPI_Bcast(&snapshotStart, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    
    if (verbosity > 0 && GlobalConfig.is_master_task) writeln("lmr run-steady: Initialise simulation.");
    initNewtonKrylovSimulation(snapshotStart, maxCPUs, threadsPerMPITask, maxWallClock);

    if (verbosity > 0 && GlobalConfig.is_master_task) writeln("lmr run-steady: Perform Newton steps.");
    performNewtonKrylovUpdates(snapshotStart, startCFL, maxCPUs, threadsPerMPITask);

    return;
}

}
