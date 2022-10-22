/**
 * Module for launching a steady-state simulation.
 *
 * Authors: RJG, KAD, NNG, PJ
 * Date: 2022-08-13
 */

module runsteady;

import std.getopt;
import std.stdio : writeln;

import globalconfig;
import command;
import newtonkrylovsolver : initNewtonKrylovSimulation, performNewtonKrylovUpdates;

Command runSteadyCmd;

static this()
{
    runSteadyCmd.main = &command.callShellCommand;
    runSteadyCmd.description = "Run a steady-state simulation with Eilmer.";
    runSteadyCmd.shortDescription = runSteadyCmd.description;
    runSteadyCmd.helpMsg =
`lmr run-steady [options]

`;
}

version(runsteady_main)
{

void main(string[] args)
{
    if (GlobalConfig.is_master_task) {
        writeln("Eilmer simulation code: steady-state solver mode.");
        writeln("Revision-id: PUT_REVISION_STRING_HERE");
        writeln("Revision-date: PUT_REVISION_DATE_HERE");
        writeln("Compiler-name: PUT_COMPILER_NAME_HERE");
        writeln("Build-date: PUT_BUILD_DATE_HERE");
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
