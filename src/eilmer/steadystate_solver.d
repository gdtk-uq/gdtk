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
 */

import core.stdc.stdlib : exit;
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
import sfluidblock;
import globaldata;
import globalconfig;
import simcore;
import fvcore;
import fileutil;
import user_defined_source_terms;
import conservedquantities;
import postprocess : readTimesFile;
import loads;
import shape_sensitivity_core : sss_preconditioner_initialisation, sss_preconditioner;

void main(string[] args)
{
    writeln("Eilmer compressible-flow simulation code -- steady state solver.");
    writeln("Revision: PUT_REVISION_STRING_HERE");
    version(mpi_parallel) {
        assert(0, "Steady-state solver is not MPI parallel, yet.");
    }

    string msg = "Usage:                               Comment:\n";
    msg       ~= "e4sss    [--job=<string>]            file names built from this string\n";
    msg       ~= "         [--verbosity=<int>]         defaults to 0\n";
    msg       ~= "\n";
    msg       ~= "         [--snapshot-start=<int>|last] defaults to 0\n";
    msg       ~= "         [--max-cpus=<int>]          defaults to ";
    msg       ~= to!string(totalCPUs) ~" on this machine\n";
    msg       ~= "         [--max-wall-clock=<int>]    in seconds\n";
    msg       ~= "\n";
    msg       ~= "         [--help]                    writes this message\n";
    if ( args.length < 2 ) {
        writeln("Too few arguments.");
        write(msg);
        exit(1);
    }
    string jobName = "";
    int verbosityLevel = 0;
    string snapshotStartStr = "0";
    int snapshotStart = 0;
    int maxCPUs = totalCPUs;
    int maxWallClock = 5*24*3600; // 5 days default
    bool helpWanted = false;
    try {
        getopt(args,
               "job", &jobName,
               "verbosity", &verbosityLevel,
               "snapshot-start", &snapshotStartStr,
               "max-cpus", &maxCPUs,
               "max-wall-clock", &maxWallClock,
               "help", &helpWanted
               );
    } catch (Exception e) {
        writeln("Problem parsing command-line options.");
        writeln("Arguments not processed: ");
        args = args[1 .. $]; // Dispose of program name in first argument
        foreach (arg; args) writeln("   arg: ", arg);
        write(msg);
        exit(1);
    }
    if (helpWanted) {
        write(msg);
        exit(0);
    }

    if (jobName.length == 0) {
        writeln("Need to specify a job name.");
        write(msg);
        exit(1);
    }

    GlobalConfig.base_file_name = jobName;
    GlobalConfig.verbosity_level = verbosityLevel;
    maxCPUs = min(max(maxCPUs, 1), totalCPUs); // don't ask for more than available

    if (verbosityLevel > 0) {
        writeln("Begin simulation with command-line arguments.");
        writeln("  jobName: ", jobName);
        writeln("  snapshotStart: ", snapshotStartStr);
        writeln("  maxWallClock: ", maxWallClock);
        writeln("  verbosityLevel: ", verbosityLevel);
        writeln("  maxCPUs: ", maxCPUs);
    }

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

    writefln("Initialising simulation from snapshot: %d", snapshotStart);
    init_simulation(snapshotStart, -1, maxCPUs, 1, maxWallClock);

    // Additional memory allocation specific to steady-state solver
    allocate_global_workspace();
    foreach (blk; localFluidBlocks) {
        blk.allocate_GMRES_workspace();
    }

    with_k_omega = (GlobalConfig.turbulence_model == TurbulenceModel.k_omega);

    if ( GlobalConfig.extrema_clipping ) {
        writeln("WARNING:");
        writeln("   extrema_clipping is set to true.");
        writeln("   This is not recommended when using the steady-state solver.");
        writeln("   Its use will likely delay or stall convergence.");
        writeln("   Continuing with simulation anyway.");
        writeln("END WARNING.");
    }

    /* Check that items are implemented. */
    bool goodToProceed = true;
    if ( GlobalConfig.gmodel_master.n_species > 1 ) {
        writeln("Steady-state solver not implemented for multiple-species calculations.");
        goodToProceed = false;
    }
    if ( !goodToProceed ) {
        writeln("One or more options are not yet available for the steady-state solver.");
        writeln("Bailing out!");
        exit(1);
    }

    iterate_to_steady_state(snapshotStart, maxCPUs);
    writeln("Done simulation.");
}
