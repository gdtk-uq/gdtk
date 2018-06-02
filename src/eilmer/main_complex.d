/** main_complex.d
 * Eilmer4 compressible-flow simulation code, 
 * stripped-down top-level function with complex flavour.
 *
 * Author: Peter J. Rowan G. anf Kyle D.
 * First code: 2015-02-05
 * Complex flavour: 2018-06-02
 */

import core.memory;
import core.thread;
import core.stdc.stdlib : exit;
import std.stdio;
import std.string;
import std.file;
import std.path;
import std.getopt;
import std.conv;
import std.parallelism;
import std.algorithm;
import std.math;
import nm.complex;
import nm.number;

import geom;
import gas;
import fvcore: FlowSolverException;
import globalconfig;
import simcore;
version(mpi_parallel) {
    import mpi;
    import mpi.util;
}

void moveFileToBackup(string fileName)
{
    if (exists(fileName)) {
        if (exists(fileName~".bak")) { remove(fileName~".bak"); }
        rename(fileName, fileName~".bak");
    }
    return;
}

int main(string[] args)
{
    int exitFlag = 0; // Presume OK in the beginning.
    version(mpi_parallel) {
        // This preamble copied directly from the OpenMPI hello-world example.
        int argc = cast(int)args.length;
        auto argv = args.toArgv();
        MPI_Init(&argc, &argv);
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        GlobalConfig.mpi_rank_for_local_task = rank;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        GlobalConfig.mpi_size = size;
        scope(exit) { MPI_Finalize(); }
        // Make note that we are in the context of an MPI task, presumably, one of many.
        GlobalConfig.in_mpi_context = true;
        GlobalConfig.is_master_task = (GlobalConfig.mpi_rank_for_local_task == 0);
    } else {
        // We are NOT in the context of an MPI task.
        GlobalConfig.in_mpi_context = false;
        GlobalConfig.is_master_task = true;
    }
    //
    string msg = "Usage:                               Comment:\n";
    msg       ~= "e4shared [--job=<string>]            file names built from this string\n";
    msg       ~= "         [--verbosity=<int>]         defaults to 0\n";
    msg       ~= "\n";
    msg       ~= "         [--prep]                    prepare config, grid and flow files\n";
    msg       ~= "\n";
    msg       ~= "         [--run]                     run the simulation over time\n";
    msg       ~= "         [--tindx-start=<int>|last|9999]  defaults to 0\n";
    msg       ~= "         [--next-loads-indx=<int>]   defaults to (final index + 1) of loads.times file\n";
    msg       ~= "         [--max-cpus=<int>]          (e4shared only) defaults to ";
    msg       ~= to!string(totalCPUs) ~" on this machine\n";
    msg       ~= "         [--threads-per-mpi-task=<int>] (e4mpi only) defaults to 1\n";
    msg       ~= "         [--max-wall-clock=<int>]    in seconds\n";
    msg       ~= "         [--report-residuals]        include residual reporting in console output\n";
    msg       ~= "\n";
    msg       ~= "         [--post]                    post-process simulation data\n";
    msg       ~= "         [--list-info]               report some details of this simulation\n";
    msg       ~= "         [--tindx-plot=<int>|all|last|9999]  defaults to last\n";
    msg       ~= "         [--add-vars=\"mach,pitot,total-h,total-p,entropy,molef,conc\"]\n";
    msg       ~= "         [--ref-soln=<filename>]     Lua file for reference solution\n";
    msg       ~= "         [--vtk-xml]                 produce XML VTK-format plot files\n";
    msg       ~= "         [--binary-format]           use binary within the VTK-XML\n";
    msg       ~= "         [--tecplot]                 write a binary szplt file in Tecplot format\n";
    msg       ~= "         [--tecplot-ascii]           write an ASCII (text) file in Tecplot format\n";
    msg       ~= "         [--plot-dir=<string>]       defaults to plot\n";
    msg       ~= "         [--output-file=<string>]    defaults to stdout\n";
    msg       ~= "         [--slice-list=\"blk-range,i-range,j-range,k-range;...\"]\n";
    msg       ~= "         [--surface-list=\"blk,surface-id;...\"]\n";
    msg       ~= "         [--extract-streamline=\"x,y,z;...\"]        streamline locus points\n";
    msg       ~= "         [--track-wave=\"x,y,z(,nx,ny,nz);..\"]      track wave from given point in given plane, default is n=(0,0,1)\n";
    msg       ~= "         [--extract-line=\"x0,y0,z0,x1,y1,z1,n;...\"]    sample along a line in flow blocks\n";
    msg       ~= "         [--extract-solid-line=\"x0,y0,z0,x1,y1,z1,n;...\"]    sample along a line in solid blocks\n";
    msg       ~= "         [--compute-loads-on-group=\"\"]    group tag\n";
    msg       ~= "         [--probe=\"x,y,z;...\"]       locations to sample flow data\n";
    msg       ~= "         [--output-format=<string>]  gnuplot|pretty\n";
    msg       ~= "         [--norms=\"varName,varName,...\"] report L1,L2,Linf norms\n";
    msg       ~= "         [--region=\"x0,y0,z0,x1,y1,z1\"]  limit norms calculation to a box\n";
    msg       ~= "\n";
    msg       ~= "         [--custom-post]             run custom post-processing script\n";
    msg       ~= "         [--script-file=<string>]    defaults to post.lua\n";
    msg       ~= "\n";
    msg       ~= "         [--help]                    writes this message\n";
    if ( args.length < 2 ) {
        if (GlobalConfig.is_master_task) {
            writeln("Too few arguments.");
            write(msg);
            stdout.flush();
        }
        exitFlag = 1;
        return exitFlag;
    }
    string jobName = "";
    int verbosityLevel = 1; // default to having a little information
    bool prepFlag = false;
    bool runFlag = false;
    string tindxStartStr = "0";
    int tindxStart = 0;
    int nextLoadsIndx = -1;
    int maxCPUs = totalCPUs;
    int threadsPerMPITask = 1;
    int maxWallClock = 5*24*3600; // 5 days default
    bool reportResiduals = false;
    bool postFlag = false;
    bool listInfoFlag = false;
    string tindxPlot = "last";
    string addVarsStr = "";
    string luaRefSoln = "";
    bool vtkxmlFlag = false;
    bool binaryFormat = false;
    bool tecplotBinaryFlag = false;
    bool tecplotAsciiFlag = false;
    string plotDir = "plot";
    string outputFileName = "";
    string sliceListStr = "";
    string surfaceListStr = "";
    string extractStreamStr = "";
    string trackWaveStr = "";
    string extractLineStr = "";
    string extractSolidLineStr = "";
    string computeLoadsOnGroupStr = "";
    string probeStr = "";
    string outputFormat = "gnuplot";
    string normsStr = "";
    string regionStr = "";
    bool customPostFlag = false;
    string scriptFile = "post.lua";
    bool helpWanted = false;
    try {
        getopt(args,
               "job", &jobName,
               "verbosity", &verbosityLevel,
               "prep", &prepFlag,
               "run", &runFlag,
               "tindx-start", &tindxStartStr,
               "next-loads-indx", &nextLoadsIndx,
               "max-cpus", &maxCPUs,
               "threads-per-mpi-task", &threadsPerMPITask,
               "max-wall-clock", &maxWallClock,
               "report-residuals", &reportResiduals,
               "post", &postFlag,
               "list-info", &listInfoFlag,
               "tindx-plot", &tindxPlot,
               "add-vars", &addVarsStr,
               "ref-soln", &luaRefSoln,
               "vtk-xml", &vtkxmlFlag,
               "binary-format", &binaryFormat,
               "tecplot", &tecplotBinaryFlag,
               "tecplot-ascii", &tecplotAsciiFlag,
               "plot-dir", &plotDir,
               "output-file", &outputFileName,
               "slice-list", &sliceListStr,
               "surface-list", &surfaceListStr,
               "extract-streamline", &extractStreamStr,
               "track-wave", &trackWaveStr,
               "extract-line", &extractLineStr,
               "extract-solid-line", &extractSolidLineStr,
               "compute-loads-on-group", &computeLoadsOnGroupStr,
               "probe", &probeStr,
               "output-format", &outputFormat,
               "norms", &normsStr,
               "region", &regionStr,
               "custom-post", &customPostFlag,
               "script-file", &scriptFile,
               "help", &helpWanted
               );
    } catch (Exception e) {
        if (GlobalConfig.is_master_task) {
            writeln("Problem parsing command-line options.");
            writeln("Arguments not processed: ");
            args = args[1 .. $]; // Dispose of program name in first argument.
            foreach (myarg; args) writeln("    arg: ", myarg);
            write(msg);
            stdout.flush();
        }
        exitFlag = 1;
        return exitFlag;
    }
    if (verbosityLevel > 0) {
        version(mpi_parallel) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (GlobalConfig.is_master_task) {
                writeln("Eilmer4 compressible-flow simulation code.");
                writeln("Revision: PUT_REVISION_STRING_HERE");
                writefln("MPI-parallel, number of tasks %d", GlobalConfig.mpi_size);
            }
            stdout.flush();
            // Give master_task a chance to be seen first.
            Thread.sleep(dur!("msecs")(100));
            MPI_Barrier(MPI_COMM_WORLD);
            // Now, get all tasks to report.
            writefln("MPI-parallel, start task %d", GlobalConfig.mpi_rank_for_local_task);
            stdout.flush();
            Thread.sleep(dur!("msecs")(100));
            MPI_Barrier(MPI_COMM_WORLD);
        } else {
            writeln("Eilmer4 compressible-flow simulation code.");
            writeln("Revision: PUT_REVISION_STRING_HERE");
            writeln("Shared-memory");
        }
    }
    if (helpWanted) {
        if (GlobalConfig.is_master_task) { write(msg); stdout.flush(); }
        exitFlag = 0;
        return exitFlag;
    }

    if (runFlag) {
        if (jobName.length == 0) {
            writeln("Need to specify a job name.");
            write(msg);
            exitFlag = 1;
            return exitFlag;
        }
        GlobalConfig.base_file_name = jobName;
        GlobalConfig.verbosity_level = verbosityLevel;
        GlobalConfig.report_residuals = reportResiduals;
        // don't ask for more CPUs than available to this process
        maxCPUs = min(max(maxCPUs, 1), totalCPUs);
        threadsPerMPITask = min(threadsPerMPITask, totalCPUs);
        switch (tindxStartStr) {
        case "9999":
        case "last":
            auto times_dict = readTimesFile(jobName);
            auto tindx_list = times_dict.keys;
            sort(tindx_list);
            tindxStart = tindx_list[$-1];
            break;
        default:
            // We assume that the command-line argument was an integer.
            tindxStart = to!int(tindxStartStr);
        } // end switch
        if (verbosityLevel > 0 && GlobalConfig.is_master_task) {
            writeln("Begin simulation with command-line arguments.");
            writeln("  jobName: ", jobName);
            writeln("  tindxStart: ", tindxStart);
            writeln("  maxWallClock: ", maxWallClock);
            writeln("  verbosityLevel: ", verbosityLevel);
            version(mpi_parallel) {
                writeln("  threadsPerMPITask: ", threadsPerMPITask);
            } else {
                writeln("  maxCPUs: ", maxCPUs, " for shared memory-parallelism");
            }
        }
        version(mpi_parallel) {
            stdout.flush();
            Thread.sleep(dur!("msecs")(100));
            MPI_Barrier(MPI_COMM_WORLD);
        }        
        init_simulation(tindxStart, nextLoadsIndx, maxCPUs, threadsPerMPITask, maxWallClock);
        if (verbosityLevel > 0 && GlobalConfig.is_master_task) {
            writeln("starting simulation time= ", simcore.sim_time);
        }
        version(mpi_parallel) {
            stdout.flush();
            Thread.sleep(dur!("msecs")(100));
            MPI_Barrier(MPI_COMM_WORLD);
        }        
        if (GlobalConfig.block_marching) {
            version(mpi_parallel) {
                if (GlobalConfig.is_master_task) {
                    writeln("Do not run a block-marching simulation with MPI parallelism.");
                    stdout.flush();
                }
                exitFlag = 1;
                return exitFlag;
            } else { // NOT mpi_parallel
                march_over_blocks();
                finalize_simulation();
            }
        } else {
            integrate_in_time(GlobalConfig.max_time);
            finalize_simulation();
        }
        if (verbosityLevel > 0 && GlobalConfig.is_master_task) {
            writeln("Done simulation.");
        }
    } // end if runFlag

    //
    return exitFlag;
} // end main()


