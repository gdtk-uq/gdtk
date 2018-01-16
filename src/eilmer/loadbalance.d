/** loadbalance.d
 * A tool to statically load balance the distribution
 * of blocks to tasks in an MPI simulation.
 *
 * Authors: Rowan G. and Peter J.
 * Date: 2017-05-21
 */

import core.stdc.stdlib : exit;
import std.stdio;
import std.file;
import std.format;
import std.conv;
import std.getopt;
import std.string;
import std.array;
import std.algorithm;
import std.algorithm.searching;
import std.typecons;

import globaldata;
import globalconfig;

void main(string[] args)
{
    writeln("Eilmer compressible-flow simulation code -- load balancing tool.");
    writeln("Revision: PUT_REVISION_STRING_HERE");

    string msg = "Usage:                              Comment:\n";
    msg       ~= "e4loadbalance  [--job=<string>]     name of job\n";
    msg       ~= "               [--ntasks=<int>]     number of tasks allocated to job\n";
    msg       ~= "               [--sweep=\"<int>,<int>\" sweep through a range of ntasks\n";
    msg       ~= "                                    reporting on load balance efficacy\n";
    msg       ~= "               [--help]             write this message\n";
    if ( args.length < 3 ) {
        writeln("Too few arguments.");
        write(msg);
        exit(1);
    }
    string jobName = "";
    int nTasks = 0;
    string sweepStr = "";
    bool helpWanted = false;
    int sweepMin, sweepMax;
    bool performSweep = false;
    try {
        getopt(args,
               "job", &jobName,
               "ntasks", &nTasks,
               "sweep", &sweepStr,
               "help", &helpWanted,
               );
    } catch (Exception e) {
        writeln("Problem parsing command-line options.");
        writeln("Arguments not processed:");
        args = args[1 .. $]; // Dispose of program in first arg
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

    if (nTasks < 1) {
        writeln("Need to specify 'ntasks' greater or equal to 1.");
        write(msg);
        exit(1);
    }

    if (sweepStr.length > 0) {
        auto sweepRange = sweepStr.split(",");
        if (sweepRange.length != 2) {
            writeln("The sweep range should consist of two integers exactly.");
            write(msg);
            exit(1);
        }
        sweepMin = to!int(sweepRange[0]);
        sweepMax = to!int(sweepRange[1]);
        if ( sweepMax <= sweepMin ) {
            writeln("The sweep range should be given in order 'min,max'");
            writeln("The provided values are:");
            writeln("   sweepMin= ", sweepMin);
            writeln("   sweepMax= ", sweepMax);
            write(msg);
            exit(1);
        }
        performSweep = true;
    }
      
    GlobalConfig.base_file_name = jobName;
    read_config_file();
    
    int[][] taskMap;
    taskMap.length = nTasks;
    int[] taskLoads;
    taskLoads.length = nTasks;

    // Determine compute load per block.
    // We simply use the cell count as an estimate of load.
    Tuple!(int,int)[] blockLoads;
    blockLoads.length = GlobalConfig.nFluidBlocks;
    foreach (iblk, blk; localFluidBlocks)
        blockLoads[iblk] = tuple(to!int(iblk), to!int(globalFluidBlocks[iblk].ncells)); 
    sort!("a[1] > b[1]")(blockLoads);
    
    // Perform load balance and write out.
    distributeBlocksToTasks(nTasks, blockLoads, taskMap, taskLoads);
    string fName = jobName~".mpimap";
    writeBlocksToTasksMap(fName, taskMap, taskLoads);

    writeln("Done.");
}



void distributeBlocksToTasks(int nTasks, Tuple!(int,int)[] blockLoads, int[][] taskMap, int[] taskLoads)
{
    taskLoads.length = nTasks;
    taskLoads[] = 0;

    foreach (bLoad; blockLoads) {
        // The algorithm is:
        // Given a list sorted from heaviest load to lightest,
        // at any stage, assign block to the processing task
        // with the minimum load currently.
        auto blkId = bLoad[0];
        auto blkLoad = bLoad[1];
        auto iMin = taskLoads.minIndex;
        taskLoads[iMin] += blkLoad;
        taskMap[iMin] ~= blkId;
    }

}

void writeBlocksToTasksMap(string fName, int[][] taskMap, int[] taskLoads)
{
    auto nTasks = taskMap.length;
    auto f = File(fName, "w");
    // Write out task map in JSON format
    f.writeln("{");
    f.writefln("\"ntasks\" :  %d,", nTasks);
    foreach (i, task; taskMap) {
        f.writefln("\"task_%d\": {", i);
        f.writefln("    \"load\" : %d,", taskLoads[i]);
        f.writef("    \"blocks\" : [ ");
        foreach (iblk; task[0 .. $-1]) {
            f.writef("%d, ", iblk);
        }
        // Write final entry WITHOUT trailing comma
        f.writefln("%d ]", task[$-1]);
        f.writeln("},");
    }
    f.writeln("\"dummy_entry_without_trailing_comma\" : 0");
    f.writeln("}");
}
