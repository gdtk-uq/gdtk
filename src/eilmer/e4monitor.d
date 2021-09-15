// e4monitor.d
// Program to monitor the progress of an Eilmer job in PBS
// and to terminate that job if progress is not being made.
//
// PJ 2019-09-19, first experiments.
//    2021-06-07, quit after the first timeout and no text.
//    2021-09-16, Let the monitor program run in foreground.
//
// Option 1:
//   From within your PBS batch script, start the monitor program
//   as a background process just before the main run command.
//   For example, in the 2D/sharp-cone-20-degrees/sg-mpi/ example:
//   $ e4monitor --job=cone20 --startup=2 --period=1 &
//   $ mpirun -np 8 --oversubscribe e4mpi --run --job=cone20 --verbosity=1
//
// Option 2:
//   From within your batch script, let the montior program run
//   as a foreground process and start the main run command as a subprocess.
//   For example:
//   $ e4monitor --job=cone20 --startup=1 --period=1 \
//               --command="mpirun -np 8 --oversubscribe e4mpi --run --job=cone20"
//

import std.conv;
import std.stdio;
import std.format;
import std.getopt;
import std.file;
import std.string;
import std.algorithm;
import core.thread;
import std.process;

int main(string[] args)
{
    string jobName = "";
    int startup = 30;
    int period = 60;
    string cmdStr = "";
    //
    string usageMsg = "Usage: e4monitor --job=<string>"
        ~ " [--startup=<int>]"
        ~ " [--period=<int>]"
        ~ " [--command=<string>]";
    if (args.length < 2) {
        writeln("Too few arguments.  Need at least --job=<string>.");
        writeln(usageMsg);
        stdout.flush();
        return 1;
    }
    try {
        getopt(args, "job", &jobName, "startup", &startup, "period", &period, "command", &cmdStr);
    } catch (Exception e) {
        writeln("Problem parsing command-line options.");
        writeln("Arguments not processed: ");
        args = args[1 .. $]; // Dispose of program name in first argument.
        foreach (myarg; args) writeln("    arg: ", myarg);
        writeln(usageMsg);
        stdout.flush();
        return 1;
    }
    //
    writefln("monitor begin, job=\"%s\" startup=%d period=%d command=\"%s\"",
             jobName, startup, period, cmdStr);
    string progressFile = format("config/%s-progress.txt", jobName);
    writefln("monitor file=\"%s\"", progressFile);
    //
    // We may have decided to run the main command in a sub-shell.
    Pid cmdPid;
    if (cmdStr.length > 0) {
        writeln("Starting the main run command as a sub-shell.");
        cmdPid = spawnShell(cmdStr);
        scope(exit) { wait(cmdPid); }
    }
    //
    // Start the monitoring activities.
    int myCount = 0;
    Thread.sleep(dur!("seconds")(startup));
    string content = chomp(readTextWithRetries(progressFile));
    writefln("monitor %d content \"%s\"", myCount, content);
    if (content.length == 0) {
        writeln("Quitting the monitor program; progress text is empty at the first timeout.");
        stdout.flush();
        return 0;
    }
    bool simulationDone = canFind(content, "done");
    int oldStep = 0;
    int newStep = 0;
    int stalledCount = 0;
    if (isNumeric(content)) { oldStep = to!int(content); }
    while (!simulationDone) {
        Thread.sleep(dur!("seconds")(period));
        myCount += 1;
        content = chomp(readTextWithRetries(progressFile));
        writef("monitor %d content \"%s\"", myCount, content);
        simulationDone = canFind(content, "done");
        if (simulationDone) {
            writeln(" simulation done, breaking loop.");
            break;
        }
        //
        if (isNumeric(content)) { newStep = to!int(content); }
        bool progressing = newStep > oldStep;
        if (!progressing) { stalledCount += 1; } else { stalledCount = 0; }
        writefln(" %s", (progressing) ? "progressing" : "stalled");
        stdout.flush();
        //
        if (stalledCount > 10) {
            writeln("monitor delete job");
            auto which_cmd = execute(["which", "qdel"]);
            string pbs_jobid = environment.get("PBS_JOBID");
            if ((which_cmd.status == 0) && (pbs_jobid !is null)) {
                auto qdel_cmd = execute(["qdel", pbs_jobid]);
                if (qdel_cmd.status == 0) {
                    writeln("monitor qdel command succeeded");
                } else {
                    writeln("monitor qdel command failed");
                    writeln(qdel_cmd.output);
                }
            } else {
                writeln("monitor does not know how to delete job");
            }
            stdout.flush();
        }
        oldStep = newStep;
    }
    writeln("monitor done.");
    stdout.flush();
    return 0;
} // end main

string readTextWithRetries(string fileName)
{
    string content = "";
    foreach (i; 0 .. 3) {
        try {
            content = readText(fileName);
        } catch (FileException e) {
            writefln("monitor could not read file \"%s\"", fileName);
            content = "";
        }
        if (content.length > 0) { break; }
        Thread.sleep(dur!("seconds")(1));
    }
    return content;
} // end readTextWithRetries
