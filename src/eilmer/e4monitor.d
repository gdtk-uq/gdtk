// e4monitor.d
// Program to monitor the progress of an Eilmer job
// and to terminate that job if progress is not being made.
//
// PJ 2019-09-19, first experiments.
//

import std.conv;
import std.stdio;
import std.format;
import std.getopt;
import std.file;
import std.string;
import std.algorithm;
import core.thread;

int main(string[] args)
{
    string jobName = "";
    int startup = 30;
    int period = 60;
    //
    string usageMsg = "Usage: e4monitor --job=<string>"
        ~ " [--startup=<int>]"
        ~ " [--period=<int>]";
    if (args.length < 2) {
        writeln("Too few arguments.  Need at least --job=<string>.");
        writeln(usageMsg);
        stdout.flush();
        return 1;
    }
    try {
        getopt(args, "job", &jobName, "startup", &startup, "period", &period);
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
    writefln("monitor begin, job=\"%s\" startup=%d period=%d",
             jobName, startup, period);
    string progressFile = format("config/%s-progress.txt", jobName);
    writefln("monitor file=\"%s\"", progressFile);
    int myCount = 0;
    Thread.sleep(dur!("seconds")(startup));
    string content = chomp(readTextWithRetries(progressFile));
    writefln("monitor %d content \"%s\"", myCount, content);
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
            writeln(" breaking loop.");
            break;
        }
        if (isNumeric(content)) { newStep = to!int(content); }
        bool progressing = newStep > oldStep;
        if (!progressing) { stalledCount += 1; } else { stalledCount = 0; }
        writefln(" %s", (progressing) ? "progressing" : "stalled");
        if (stalledCount > 10) {
            writeln("monitor kill job");
            // FIX ME with some real action here.
        }
        oldStep = newStep;
    }
    writeln("monitor done.");
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
}
