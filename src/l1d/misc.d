// misc.d
// Bits and pieces for the 1D Lagrangian flow solver.
//
// PAJ, 2020-04-09
//
module misc;

import std.conv;
import std.stdio;
import std.string;
import std.format;
import std.algorithm;

import config;

void skip_to_data_at_tindx(File fp, int tindx=0)
{
    string text;
    int my_tindx;
    bool found = false;
    // Let's befin our search for the correct tindx.
    text = fp.readln().chomp();
    while (!found) {
        while (!canFind(text, "tindx")) { text = fp.readln().chomp(); }
        text.formattedRead!"# tindx %d"(my_tindx);
        if (my_tindx == tindx) {
            // Stop the search.
            found = true;
        } else {
            // Continue the search.
            text = fp.readln().chomp();
        }
    }
    // At this point the file pointer should be just after the tindx line.
    return;
} // end skip_to_data_at_tindx()


double get_time_from_times_file(int tindx=0)
{
    // Returns the value of time at tindx.
    string fileName = L1dConfig.job_name ~ "/times.data";
    File fp = File(fileName, "r");
    string text = fp.readln().chomp();
    while (text.canFind("#")) { text = fp.readln().chomp(); }
    string[] items = text.split();
    int myTindx = to!int(items[0]);
    while (myTindx < tindx) {
        text = fp.readln().chomp();
        items = text.split();
        myTindx = to!int(items[0]);
    }
    // We should be at the line that contains the requested tindx.
    return to!double(items[1]);
} // end get_time_from_times_file()


double[int] readTimesFile()
{
    // Returns the associative array of all time values.
    double[int] times;
    string fileName = L1dConfig.job_name ~ "/times.data";
    File fp = File(fileName, "r");
    string txt = fp.readln().chomp(); // Discard header line
    double previous_time = 0.0;
    while (!fp.eof()) {
        txt = fp.readln().chomp();
        if (txt.length > 0) {
            int tindx; double tme;
            txt.formattedRead!"%d %e"(tindx, tme);
            if (tme < previous_time) {
                writeln("Warning: at tindx=%d, time=%e but previous=%e",
                        tindx, tme, previous_time);
            }
            times[tindx] = tme;
            previous_time = tme;
        }
    }
    fp.close();
    return times;
} // end readTimesFile()

