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
