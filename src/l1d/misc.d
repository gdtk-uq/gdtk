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
