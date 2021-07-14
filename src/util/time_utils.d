/**
 * time_utils.d
 * Some routines for converting and handling times
 *
 * @author: Nick Gibbons
 */

module util.time_utils;

import std.string;
import std.format;
import std.conv;

pure int timeStringToSeconds(string timeString) {
/*
    Convert the command line argument max-wall-clock from a string to an amount of seconds.

    Inputs:
     - timeString: may be an amount of seconds or HH:MM:SS format.

    Outputs:
     - The max wall clock time, in seconds, as an integer.

    Notes:
    @author: Nick Gibbons
*/
    int outseconds;
    auto tokens = timeString.split(":");

    if (tokens.length==3) {
        int hours = to!int(tokens[0]);
        int minutes = to!int(tokens[1]);
        int seconds = to!int(tokens[2]);
        outseconds = 3600*hours + 60*minutes + seconds;

    } else if (tokens.length==1) {
        outseconds = to!int(tokens[0]);
    } else {
        throw new Error(format("Invalid maxWallClock setting %s", timeString));
    }
    
    return outseconds;
}

