/**
 * Module to deal with errors and exiting.
 *
 * Authors: RJG
 * Date: 2024-04-27
 */

module lmr.lmrerrors;

import core.stdc.stdlib : exit;
import std.stdio : writeln, writefln;

enum LmrError {
    initialisation = 1,
    inputOutput = 2,
    preprocessing = 3,
    unrecoverableUpdate = 4,
    hitNumericalRailGuard = 5,
    commandLineError = 6,
}

void lmrErrorExit(LmrError err, string msg)
{
    writeln("Exiting lmr program due to error.");
    writeln("Error message is:");
    writefln("%s", msg);
    exit(err);
}

