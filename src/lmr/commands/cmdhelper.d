/**
 * This module houses some frequently used service functions
 * for the various commands. Hence it is a helper module.
 *
 * Author: Rowan J. Gollan
 * Date: 2023-07-14
 */

module cmdhelper;

import std.stdio;
import std.file;
import std.string;
import std.format : format;
import std.algorithm;
import std.range : array;
import std.regex;
import std.path;
import std.conv : to;

import dyaml;

import lmrconfig : lmrCfg;

string[] determineAvailableSnapshots(string dir="")
{
    string snapshotDir = (dir == "") ? lmrCfg.snapshotDir : dir;

    string[] availSnapshots;
    // Build list of all available snapshots
    auto dirs = dirEntries(snapshotDir, SpanMode.shallow).
        filter!(a => a.isDir()).
        map!(a => a.name.baseName)().
        array();
    sort(dirs);
    /*
     * Keep up to date with lmr.cfg::snapshotIdxFmrStr
     */
    auto fourDigits = regex(r"[0-9][0-9][0-9][0-9]");
    foreach (string d; dirs) {
        auto rgxMtch = matchAll(d, fourDigits);
        if (!rgxMtch.empty()) {
            availSnapshots ~= d;
        }
    }
    return availSnapshots;
}

/**
 * Determine which snapshots to process.
 *
 * The first parameter gives a list of available snapshots.
 * The remaining parameters would typically be picked up by the command line parser.
 * If the parameters indicate no selection, then default to selecting the FINAL snapshot.
 *
 * Author: Rowan J. Gollan
 * Date: 2023-07-14
 */

string[] determineSnapshotsToProcess(string[] availSnapshots, int[] snapshots,
    bool allSnapshots, bool finalSnapshot)
{
    string[] snaps;
    if (allSnapshots) {
        snaps.length = availSnapshots.length;
        snaps[] = availSnapshots[];
    }
    if (snapshots && !allSnapshots) {
        foreach (i; snapshots) snaps ~= format(lmrCfg.snapshotIdxFmt, i);
    }
    if (finalSnapshot && !allSnapshots) {
        snaps ~= availSnapshots[$-1];
    }
    if (snaps.length == 0) { // We've picked up nothing so far from the previous options,
                             // so do the final snapshot.
        snaps ~= availSnapshots[$-1];
    }
    return uniq(sort(snaps)).array;
}

/**
 * Find timesteps associated with snapshots from times file.
 *
 * Author: Rowan J. Gollan
 * Date: 2024-02-24
 */
auto mapTimesToSnapshots(string[] snaps)
{
    double[] times;
    if (snaps.length == 1 && snaps[0] == format(lmrCfg.snapshotIdxFmt, 0)) {
        // This means we only have an initial snapshot
        // and no times file yet.
        times ~= 0.0;
        return times;
    }
    Node timesData = dyaml.Loader.fromFile(lmrCfg.timesFile).load();
    foreach (snap; snaps) {
        Node snapData = timesData[snap];
        times ~= snapData["time"].as!double;
    }
    return times;
}

/**
 * For structured block processing, decode user-supplied range indices.
 *
 * Author: PAJ
 * Date: 2024-05-22 (moved into this module)
 */
size_t[] decode_range_indices(string rangeStr, size_t first, size_t endplus1)
// Decode strings such as "0:$", ":", "0:3", "$"
// On input, first and endplus1 represent the largest, available range.
// Return the pair of numbers that can be used in a foreach loop range.
{
    if (rangeStr == ":") {
        return [first, endplus1];
    }
    if (canFind(rangeStr, ":")) {
        // We have a range specification to pull apart.
        auto items = rangeStr.split(":");
        first = to!size_t(items[0]);
        if (items.length > 1 && items[1] != "$") {
            // Presume that we have a second integer.
            size_t new_endplus1 = to!size_t(items[1]);
            if (new_endplus1 < endplus1) endplus1 = new_endplus1;
        }
    } else if (rangeStr == "$") {
        // With just a single "$" specified, we want only the last index.
        first = endplus1 - 1;
    }else {
        // Presume that we have a single integer.
        first = to!size_t(rangeStr);
        if (first < endplus1) endplus1 = first+1;
    }
    return [first, endplus1];
} // end decode_range_indices()
