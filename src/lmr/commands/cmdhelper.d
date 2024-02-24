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
    Node timesData = dyaml.Loader.fromFile(lmrCfg.timesFile).load();
    foreach (snap; snaps) {
        Node snapData = timesData[snap];
        times ~= snapData["time"].as!double;
    }
    return times;
}


