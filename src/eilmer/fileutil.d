// fileutil.d
// A small set of functions to handle the making of filenames and directories.
//
// Extracted from e4_core.d 2015-02-28 so we can reuse them.

import std.stdio;
import std.file;
import std.conv;
import std.array;
import std.format;
import std.string;
import globalconfig: FlowSolverException;
import core.thread;

string make_path_name(string mytype)(int tindx)
// Build a pathname for "grid" and "flow" subdirectories.
// These directories are "labelled" with time index.
{
    auto writer = appender!string();
    formattedWrite(writer, "%s/t%04d", mytype, tindx);
    return writer.data;
}

string make_path_name(string mytype, int tindx)
{
    auto writer = appender!string();
    formattedWrite(writer, "%s/t%04d", mytype, tindx);
    return writer.data;
}

string make_file_name(string mytype)(string base_file_name, int blk_id, int tindx, string myext)
// Build a pathname for "grid" and "flow" files which are stored in
// subdirectories and labelled with the block id and time index counters.
{
    auto writer = appender!string();
    formattedWrite(writer, "%s/t%04d/%s.%s.b%04d.t%04d.%s",
                   mytype, tindx, base_file_name,
                   mytype, blk_id, tindx, myext);
    return writer.data;
}

string make_file_name(string folder, string mytype, string base_file_name, int blk_id, int tindx, string myext)
// Build a pathname for "grid" and "flow" files which are stored in
// subdirectories and labelled with the block id and time index counters.
{
    auto writer = appender!string();
    formattedWrite(writer, "%s/%s/t%04d/%s.%s.b%04d.t%04d.%s",
                   folder, mytype, tindx, base_file_name,
                   mytype, blk_id, tindx, myext);
    return writer.data;
}

string make_snapshot_path_name(string type, int idx)
{
    return format!"%s/snapshot-%04d"(type, idx);
}

string make_snapshot_file_name(string type, string base_file_name, int blk_id, int idx, string ext)
{
    return format!"%s/snapshot-%04d/%s.%s.b%04d.snapshot.%04d.%s"(type, idx, base_file_name, type, blk_id, idx, ext);
}

string make_DFT_file_name(string base_file_name, int blk_id, string myext) {
// Based on make_file_name, but specifically for building the DFT file. This has no time indices, and no block subdirs.
    
    auto writer = appender!string();
    formattedWrite(writer, "DFT/%s.DFT.b%04d.%s",
                    base_file_name, blk_id, myext);
    return writer.data;
}

void ensure_directory_is_present(string dir_name)
{
    if (exists(dir_name) && isDir(dir_name)) return;
    try {
        mkdirRecurse(dir_name);
    } catch (FileException e) {
        string msg = text("Failed to ensure directory is present: ", dir_name);
        throw new FlowSolverException(msg);
    }
}

void wait_for_directory_to_be_present(string dir_name)
{
    int nChecks = 30;
    int waitTime = 100; // msecs
    foreach (i; 0 .. nChecks) {
        if (exists(dir_name) && isDir(dir_name)) return;
        Thread.sleep(dur!("msecs")( waitTime ));
    }
    string msg = text("Failed to ensure directory is present: ", dir_name);
    throw new FlowSolverException(msg);
}
