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


string make_path_name(string mytype)(int tindx)
// Build a pathname for "grid" and "flow" subdirectories.
// These directories are "labelled" with time index.
{
    auto writer = appender!string();
    formattedWrite(writer, "%s/t%04d", mytype, tindx);
    return writer.data;
}

string make_file_name(string mytype)(string base_file_name, int blk_id, int tindx)
// Build a pathname for "grid" and "flow" files which are stored in
// subdirectories and labelled with the block id and time index counters.
{
    auto writer = appender!string();
    formattedWrite(writer, "%s/t%04d/%s.%s.b%04d.t%04d.gz",
		   mytype, tindx, base_file_name,
		   mytype, blk_id, tindx);
    return writer.data;
}

void ensure_directory_is_present(string dir_name)
{
    if (exists(dir_name) && isDir(dir_name)) return;
    try {
	mkdirRecurse(dir_name);
    } catch (FileException e) {
	throw new Error(text("Failed to ensure directory is present: ", dir_name));
    }
}
