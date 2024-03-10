/**
 * Module for calling prep-gas tool with the --list-available-species option.
 *
 * Note: This command has been separated from prepgas.d despite the fact
 * that both commands call the lua 'prep-gas' tool. The reason is that
 * for prepgas.d I wanted to ensure that an input and output file had
 * been supplied. If I had added the "--list-species" option to that
 * command, users would have had to supply dummy inputs and outputs.
 * That would get clunky and unintuitive.
 *
 * Authors: RJG
 * Date: 2024-03-10
 */

module lmr.commands.listspecies;

import std.getopt;
import std.stdio : writeln, writefln;
import std.string : toStringz;
import std.conv : to;
import std.path : dirName;
import std.file : thisExePath, exists;
import std.format : format;

import util.lua;

import command;
import lmrexceptions : LmrPreProcessingException;

Command listSpeciesCmd;
string cmdName = "list-species";

static this()
{
    listSpeciesCmd.main = &main_;
    listSpeciesCmd.description = "List available species from the database.";
    listSpeciesCmd.shortDescription = listSpeciesCmd.description;
    listSpeciesCmd.helpMsg = format(
`lmr %s

List all available species contained in the database.

options ([+] can be repeated):
  -v, --verbose [+]
      Increase verbosity.

`, cmdName);
}

int main_(string[] args)
{
    int verbosity = 0;
    string inputFile, outputFile;
    try {
        getopt(args,
               config.bundling,
               "v|verbose+", &verbosity);
    }
    catch (Exception e) {
        writefln("Eilmer %s program quitting.", cmdName);
        writeln("There was something wrong with the command-line arguments/options.");
        writeln(e.msg);
        return 1;
    }

    if (verbosity > 0) {
        writefln("lmr %s: Begin listing species.", cmdName);
    }

    auto L = luaL_newstate();
    luaL_openlibs(L);

    // arg[1] = "--list-available-species"
    lua_newtable(L);
    lua_pushstring(L, "--list-available-species");
    lua_rawseti(L, -2, 1);
    lua_setglobal(L, "arg");

    // Now call main()
    if (luaL_dofile(L, toStringz(dirName(thisExePath())~"/prep-gas")) != 0) {
        writeln("There was a problem when trying to use prep-gas.");
        string errMsg = to!string(lua_tostring(L, -1));
        throw new LmrPreProcessingException(errMsg);
    }

    if (verbosity > 0) {
        writefln("lmr %s: Done listing species.", cmdName);
    }
    return 0;
}
