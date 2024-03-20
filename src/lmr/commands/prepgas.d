/**
 * Module for calling prep-gas tool as lmr-type command.
 *
 * Authors: RJG
 * Date: 2024-03-09
 */

module lmr.commands.prepgas;

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

Command prepGasCmd;
string cmdName = "prep-gas";

static this()
{
    prepGasCmd.main = &main_;
    prepGasCmd.description = "Prepare gas file for a simulation.";
    prepGasCmd.shortDescription = prepGasCmd.description;
    prepGasCmd.helpMsg = format(
`lmr %s --input=gas-model.lua --output=model.gas

Prepare a gas file for eilmer based on a high-level input file.

required options:
  -i, --input=gasinputfile
      A gas input file, typically a Lua-type input
  -o, --output=gasoutfile
      The name of the output file. This is the file to be used
      by the simulation program. It is the file that appears
      in a setGasModel() call in a simulation input script.

options ([+] can be repeated):
  -v, --verbose [+]
      Increase verbosity during gas file creation.

`, cmdName);
}

int main_(string[] args)
{
    int verbosity = 0;
    string inputFile, outputFile;
    try {
        getopt(args,
               config.bundling,
               "v|verbose+", &verbosity,
               std.getopt.config.required,
               "i|input", &inputFile,
               std.getopt.config.required,
               "o|output", &outputFile);
    }
    catch (Exception e) {
        writefln("Eilmer %s program quitting.", cmdName);
        writeln("There was something wrong with the command-line arguments/options.");
        writeln(e.msg);
        return 1;
    }

    if (inputFile == outputFile) {
        writefln("Eilmer %s program quitting because the output file (%s) would overwrite an input file", cmdName, outputFile);
        writefln("Try changing the name of your output file");
        return 1;
    }

    if (verbosity > 0) {
        writefln("lmr %s: Begin creation of gas file.", cmdName);
    }
    if (verbosity > 1) {
        writefln("Input file: %s", inputFile);
    }

    auto L = luaL_newstate();
    luaL_openlibs(L);

    // arg[1] = inputFile; arg[2] = outputFile
    lua_newtable(L);
    lua_pushstring(L, inputFile.toStringz);
    lua_rawseti(L, -2, 1);
    lua_pushstring(L, outputFile.toStringz);
    lua_rawseti(L, -2, 2);
    lua_setglobal(L, "arg");

    // Now call main()
    if (luaL_dofile(L, toStringz(dirName(thisExePath())~"/prep-gas")) != 0) {
        writeln("There was a problem when trying to use prep-gas.");
        string errMsg = to!string(lua_tostring(L, -1));
        throw new LmrPreProcessingException(errMsg);
    }

    if (verbosity > 1) {
        writefln("Output file: %s", outputFile);
    }

    if (verbosity > 0) {
        writefln("lmr %s: Done creating gas file.", cmdName);
    }
    return 0;
}
