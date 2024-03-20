/**
 * Module for calling prep-chem tool as lmr-type command.
 *
 * I'm going to call this the prep-reactions command because
 * it prepares a "reactions_file" as known in Eilmer config.
 * We can have an alias to prep-chem.
 *
 * Authors: RJG
 * Date: 2024-03-10
 */

module lmr.commands.prepreactions;

import std.getopt;
import std.stdio : writeln, writefln;
import std.string : toStringz;
import std.conv : to;
import std.path : dirName;
import std.file : thisExePath, exists;
import std.format : format;
version(macosx) {
    import core.stdc.stdlib : system;
}

import util.lua;

import command;
import lmrexceptions : LmrPreProcessingException;

Command prepReacCmd;
string cmdName = "prep-reactions";

static this()
{
    prepReacCmd.main = &main_;
    prepReacCmd.description = "Prepare a reactions file for a simulation.";
    prepReacCmd.shortDescription = prepReacCmd.description;
    prepReacCmd.helpMsg = format(
`lmr %s [options] --gasfile=model.gas --input=reactions.lua --output=reactions.chem

Prepare a reactions file for eilmer based on a high-level input file.

required options:
  -g, --gasfile=gasmodelfile
      where 'gasmodelfile' has been prepared with the prep-gas command.
      (see: lmr help prep-gas)

  -i, --input=reacinputfile
      where 'reacinputfile' is a Lua-type input that specifies the
      chemical reactions and the parameters for computing rate constants.
      The species that appear in 'reacinputfile' should be consistent
      with the set provided in 'gasmodelfile'.

  -o, --output=reacoutputfile
      where 'reacoutputfile' is the name of the output file. This is the file to be used
      by the simulation program. It is the file that appears when setting 'config.reactions_file'
      in an Eilmer input script.

options ([+] can be repeated):
  -v, --verbose [+]
      Increase verbosity during gas file creation.
  -c, --compact
      Also produce a compact form of the reactions file called 'chem-compact-notation.inp'

For more information on how to construct a 'reacinputfile', see: 

      https://gdtk.uqcloud.net/pdfs/reacting-gas-guide.pdf

`, cmdName);
}

int main_(string[] args)
{
    int verbosity = 0;
    string gasFile, inputFile, outputFile;
    bool withCompact = false;
    try {
        getopt(args,
               config.bundling,
               "v|verbose+", &verbosity,
               "c|compact", &withCompact,
               std.getopt.config.required,
               "g|gasfile", &gasFile,
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

    if (outputFile == inputFile || outputFile == gasFile) {
        writefln("Eilmer %s program quitting because the output file (%s) would overwrite an input file", cmdName, outputFile);
        writefln("Try changing the name of your output file");
        return 1;
    }

    if (verbosity > 0) {
        writefln("lmr %s: Begin creation of reactions file.", cmdName);
    }
    if (verbosity > 1) {
        writefln("Gas file: %s", gasFile);
        writefln("Reactions input file: %s", inputFile);
    }

    version(macosx) {
        // RJG, 2024-03-10
        // On macosx, I can't get lmr program to play nicely with lpeg.so
        // So we'll do raw system call to the Lua program prep-chem
        string cmd = "prep-chem ";
        if (withCompact) cmd ~= "--compact ";
        cmd ~= format("%s %s %s", gasFile, inputFile, outputFile);
        auto flag = system(cmd.toStringz);
        if (flag != 0) return flag;
    }
    else {
        // On linux, we can execute the Lua environment from this program
        auto L = luaL_newstate();
        luaL_openlibs(L);

        // Set arg table for prep-chem.
        lua_newtable(L);
        int index = 1;
        if (withCompact) {
            lua_pushstring(L, "--compact");
            lua_rawseti(L, -2, index++);
        }
        lua_pushstring(L, gasFile.toStringz);
        lua_rawseti(L, -2, index++);
        lua_pushstring(L, inputFile.toStringz);
        lua_rawseti(L, -2, index++);
        lua_pushstring(L, outputFile.toStringz);
        lua_rawseti(L, -2, index++);
        lua_setglobal(L, "arg");

        // Now call main()
        if (luaL_dofile(L, toStringz(dirName(thisExePath())~"/prep-chem")) != 0) {
            writeln("There was a problem when trying to use the Lua program prep-chem.");
            string errMsg = to!string(lua_tostring(L, -1));
            throw new LmrPreProcessingException(errMsg);
        }
    } // end else-version

    if (verbosity > 1) {
        writefln("Reactions output file: %s", outputFile);
    }

    if (verbosity > 0) {
        writefln("lmr %s: Done creating reactions file.", cmdName);
    }
    return 0;
}
