/**
 * Module for calling prep-kinetics tool as lmr-type command.
 *
 * I'm going to call this the prep-energy-exchange command because
 * it prepares a "energy_exchange_file" as known in Eilmer config.
 * We can have an alias to prep-kinetics.
 *
 * Authors: RJG
 * Date: 2024-03-10
 */

module lmr.commands.prepenergyexchange;

import std.getopt;
import std.stdio : writeln, writefln;
import std.string : toStringz;
import std.conv : to;
import std.path : dirName;
import std.file : thisExePath, exists;
import std.format : format;
import std.range : empty;
version(macosx) {
    import core.stdc.stdlib : system;
}

import util.lua;

import command;
import lmrexceptions : LmrPreProcessingException;

Command prepExchCmd;
string cmdName = "prep-energy-exchange";

static this()
{
    prepExchCmd.main = &main_;
    prepExchCmd.description = "Prepare an energy exchange file for a simulation.";
    prepExchCmd.shortDescription = prepExchCmd.description;
    prepExchCmd.helpMsg = format(
`lmr %s [options] --gasfile=model.gas --input=energy-exchange.lua --output=mechanisms.exch [--reacfile=reactions.chem]

Prepare an energy exchange file for eilmer based on a high-level input file.

required options:
  -g, --gasfile=gasmodelfile
      where 'gasmodelfile' has been prepared with the prep-gas command.
      (see: lmr help prep-gas)

  -i, --input=eeinputfile
      where 'eeinputfile' is a Lua-type input that specifies the
      energy exchange mechanisms and the parameters for computing relaxation times.
      The species and energy modes that appear in 'eeinputfile' should be consistent
      with the set provided in 'gasmodelfile'.

  -o, --output=eeoutputfile
      where 'eeoutputfile' is the name of the output file. This is the file to be used
      by the simulation program. It is the file that appears when setting 'config.energy_exchange_file'
      in an Eilmer input script.

options ([+] can be repeated):
  -v, --verbose [+]
      Increase verbosity during gas file creation.
  -r, --reacfile=reacfile
      where 'reacfile' has been prepared with the prep-reactions command (see: lmr help prep-reactions).
      A 'reacfile' is required when there is chemistry coupling specified in the 'eeinputfile'.

`, cmdName);
}

int main_(string[] args)
{
    int verbosity = 0;
    string gasFile, reacFile, inputFile, outputFile;
    bool withCompact = false;
    try {
        getopt(args,
               config.bundling,
               "v|verbose+", &verbosity,
               "r|reacfile", &reacFile, 
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

    if (outputFile == inputFile || outputFile == gasFile || outputFile == reacFile ) {
        writefln("Eilmer %s program quitting because the output file (%s) would overwrite an input file", cmdName, outputFile);
        writefln("Try changing the name of your output file");
        return 1;
    }

    if (verbosity > 0) {
        writefln("lmr %s: Begin creation of energy exchange file.", cmdName);
    }
    if (verbosity > 1) {
        writefln("Gas file: %s", gasFile);
        if (!reacFile.empty) writefln("Reactions file: %s", reacFile);
        writefln("Energy exchange input file: %s", inputFile);
    }

    version(macosx) {
        // RJG, 2024-03-10
        // On macosx, I can't get lmr program to play nicely with lpeg.so
        // So we'll do raw system call to the Lua program prep-kinetics
        string cmd = format("prep-kinetics %s ", gasFile);
        if (!reacFile.empty) cmd ~= format("%s ", reacFile);
        cmd ~= format("%s %s", inputFile, outputFile);
        auto flag = system(cmd.toStringz);
        if (flag != 0) return flag;
    }
    else {
        auto L = luaL_newstate();
        luaL_openlibs(L);

        // Set arg table for prep-chem.
        lua_newtable(L);
        int index = 1;
        lua_pushstring(L, gasFile.toStringz);
        lua_rawseti(L, -2, index++);
        if (!reacFile.empty) {
            lua_pushstring(L, reacFile.toStringz);
            lua_rawseti(L, -2, index++);
        }
        lua_pushstring(L, inputFile.toStringz);
        lua_rawseti(L, -2, index++);
        lua_pushstring(L, outputFile.toStringz);
        lua_rawseti(L, -2, index++);
        lua_setglobal(L, "arg");

        // Now call main()
        if (luaL_dofile(L, toStringz(dirName(thisExePath())~"/prep-kinetics")) != 0) {
            writeln("There was a problem when trying to use the Lua program prep-kinetics.");
            string errMsg = to!string(lua_tostring(L, -1));
            throw new LmrPreProcessingException(errMsg);
        }
    }
    if (verbosity > 1) {
        writefln("Energy exchange output file: %s", outputFile);
    }

    if (verbosity > 0) {
        writefln("lmr %s: Done creating energy exchange file.", cmdName);
    }
    return 0;
}
