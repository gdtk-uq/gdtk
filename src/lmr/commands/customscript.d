/**
 * Module for doing miscelleneous tasks, specified in a lua script.
 *
 * Authors: RJG, PJ, KAD, NNG
 * Date: 2022-08-09
 * History:
 *   2024-02-26 -- Ported from eilmer4 by Nick
 */

module customscript;

import std.getopt;
import std.stdio : writeln, writefln;
import std.string : toStringz, strip, split, format;
import std.conv : to;
import std.path : dirName;
import std.file : thisExePath, exists;
import std.json : JSONValue;
import std.algorithm : sort, uniq;

import util.lua;
import util.lua_service : getString;
import geom.luawrap;
import gas;
import gas.luagas_model;
import nm.luabbla;

import util.json_helper;
import lua_helper : initLuaStateForPrep;
import lmrconfig : lmrCfg;
import command;
import globalconfig;
import luaflowsolution;
import luaflowstate;
import gasdyn.luaidealgasflow;
import gasdyn.luagasflow;
import blockio : luafn_writeFluidMetadata, luafn_writeInitialFluidFile;

Command customScriptCmd;
string cmdName = "custom-script";

static this()
{
    customScriptCmd.main = &main_;
    customScriptCmd.description = "Run a user configured lua script.";
    customScriptCmd.shortDescription = customScriptCmd.description;
    customScriptCmd.helpMsg =
`lmr custom-script [options]

Prepare an initial flow field and simulation parameters based on a file called job.lua.

options ([+] can be repeated):

 -v, --verbose [+]
     Increase verbosity during preparation of the simulation files.

 -j, --job=script.lua
     Specify the input file to be something other than job.lua.
     default: job.lua
`;

}

int main_(string[] args)
{
    int verbosity = 0;
    string scriptFile = "script.lua";
    try {
        getopt(args,
               config.bundling,
               "v|verbose+", &verbosity,
               "j|job", &scriptFile,
               );
    } catch (Exception e) {
        writefln("Eilmer %s program quitting.", cmdName);
        writeln("There is something wrong with the command-line arguments/options.");
        writeln(e.msg);
        return 1;
    }

    if (verbosity > 1) { writeln("lmr custom-script: Start lua connection."); }

    auto L = initLuaStateForPrep();
    lua_pushinteger(L, verbosity);
    lua_setglobal(L, "verbosity");

    // We are ready for the user's input script.
    if (!exists(scriptFile)) {
        writefln("The file %s does not seems to exist.", scriptFile);
        writeln("Did you mean to specify a different file name?");
        return 1;
    }

    if (luaL_dofile(L, toStringz(scriptFile)) != 0) {
        writeln("There was a problem in the user-supplied Lua script: ", scriptFile);
        string errMsg = to!string(lua_tostring(L, -1));
        throw new FlowSolverException(errMsg);
    }

    if (verbosity > 0) { writeln("lmr custom-script: Done."); }

    return 0;
}


