/**
 * Module for doing preparation of flow fields.
 *
 * Authors: RJG, PJ, KAD, NNG
 * Date: 2022-08-09
 */

module prepflow;

import std.getopt;
import std.stdio : writeln, writefln;
import std.string : toStringz, strip, split, format;
import std.conv : to;
import std.path : dirName;
import std.file : thisExePath;
import std.json : JSONValue;
import std.algorithm : sort, uniq;

import util.lua;
import geom.luawrap;
import gas;
import gas.luagas_model;
import nm.luabbla;

import json_helper;
import lua_helper : initLuaStateForPrep;
import lmrconfig : lmrCfg;
import command;
import globalconfig;
import luaflowsolution;
import luaflowstate;
import luaidealgasflow;
import luagasflow;
import blockio : luafn_writeFlowMetadata, luafn_writeInitialFlowFile;

Command prepFlowCmd;

static this()
{
    prepFlowCmd.main = &main_;
    prepFlowCmd.description = "Prepare initial flow fields for an Eilmer simulation.";
    prepFlowCmd.shortDescription = prepFlowCmd.description;
    prepFlowCmd.helpMsg =
`lmr prep-flow [options]

Prepare an initial flow field based on a job called lmr-flow.lua.

options ([+] can be repeated):

 -v, --verbose [+]
     Increase verbosity during grid generation.

 -m, --mode=<mode>
     Select mode for file preparation. Mode can be "steady" or "transient".

`;

}

enum Mode { steady, transient };

void main_(string[] args)
{
    string blocksForPrep = "";
    int verbosity = 0;
    Mode mode;
    getopt(args,
           config.bundling,
           "v|verbose+", &verbosity,
           "mode", &mode,
           );

    if (mode == Mode.transient) {
        writeln("lmr prep-flow: You have selected transient mode. Unfortunately, that's not yet implemented.");
        writeln("lmr prep-flow: Exiting.");
        return;
    }

    if (verbosity > 0) {
        writeln("lmr prep-flow: Begin preparation of initial flow field files.");
        if (mode == Mode.steady) {
            writeln("lmr prep-flow: mode=steady selected.");
        }
    }

    if (verbosity > 1) { writeln("lmr prep-flow: Start lua connection."); }

    auto L = initLuaStateForPrep();
    lua_pushinteger(L, verbosity);
    lua_setglobal(L, "verbosity");
    // RJG, 2023-06-27
    // Add a few more lua-wrapped functions for use in prep.
    // These functions are not backported into Eilmer 4, and
    // I don't want to hijack initLuaStateForPrep() just yet.
    // At some point in the future, this can be handled inside
    // initLuaStateForPrep().
    lua_pushcfunction(L, &luafn_writeFlowMetadata);
    lua_setglobal(L, "writeFlowMetadata");
    lua_pushcfunction(L, &luafn_writeInitialFlowFile);
    lua_setglobal(L, "writeInitialFlowFile");

    // Determine which fluidBlocks we need to process.
    int[] blockIdList;
    blocksForPrep = blocksForPrep.strip();
    foreach (blkStr; blocksForPrep.split(",")) {
        blkStr = blkStr.strip();
        auto blkRange = blkStr.split("..<");
        if (blkRange.length == 1) {
            blockIdList ~= to!int(blkRange[0]);
        }
        else if (blkRange.length == 2) {
            auto start = to!int(blkRange[0]);
            auto end = to!int(blkRange[1]);
            if (end < start) {
                string errMsg = "Supplied block list is in error. Range given is not allowed.";
                errMsg ~= format("Bad supplied range is: %s", blkStr);
                throw new UserInputError(errMsg);
            }
            foreach (i; start .. end) {
                blockIdList ~= i;
            }
        }
        else {
            string errMsg = "Supplied block list is in error. Range given is not allowed.";
            errMsg ~= format("Bad supplied range is: %s", blkStr);
            throw new UserInputError(errMsg);
        }
    }
    // Let's sort blocks in ascending order
    blockIdList.sort();
    lua_newtable(L);
    lua_setglobal(L, "fluidBlockIdsForPrep");
    lua_getglobal(L, "fluidBlockIdsForPrep");
    // Use uniq so that we remove any duplicates the user might have supplied
    import std.range;
    foreach (i, blkId; blockIdList.uniq().enumerate(1)) {
        lua_pushinteger(L, blkId);
        lua_rawseti(L, -2, to!int(i));
    }
    lua_pop(L, 1);
    // Now that we have set the Lua interpreter context,
    // process the Lua scripts.
    if (luaL_dofile(L, toStringz(dirName(thisExePath())~"/../lib/prepflow.lua")) != 0) {
        writeln("There was a problem in the Eilmer Lua code: prepflow.lua");
        string errMsg = to!string(lua_tostring(L, -1));
        throw new FlowSolverException(errMsg);
    }
    if (luaL_dostring(L, toStringz("readGridMetadata()")) != 0) {
        writeln("There was a problem in the Eilmer build function readGridMetadata() in prepflow.lua");
        string errMsg = to!string(lua_tostring(L, -1));
        throw new FlowSolverException(errMsg);
    }
    // We are ready for the user's input script.
    string userFlowName = lmrCfg.jobFile;
    if (luaL_dofile(L, toStringz(userFlowName)) != 0) {
        writeln("There was a problem in the user-supplied input lua script: ", userFlowName);
        string errMsg = to!string(lua_tostring(L, -1));
        throw new FlowSolverException(errMsg);
    }
    if (luaL_dostring(L, toStringz("buildRuntimeConfigFiles()")) != 0) {
        writeln("There was a problem in the Eilmer build function buildRuntimeConfigFiles() in prepflow.lua");
        string errMsg = to!string(lua_tostring(L, -1));
        throw new FlowSolverException(errMsg);
    }
    JSONValue jsonData = readJSONfile(lmrCfg.cfgFile);
    set_config_for_core(jsonData);
    // We may not proceed to building of block files if the config parameters are incompatible.
    checkGlobalConfig();
    if (luaL_dostring(L, toStringz("buildFlowAndGridFiles()")) != 0) {
        writeln("There was a problem in the Eilmer build function buildFlowAndGridFiles() in prepflow.lua");
        string errMsg = to!string(lua_tostring(L, -1));
        throw new FlowSolverException(errMsg);
    }
    if (verbosity > 0) { writeln("lmr prep-flow: Done."); }

    return;
}


