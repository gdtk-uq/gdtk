/**
 * Module for doing preparation of flow fields and simulation parameters.
 *
 * Authors: RJG, PJ, KAD, NNG
 * Date: 2022-08-09
 * History:
 *   2024-02-11 -- Renamed module: prepflow --> prepsim
 *                 to better reflect its more general role
 */

module prepsim;

import core.stdc.stdlib : system;
import std.getopt;
import std.stdio : writeln, writefln;
import std.string : toStringz, strip, split, format;
import std.conv : to;
import std.path : dirName;
import std.file : thisExePath, exists, rmdirRecurse;
import std.json : JSONValue;
import std.algorithm : sort, uniq;

import util.lua;
import util.lua_service : getString;
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

Command prepSimCmd;
string cmdName = "prep-sim";

static this()
{
    prepSimCmd.main = &main_;
    prepSimCmd.description = "Prepare initial flow fields and parameters for a simulation.";
    prepSimCmd.shortDescription = prepSimCmd.description;
    prepSimCmd.helpMsg =
`lmr prep-sim [options]

Prepare an initial flow field and simulation parameters based on a file called job.lua.

options ([+] can be repeated):

 -v, --verbose [+]
     Increase verbosity during preparation of the simulation files.

 -j, --job=flow.lua
     Specify the input file to be something other than job.lua.
     default: job.luare

 --with-cea-usage
     If using a CEAGas as part of the simulation setup, then the debug-version
     of the executable needs to be called. This flag lets the program know
     in advance that cea will be used, and so the debug executable can be called automatically.
     (Alternatively, one could call 'lmr-debug' directly at prep stage.)
`;

}

int main_(string[] args)
{
    string blocksForPrep = "";
    int verbosity = 0;
    string userFlowName = lmrCfg.jobFile;
    bool withCEAGas = false;

    // First only look for with-cea-usage. If we find it, we delegate the rest.
    try {
        getopt(args,
               config.bundling,
               config.passThrough,
               "with-cea-usage", &withCEAGas
               );
    } catch (Exception e) {
        writefln("Eilmer %s program quitting.", cmdName);
        writeln("There is something wrong with the command-line arguments/options.");
        writeln(e.msg);
        return 1;
    }
    if (withCEAGas) {
        string shellStr = "lmr-debug prep-sim";
        if (args.length >= 3) {
            foreach (s; args[2 .. $]) {
                shellStr ~= " " ~ s;
            }
        }
        return system(shellStr.toStringz);
    }

    // From here on we do a regular execution of prep-sim
    try {
        getopt(args,
               config.bundling,
               "v|verbose+", &verbosity,
               "j|job", &userFlowName,
               );
    } catch (Exception e) {
        writefln("Eilmer %s program quitting.", cmdName);
        writeln("There is something wrong with the command-line arguments/options.");
        writeln(e.msg);
        return 1;
    }

    if (verbosity > 1) { writeln("lmr prep-sim: Start lua connection."); }
    auto L = initLuaStateForPrep();
    lua_pushinteger(L, verbosity);
    lua_setglobal(L, "verbosity");

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
    if (luaL_dofile(L, toStringz(dirName(thisExePath())~"/../lib/prepsim.lua")) != 0) {
        writeln("There was a problem in the Eilmer Lua code: prepsim.lua");
        string errMsg = to!string(lua_tostring(L, -1));
        throw new FlowSolverException(errMsg);
    }
    if (luaL_dostring(L, toStringz("readGridMetadata()")) != 0) {
        writeln("There was a problem in the Eilmer build function readGridMetadata() in prepsim.lua");
        string errMsg = to!string(lua_tostring(L, -1));
        throw new FlowSolverException(errMsg);
    }
    //
    // We are ready for the user's input script.
    if (!exists(userFlowName)) {
        writefln("The file %s does not seems to exist.", userFlowName);
        writeln("Did you mean to specify a different job name?");
        return 1;
    }
    if (luaL_dofile(L, toStringz(userFlowName)) != 0) {
        writeln("There was a problem in the user-supplied input lua script: ", userFlowName);
        string errMsg = to!string(lua_tostring(L, -1));
        throw new FlowSolverException(errMsg);
    }
    //
    if (luaL_dostring(L, toStringz("buildRuntimeConfigFiles()")) != 0) {
        writeln("There was a problem in the Eilmer build function buildRuntimeConfigFiles() in prepsim.lua");
        string errMsg = to!string(lua_tostring(L, -1));
        throw new FlowSolverException(errMsg);
    }
    JSONValue jsonData = readJSONfile(lmrCfg.cfgFile);
    set_config_for_core(jsonData);
    // We may not proceed to building of block files if the config parameters are incompatible.
    checkGlobalConfig();

    // Clean out anything in snapshots area if present from earlier run
    if (lmrCfg.snapshotDir.exists) {
        if (verbosity > 1) { writeln("lmr prep-sim: Removing old snapshots."); }
        lmrCfg.snapshotDir.rmdirRecurse;
    }
    // and now we're ready to build new stuff!
    if (luaL_dostring(L, toStringz("buildGridAndFieldFiles()")) != 0) {
        writeln("There was a problem in the Eilmer build function buildGridAndFieldFiles() in prepsim.lua");
        string errMsg = to!string(lua_tostring(L, -1));
        throw new FlowSolverException(errMsg);
    }
    if (verbosity > 0) { writeln("lmr prep-sim: Done."); }

    return 0;
}
