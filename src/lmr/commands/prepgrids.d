/**
 * Module for doing preparation of grids.
 *
 * Authors: RJG, PJ, KAD, NNG
 * Date: 2022-08-08
 */

module prepgrids;

import std.getopt;
import std.stdio : writeln, writefln;
import std.string : toStringz;
import std.conv : to;
import std.path : dirName;
import std.file : thisExePath;

import util.lua;
import geom.luawrap;

import lua_helper : initLuaStateForPrep;
import lmrconfig : lmrCfg;
import command;
import globalconfig;
import luaflowsolution;
import luaflowstate;

Command prepGridCmd;

static this()
{
    prepGridCmd.main = &main_;
    prepGridCmd.description = "Prepare grids for an Eilmer simulation.";
    prepGridCmd.shortDescription = prepGridCmd.description;
    prepGridCmd.helpMsg =
`lmr prep-grids [options]

Prepare a grid based on a job called lmr-grid.lua.

options ([+] can be repeated):

 -v, --verbose [+]
     Increase verbosity during grid generation.

`;

}


void main_(string[] args)
{
    int verbosity = 0;
    getopt(args,
           config.bundling,
           "v|verbose+", &verbosity
           );

    if (verbosity > 0) {
        writeln("lmr prep-grids: Begin preparation of grid files.");
    }

    if (verbosity > 1) writeln("lmr prep-grids: Start lua connection.");
    auto L = initLuaStateForPrep();
    lua_pushinteger(L, verbosity);
    lua_setglobal(L, "verbosity");
    // Now that we have set the Lua interpreter context,
    // process the Lua scripts.
    if ( luaL_dofile(L, toStringz(dirName(thisExePath())~"/../lib/prepgrids.lua")) != 0 ) {
        writeln("There was a problem in the Eilmer Lua code: prepgrids.lua");
        string errMsg = to!string(lua_tostring(L, -1));
        throw new FlowSolverException(errMsg);
    }
    string userGridName = lmrCfg.jobFile;
    if ( luaL_dofile(L, toStringz(userGridName)) != 0 ) {
        writeln("There was a problem in the user-supplied input lua script: ", userGridName);
        string errMsg = to!string(lua_tostring(L, -1));
        throw new FlowSolverException(errMsg);
    }
    if ( luaL_dostring(L, "writeGridFiles()") != 0 ) {
        writeln("There was a problem in the Lua function writeGridFiles() in prepgrids.lua");
        string errMsg = to!string(lua_tostring(L, -1));
        throw new FlowSolverException(errMsg);
    }
    if (verbosity > 0) { writeln("lmr prep-grids: Done."); }
    return;
}
