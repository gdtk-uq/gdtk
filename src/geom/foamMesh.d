/**
 * foamMesh.d
 *
 * Authors: Rowan G., Ingo J. and Peter J.
 * Date: 2017-07-11
 *
 */
import core.stdc.stdlib : exit;
import std.stdio;
import std.getopt;
import std.conv;
import std.string;
import std.path;
import std.file;

import util.lua;

import geom;
import geom.luawrap;

class FoamMeshException : Exception {
    @nogc
    this(string message, string file=__FILE__, size_t line=__LINE__,
         Throwable next=null)
    {
        super(message, file, line, next);
    }
}


int main(string[] args)
{
    int exitFlag = 0;

    string msg = "Usage:                               Comment:\n";
    msg       ~= "foamMesh [--job=<string>]            Job/filename (without .lua extension)\n";
    msg       ~= "         [--verbosity=<int>]         defaults to 1\n";
    msg       ~= "         [--help]                    writes this message\n";
    msg       ~= "\n";

    if (args.length < 2) {
        writeln("Too few arguments.");
        write(msg);
        exitFlag = 1;
        return exitFlag;
    }

    string jobName = "";
    int verbosityLevel = 1; // default to having a little information printed
    bool helpWanted = false;
    try {
        getopt(args,
               "job", &jobName,
               "verbosity", &verbosityLevel,
               "help", &helpWanted
               );
    }
    catch (Exception e) {
        writeln("Problem parsing command-line options.");
        writeln("Arguments not processed: ");
        args = args[1 .. $];
        foreach (arg; args) writeln("    arg: ", arg);
        write(msg);
        exitFlag = 1;
        return exitFlag;
    }
    if (verbosityLevel > 0) {
        writeln("foamMesh grid generator.");
        writeln("Revision: PUT_REVISION_STRING_HERE");
    }
    if (helpWanted) {
        write(msg);
        exitFlag = 0;
        return exitFlag;
    }

    // Check that file actually exists in directory.
    auto fileName = jobName~".lua";
    if (!exists(fileName)) {
        writeln(format("ERROR: The job file '%s' does not appear to exist in the working directory.", fileName));
        exitFlag = 1;
        return exitFlag;
    }

    auto L = luaL_newstate();
    luaL_openlibs(L);
    registerVector3(L);
    registerPaths(L);
    registerSurfaces(L);
    registerVolumes(L);
    registerUnivariateFunctions(L);
    registerStructuredGrid(L);
    registerUnstructuredGrid(L);
    registerSketch(L);

    if (luaL_dofile(L, toStringz(dirName(thisExePath())~"/foam-mesh.lua")) != 0) {
        writeln("There was a problem in the foamMesh lua code: foam-mesh.lua");
        writeln("------- Error string: --------");
        writeln(to!string(lua_tostring(L, -1)));
        writeln("Exiting foamMesh.");
        exitFlag = 1;
        return exitFlag;
    }

    if (luaL_dofile(L, toStringz(jobName~".lua")) != 0) {
        writeln("There was a problem in the user-supplied job script: ", jobName~".lua");
        writeln("------- Error string: --------");
        writeln(to!string(lua_tostring(L, -1)));
        writeln("Exiting foamMesh.");
        exitFlag = 1;
        return exitFlag;
    }

    // Execute "main()" in foam-mesh.lua
    lua_getglobal(L, "main");
    lua_pushnumber(L, verbosityLevel);
    if (lua_pcall(L, 1, 0, 0) != 0) {
        writeln("There was a problem in the foamMesh main() function in foam-mesh.lua");
        writeln("------- Error string: --------");
        writeln(to!string(lua_tostring(L, -1)));
        writeln("Exiting foamMesh.");
        exitFlag = 1;
        return exitFlag;
    }

    if (verbosityLevel > 0) writeln("foamMesh done.");
    return exitFlag;
}
