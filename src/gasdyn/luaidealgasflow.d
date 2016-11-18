/**
 * A Lua interface for the D idealgasflow module.
 *
 * Authors: Peter J. and Rowan G.
 * Date: 2016-10-16, just enough for the Billig shock shape correlation
 */

module luaidealgasflow;

import std.stdio;
import std.string;
import std.conv;
import std.algorithm;
import util.lua;
import util.lua_service;
import idealgasflow;

// Name of metatable
immutable string idealgasflowMT = "idealgasflow";

/+ 
Here is our prototype wrapped function, kept for reference.

extern(C) int idealgasflow_A_Astar(lua_State* L)
{
    if (!lua_isnumber(L, 1)) {
	string errMsg = "Expected a number for M";
	luaL_error(L, errMsg.toStringz);
    }
    double mach = to!double(luaL_checknumber(L, 1));
    double g = 1.4; // default value
    if (lua_isnumber(L, 2)) {
	g = to!double(luaL_checknumber(L, 2));
    }
    lua_pushnumber(L, A_Astar(mach, g));
    return 1;
}
+/

string wrapfn(string fname, string[] args)
// Generate the code to wrap D function for access from the idealgasflow table.
// The arguments are assumed to be of type double and may have a default value.
{
    char[] code;
    code ~= "extern(C) int idealgasflow_"~fname~"(lua_State* L)\n";
    code ~= "{\n";
    // Determine our list of functions and their default values, if any.
    string[] anames; string[] values; bool[] default_flags;
    foreach(a; args) {
	bool with_default = canFind(a, "=");
	if (with_default) {
	    auto items = split(a, "=");
	    anames ~= strip(items[0]);
	    values ~= strip(items[1]);
	} else {
	    anames ~= strip(a);
	    values ~= "0.0";
	}
	default_flags ~= with_default;
    }
    // Generate Lua code to get each argument value.
    foreach(i,a; anames) {
	// Lua argument is i+1
	code ~= "    double "~a~" = "~values[i]~";\n";
	if (default_flags[i]) {
	    code ~= "    if (lua_isnumber(L, "~format("%d", i+1)~")) {\n";
	    code ~= "        "~a~" = to!double(luaL_checknumber(L, "~format("%d", i+1)~"));\n";
	    code ~= "    }\n";
	} else {
	    code ~= "    if (!lua_isnumber(L, "~format("%d", i+1)~")) {\n";
	    code ~= "        string errMsg = \"Expected a number for "~a~"\";\n";
	    code ~= "        luaL_error(L, errMsg.toStringz);\n";
	    code ~= "    }\n";
	    code ~= "    "~a~" = to!double(luaL_checknumber(L, "~format("%d", i+1)~"));\n";
	}
    }
    // Now call the function and pack up the result.
    code ~= "    lua_pushnumber(L, "~fname~"(";
    foreach(a; anames) { code ~= a~", "; } // a trailing comma is ok
    code ~= "));\n";
    code ~= "    return 1;\n";
    code ~= "}\n";
    return to!string(code);
} // end wrapfn()

mixin(wrapfn("A_Astar", ["M", "g=1.4"]));
mixin(wrapfn("beta_obl", ["M1", "theta", "g=1.4", "tol=1.0e-6"]));
mixin(wrapfn("M2_obl", ["M1", "beta", "theta", "g=1.4"]));
mixin(wrapfn("p2_p1_obl", ["M1", "beta", "g=1.4"]));

string registerfn(string fname)
{
    return "    lua_pushcfunction(L, &idealgasflow_"~fname~");\n" ~
	"    lua_setfield(L, -2, \""~fname~"\");";
}

void registeridealgasflowFunctions(lua_State* L)
{
    // Register the idealgasflow table of functions.
    luaL_newmetatable(L, idealgasflowMT.toStringz);
    
    /* metatable.__index = metatable */
    lua_pushvalue(L, -1); // duplicates the current metatable
    lua_setfield(L, -2, "__index");

    mixin(registerfn("A_Astar"));
    mixin(registerfn("beta_obl"));
    mixin(registerfn("M2_obl"));
    mixin(registerfn("p2_p1_obl"));

    lua_setglobal(L, idealgasflowMT.toStringz);
} // end registeridealgasflowFunctions()
