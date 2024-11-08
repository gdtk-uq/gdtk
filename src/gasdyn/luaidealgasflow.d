/**
 * A Lua interface for the D idealgasflow module.
 *
 * Authors: Peter J. and Rowan G.
 * Date: 2016-10-16, just enough for the Billig shock shape correlation
 */

module gasdyn.luaidealgasflow;

import std.stdio;
import std.string;
import std.conv;
import std.algorithm;
import util.lua;
import util.lua_service;
import gasdyn.idealgasflow;

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

string wrapfn(string fname, string[] args, int nreturn=1)
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
    // Now call the function and stack up the result(s).
    code ~= "    try {\n";
    if (nreturn == 1) {
        code ~= "        double result = "~fname~"(";
        foreach(a; anames) { code ~= a~", "; } // a trailing comma is ok
        code ~= ");\n";
        code ~= "        lua_pushnumber(L, result);\n";
        code ~= "        return 1;\n";
    } else {
        assert(nreturn>1, "oops, expected a positive number of return values.");
        code ~= "        double[] results = "~fname~"(";
        foreach(a; anames) { code ~= a~", "; } // a trailing comma is ok
        code ~= ");\n";
        code ~= "        foreach(res; results) { lua_pushnumber(L, res); }\n";
        code ~= "        return to!int(results.length);\n";
    }
    code ~= "    } catch(Exception e) {\n";
    code ~= "        return luaL_error(L, e.msg.toStringz);\n";
    code ~= "    }\n";
    code ~= "}\n";
    return to!string(code);
} // end wrapfn()

mixin(wrapfn("A_Astar", ["M", "g=1.4"]));
mixin(wrapfn("T0_T", ["M", "g=1.4"]));
mixin(wrapfn("p0_p", ["M", "g=1.4"]));
mixin(wrapfn("r0_r", ["M", "g=1.4"]));

mixin(wrapfn("m2_shock", ["M1", "g=1.4"]));
mixin(wrapfn("r2_r1", ["M1", "g=1.4"]));
mixin(wrapfn("u2_u1", ["M1", "g=1.4"]));
mixin(wrapfn("p2_p1", ["M1", "g=1.4"]));
mixin(wrapfn("T2_T1", ["M1", "g=1.4"]));
mixin(wrapfn("p02_p01", ["M1", "g=1.4"]));
mixin(wrapfn("ds_Cv", ["M1", "g=1.4"]));
mixin(wrapfn("pitot_p", ["p1", "M1", "g=1.4"]));

mixin(wrapfn("T0_T0star", ["M", "g=1.4"]));
mixin(wrapfn("M_Rayleigh", ["T0T0star", "g=1.4"]));
mixin(wrapfn("T_Tstar", ["M", "g=1.4"]));
mixin(wrapfn("p_pstar", ["M", "g=1.4"]));
mixin(wrapfn("r_rstar", ["M", "g=1.4"]));
mixin(wrapfn("p0_p0star", ["M", "g=1.4"]));

mixin(wrapfn("PM1", ["M", "g=1.4"]));
mixin(wrapfn("PM2", ["nu", "g=1.4", "tol=1.0e-6"]));
mixin(wrapfn("MachAngle", ["M"]));

mixin(wrapfn("beta_obl", ["M1", "theta", "g=1.4", "tol=1.0e-6"]));
mixin(wrapfn("beta_obl2", ["M1", "p2_p1", "g=1.4"]));
mixin(wrapfn("theta_obl", ["M1", "beta", "g=1.4"]));
mixin(wrapfn("dtan_theta", ["M1", "beta", "g=1.4"]));
mixin(wrapfn("M2_obl", ["M1", "beta", "theta", "g=1.4"]));
mixin(wrapfn("r2_r1_obl", ["M1", "beta", "g=1.4"]));
mixin(wrapfn("Vn2_Vn1_obl", ["M1", "beta", "g=1.4"]));
mixin(wrapfn("V2_V1_obl", ["M1", "beta", "g=1.4"]));
mixin(wrapfn("p2_p1_obl", ["M1", "beta", "g=1.4"]));
mixin(wrapfn("T2_T1_obl", ["M1", "beta", "g=1.4"]));
mixin(wrapfn("p02_p01_obl", ["M1", "beta", "g=1.4"]));

mixin(wrapfn("theta_cone", ["V1", "p1", "T1", "beta", "R=287.1", "g=1.4"], 4));
mixin(wrapfn("beta_cone", ["V1", "p1", "T1", "theta", "R=287.1", "g=1.4"]));
mixin(wrapfn("beta_cone2", ["M1", "theta", "R=287.1", "g=1.4"]));

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
    mixin(registerfn("T0_T"));
    mixin(registerfn("p0_p"));
    mixin(registerfn("r0_r"));

    mixin(registerfn("m2_shock"));
    mixin(registerfn("r2_r1"));
    mixin(registerfn("u2_u1"));
    mixin(registerfn("p2_p1"));
    mixin(registerfn("T2_T1"));
    mixin(registerfn("p02_p01"));
    mixin(registerfn("ds_Cv"));
    mixin(registerfn("pitot_p"));

    mixin(registerfn("T0_T0star"));
    mixin(registerfn("M_Rayleigh"));
    mixin(registerfn("T_Tstar"));
    mixin(registerfn("p_pstar"));
    mixin(registerfn("r_rstar"));
    mixin(registerfn("p0_p0star"));

    mixin(registerfn("PM1"));
    mixin(registerfn("PM2"));
    mixin(registerfn("MachAngle"));

    mixin(registerfn("beta_obl"));
    mixin(registerfn("beta_obl2"));
    mixin(registerfn("theta_obl"));
    mixin(registerfn("dtan_theta"));
    mixin(registerfn("M2_obl"));
    mixin(registerfn("r2_r1_obl"));
    mixin(registerfn("Vn2_Vn1_obl"));
    mixin(registerfn("V2_V1_obl"));
    mixin(registerfn("p2_p1_obl"));
    mixin(registerfn("T2_T1_obl"));
    mixin(registerfn("p02_p01_obl"));

    mixin(registerfn("theta_cone"));
    mixin(registerfn("beta_cone"));
    mixin(registerfn("beta_cone2"));

    lua_setglobal(L, idealgasflowMT.toStringz);
} // end registeridealgasflowFunctions()
