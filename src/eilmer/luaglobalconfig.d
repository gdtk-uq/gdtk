// luaglobalconfig.d
// Lua access to the GlobalConfig class data, for use in the preparation script.
//
// Peter J. and Rowan G.
// 2015-03-02: First code adapted from the other lua wrapper modules.

module luaglobalconfig;

import core.stdc.stdlib : exit;

import std.stdio;
import std.string;
import std.conv;
import util.lua;
import util.lua_service;

import gas;
import gas.luagas_model;
import globalconfig;


//------------------------------------------------------------------------
