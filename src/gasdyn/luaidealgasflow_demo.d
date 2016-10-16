/**
 * luaidealgasflow_demo.d
 * Shows the wrapped idealgasflow functions in action.
 *
 * Author: Peter J. and Rowan G.
 * Date: 2016-10-16, just enough to get the Billig correlation working
 */

import std.stdio;
import std.conv;
import std.string;
import util.lua;
import luaidealgasflow;

void main()
{
    writeln("Begin demo of idealgasflow module for use in Lua.");
    auto L = luaL_newstate();
    luaL_openlibs(L);
    registeridealgasflowFunctions(L);
    string test_code = `
-- Try out some functions.
print("ideal gas flow functions.")
a = idealgasflow.A_Astar(4.0,1.4)
print("A_Astar(4.0,1.4)= ", a)
print("beta_obl(4.0, 0.1)=", idealgasflow.beta_obl(4.0, 0.1))
    `;
    if ( luaL_dostring(L, toStringz(test_code)) != 0 ) {
	writeln("There was a problem interpreting the test code.");
	writeln(to!string(lua_tostring(L, -1)));
    }
    writeln("End luaidealgasflow_demo.");
}

