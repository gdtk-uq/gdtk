/**
 * luagasflow_demo.d
 * Shows the wrapped gasflow functions in action.
 *
 * Author: Peter J. and Rowan G.
 * Date: 2017-04-02, just enough to get nenzfr2 demo going
 */

import std.stdio;
import std.conv;
import std.string;
import util.lua;
import gas.luagas_model;
import luagasflow;

void main()
{
    writeln("Begin demo of gasflow module for use in Lua.");
    auto L = luaL_newstate();
    luaL_openlibs(L);
    registerGasModel(L, LUA_GLOBALSINDEX);
    registergasflowFunctions(L);
    string test_code = `
print("Try out gasflow functions.")
-- Initialise an ideal gas model
gm = GasModel:new{'../gas/sample-data/cea-air13species-gas-model.lua'}
state1 = GasState:new{gm}
state1.p = 125.0e3; state1.T = 300.0
gm:updateThermoFromPT(state1)
print("state1:"); printValues(state1)
state2 = GasState:new{gm}

print("Test some assertions")
function approxEqual(a, b, relTol, absTol)
   relTol = relTol or 1.0e-2
   absTol = absTol or 1.0e-5
   local diff = math.abs(a - b);
   local result = false
   if diff < absTol then
      result = true
   else
      mag = math.max(math.abs(a), math.abs(b))
      result = (diff/mag < relTol)
   end
   -- print("a=", a, "b=", b, "relTol=", relTol, "absTol=", absTol, "diff=", diff)
   return result
end

print("normal shock")
Vs = 2414.0
V2, Vg = gasflow.normal_shock(state1, Vs, state2, gm)
assert(approxEqual(V2, 361.9), "V2 number after shock fail")
assert(approxEqual(Vg, 2052.1), "Vg number after shock fail")

print("Done.")
    `;
    if ( luaL_dostring(L, toStringz(test_code)) != 0 ) {
	writeln("There was a problem interpreting the test code.");
	writeln(to!string(lua_tostring(L, -1)));
    }
    writeln("End luagasflow_demo.");
}
