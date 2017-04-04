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

print("normal shock, given shock speed")
Vs = 2414.0
state2, V2, Vg = gasflow.normal_shock(state1, Vs)
print("    V2=", V2, "Vg=", Vg)
print("    state2:"); printValues(state2)
assert(approxEqual(V2, 361.9), "V2 number after shock fail")
assert(approxEqual(Vg, 2052.1), "Vg number after shock fail")
assert(approxEqual(state2.p, 7.314e6), "p2 number after shock fail")
assert(approxEqual(state2.T, 2630.0), "T2 number after shock fail")

print("normal shock computed from pressure ratio")
V1, V2, Vg = gasflow.normal_shock_p2p1(state1, 7.314e6/125.0e3)
print("    V1=", V1, " V2=", V2, " Vg=", Vg)
assert(approxEqual(V1, 2414.0), "V1 number after p2p1 shock fail")
assert(approxEqual(V2, 361.9), "V2 number after p2p1 shock fail")
assert(approxEqual(Vg, 2052.1), "Vg number after p2p1 shock fail")

print("reflected shock")
state5, Vr = gasflow.reflected_shock(state2, Vg)
print("    Vr=", Vr)
print("    state5:"); printValues(state5)
assert(approxEqual(Vr, 573.9), "Vr number after reflected shock fail")
assert(approxEqual(state5.p, 59.47e6), "p5 number after reflected shock fail")
assert(approxEqual(state5.T, 4551.8), "T5 number after reflected shock fail")

print("Expand from stagnation (with ratio of pressure to match observation)")
state5s, V5s = gasflow.expand_from_stagnation(state5, 34.37/59.47)
print("    V5s=", V5s, " Mach=", V5s/state5s.a)
print("    state5s:"); printValues(state5s)
print("    (h5s-h1)=", gm:enthalpy(state5s) - gm:enthalpy(state1)) 
assert(approxEqual(V5s, 1184.7), "V5s number after expand_from_stagnation fail")
assert(approxEqual(state5s.p, 34.37e6), "p5s number after expand_from_stagnation fail")
assert(approxEqual(state5s.T, 4161.8), "T5s number after expand_from_stagnation fail")

print("Expand to throat condition (Mach 1.0001)")
state6, V6 = gasflow.expand_to_mach(state5s, 1.0001)
print("    V6=", V6, " Mach=", V6/state6.a)
print("    state6:"); printValues(state6)
assert(approxEqual(V6, 1155.8), "V6 number after expand_to_mach fail")
assert(approxEqual(V6/state6.a, 1.0001), "mach number after expand_to_mach fail")
assert(approxEqual(state6.p, 19.32e6), "p6 number after expand_to_mach fail")
assert(approxEqual(state6.T, 3788.0), "T6 number after expand_to_mach fail")

print("Something like Mach 4 nozzle")
state7, V7 = gasflow.steady_flow_with_area_change(state6, V6, 27.0)
print("    V7=", V7, " Mach=", V7/state7.a)
print("    state7:"); printValues(state7)
assert(approxEqual(V7, 2950.7), "V7 number after steady_flow_with_area_change fail")
assert(approxEqual(V7/state7.a, 4.238), "mach number after steady_flow_with_area_change fail")
assert(approxEqual(state7.p, 93.598e3), "p7 number after steady_flow_with_area_change fail")
assert(approxEqual(state7.T, 1282.5), "T7 number after steady_flow_with_area_change fail")

print("Total condition")
state8 = gasflow.total_condition(state7, V7)
print("    state8:"); printValues(state8)

print("Pitot condition")
state9 = gasflow.pitot_condition(state7, V7)
print("    state9:"); printValues(state9)
assert(approxEqual(state9.p/state8.p, 0.06253), "pitot/total pressure ratio fail")

print("Done.")
    `;
    if ( luaL_dostring(L, toStringz(test_code)) != 0 ) {
	writeln("There was a problem interpreting the test code.");
	writeln(to!string(lua_tostring(L, -1)));
    }
    writeln("End luagasflow_demo.");
}
