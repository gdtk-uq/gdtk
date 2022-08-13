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
import luaidealgasflow;
import luagasflow;

void main()
{
    writeln("Begin demo of gasflow module for use in Lua.");
    auto L = luaL_newstate();
    luaL_openlibs(L);
    lua_pushglobaltable(L);
    registerGasModel(L);
    registeridealgasflowFunctions(L);
    registergasflowFunctions(L);
    string test_code = `
print("Try out gasflow functions.")
-- Initialise a gas model for reacting air in chemical equilibrium.
gm = GasModel:new{'../gas/sample-data/cea-air13species-gas-model.lua'}
state1 = GasState:new{gm}
state1.p = 125.0e3; state1.T = 300.0
gm:updateThermoFromPT(state1)
print("state1:"); printValues(state1)

print("Test some assertions")
function isClose(a, b, relTol, absTol)
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

print("normal shock, given shock speed (rho, T iterations)")
Vs = 2414.0
state2, V2, Vg = gasflow.normal_shock(state1, Vs)
print("    V2=", V2, "Vg=", Vg)
print("    state2:"); printValues(state2)
assert(isClose(V2, 361.9), "V2 number after shock fail")
assert(isClose(Vg, 2052.1), "Vg number after shock fail")
assert(isClose(state2.p, 7.314e6), "p2 number after shock fail")
assert(isClose(state2.T, 2630.0), "T2 number after shock fail")

print("normal shock, given shock speed (p, T iterations)")
Vs = 2414.0
state2, V2, Vg = gasflow.normal_shock_1(state1, Vs)
print("    V2=", V2, "Vg=", Vg)
print("    state2:"); printValues(state2)
assert(isClose(V2, 361.9), "V2 number after shock fail")
assert(isClose(Vg, 2052.1), "Vg number after shock fail")
assert(isClose(state2.p, 7.314e6), "p2 number after shock fail")
assert(isClose(state2.T, 2630.0), "T2 number after shock fail")

print("normal shock computed from pressure ratio")
V1, V2, Vg = gasflow.normal_shock_p2p1(state1, 7.314e6/125.0e3)
print("    V1=", V1, " V2=", V2, " Vg=", Vg)
assert(isClose(V1, 2414.0), "V1 number after p2p1 shock fail")
assert(isClose(V2, 361.9), "V2 number after p2p1 shock fail")
assert(isClose(Vg, 2052.1), "Vg number after p2p1 shock fail")

print("reflected shock")
state5, Vr = gasflow.reflected_shock(state2, Vg)
print("    Vr=", Vr)
print("    state5:"); printValues(state5)
assert(isClose(Vr, 573.9), "Vr number after reflected shock fail")
assert(isClose(state5.p, 59.47e6), "p5 number after reflected shock fail")
assert(isClose(state5.T, 4551.8), "T5 number after reflected shock fail")

print("Expand from stagnation (with ratio of pressure to match observation)")
state5s, V5s = gasflow.expand_from_stagnation(state5, 34.37/59.47)
print("    V5s=", V5s, " Mach=", V5s/state5s.a)
print("    state5s:"); printValues(state5s)
print("    (h5s-h1)=", gm:enthalpy(state5s) - gm:enthalpy(state1))
assert(isClose(V5s, 1184.7), "V5s number after expand_from_stagnation fail")
assert(isClose(state5s.p, 34.37e6), "p5s number after expand_from_stagnation fail")
assert(isClose(state5s.T, 4161.8), "T5s number after expand_from_stagnation fail")

print("Expand to throat condition (Mach 1.0001)")
state6, V6 = gasflow.expand_to_mach(state5s, 1.0001)
print("    V6=", V6, " Mach=", V6/state6.a)
print("    state6:"); printValues(state6)
assert(isClose(V6, 1155.8), "V6 number after expand_to_mach fail")
assert(isClose(V6/state6.a, 1.0001), "mach number after expand_to_mach fail")
assert(isClose(state6.p, 19.32e6), "p6 number after expand_to_mach fail")
assert(isClose(state6.T, 3788.0), "T6 number after expand_to_mach fail")

print("Something like Mach 4 nozzle")
state7, V7 = gasflow.steady_flow_with_area_change(state6, V6, 27.0)
print("    V7=", V7, " Mach=", V7/state7.a)
print("    state7:"); printValues(state7)
assert(isClose(V7, 2950.7), "V7 number after steady_flow_with_area_change fail")
assert(isClose(V7/state7.a, 4.238), "mach number after steady_flow_with_area_change fail")
assert(isClose(state7.p, 93.598e3), "p7 number after steady_flow_with_area_change fail")
assert(isClose(state7.T, 1282.5), "T7 number after steady_flow_with_area_change fail")

print("Total condition")
state8 = gasflow.total_condition(state7, V7)
print("    state8:"); printValues(state8)

print("Pitot condition")
state9 = gasflow.pitot_condition(state7, V7)
print("    state9:"); printValues(state9)
assert(isClose(state9.p/state8.p, 0.06253), "pitot/total pressure ratio fail")

print("\nFinite wave process along a cplus characteristic, stepping in pressure.")
V1 = 0.0; state1.p = 1.0e5; state1.T = 320.0 -- ideal air, not high T
gm:updateThermoFromPT(state1)
gm:updateSoundSpeed(state1)
Jplus = V1 + 2*state1.a/(1.4-1)
state2, V2 = gasflow.finite_wave_dp(state1, V1, "cplus", 60.0e3, 500)
print("    V2=", V2)
print("    state2:"); printValues(state2)
print("    ideal V2=", Jplus - 2*state2.a/(1.4-1))
assert(isClose(V2, 126.2), "velocity after finite_wave_dp fail")
assert(isClose(state2.p, 60.0e3), "pressure after finite_wave_dp fail")
assert(isClose(state2.T, 276.5), "temperature after finite_wave_dp fail")

print("\nFinite wave process along a cplus characteristic, stepping in velocity.")
V1 = 0.0; state1.p = 1.0e5; state1.T = 320.0 -- ideal air, not high T
gm:updateThermoFromPT(state1)
gm:updateSoundSpeed(state1)
Jplus = V1 + 2*state1.a/(1.4-1)
state2, V2 = gasflow.finite_wave_dv(state1, V1, "cplus", 125.0)
print("    V2=", V2)
print("    state2:"); printValues(state2)
print("    ideal V2=", Jplus - 2*state2.a/(1.4-1))
assert(isClose(V2, 125.0), "velocity after finite_wave_dv fail")
assert(isClose(state2.p, 60.3e3), "pressure after finite_wave_dv fail")
assert(isClose(state2.T, 276.9), "temperature after finite_wave_dv fail")

M1 = 1.5
print("\nOblique-shock demo for M1=", M1)
state1.p = 1.0e5; state1.T = 300.0 -- ideal air, not high T
gm:updateThermoFromPT(state1)
gm:updateSoundSpeed(state1)
beta = 45.0 * math.pi/180.0
print("    given beta(degrees)=", beta*180/math.pi)
V1 = 1.5 * state1.a
print("    state1:"); printValues(state1)
state2, theta, V2 = gasflow.theta_oblique(state1, V1, beta)
print("    theta=", theta)
print("    V2=", V2)
print("    state2:"); printValues(state2)
print("    c.f. ideal gas angle=", idealgasflow.theta_obl(M1, beta))

print("Oblique shock angle from deflection.")
beta2 = gasflow.beta_oblique(state1, V1, theta)
print("    beta2(degrees)=", beta2*180/math.pi)
assert(isClose(beta, beta2), "shock wave angle fail")

print("\nCatch an error")
-- Note that we must wrap the code (that might throw errors
-- that we want to catch and handle in Lua)
-- in a function before calling it with pcall.
-- bad_fun = function() error("deliberate oops") end
bad_fun = function() myBeta = gasflow.beta_oblique(state1, 0.5*V1, theta) end
local ok, msg = pcall(bad_fun)
if ok then
   print("Did not correctly catch the subsonic Mach number.")
else
   print("Correctly caught subsonic Mach number, ", msg)
end

print("Done.")
    `;
    if ( luaL_dostring(L, toStringz(test_code)) != 0 ) {
        writeln("There was a problem interpreting the test code.");
        writeln(to!string(lua_tostring(L, -1)));
    }
    writeln("End luagasflow_demo.");
}
