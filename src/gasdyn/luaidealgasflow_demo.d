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
print("Try out ideal gas flow functions.")
local a = idealgasflow.A_Astar(4.0,1.4)
print("A_Astar(4.0,1.4)= ", a)
print("beta_obl(4.0, 0.1)=", idealgasflow.beta_obl(4.0, 0.1))
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

print("isentropic flow")
local M = 2.4
local g = 1.4
assert(approxEqual(idealgasflow.T0_T(M,g), 2.152), "Total temperature fail")
assert(approxEqual(idealgasflow.p0_p(M,g), 14.620), "Total pressure fail")
assert(approxEqual(idealgasflow.r0_r(M,g), 6.7937), "Total density fail")

print("normal shock")
M = 2.0
g = 1.4
assert(approxEqual(idealgasflow.m2_shock(M,g), 0.5774),
       "Mach number after shock fail")
assert(approxEqual(idealgasflow.p2_p1(M,g), 4.50),
       "Pressure ratio across shock fail")
assert(approxEqual(idealgasflow.T2_T1(M,g), 1.687),
       "Temperature ratio across shock fail")
assert(approxEqual(idealgasflow.r2_r1(M,g), 2.667),
       "Density ratio across shock fail")

print("Rayleigh-line flow")
assert(approxEqual(idealgasflow.T0_T0star(M,g), 0.7934),
       "Rayleigh-line total T0_T0star fail")
assert(approxEqual(idealgasflow.T_Tstar(M,g), 0.5289),
       "Rayleigh-line static T_Tstar fail")
assert(approxEqual(idealgasflow.p_pstar(M,g), 0.3636),
       "Rayleigh-line static p_pstar fail")
assert(approxEqual(idealgasflow.r_rstar(M,g), 0.6875),
       "Rayleigh-line static p_pstar fail")
assert(approxEqual(idealgasflow.M_Rayleigh(idealgasflow.T0_T0star(M,g),g), M),
       "Rayleigh-line inverse fail")

print("isentropic flow turning")
M = 2.4
g = 1.4
assert(approxEqual(idealgasflow.PM1(M,g), 0.6413), "Prandtl-Meyer fail")
local nu = 0.6413479572
assert(approxEqual(idealgasflow.PM2(nu,g), 2.4), "Inverse Prandtl-Meyer fail");
assert(approxEqual(idealgasflow.MachAngle(M), math.asin(1/2.4)),
       "Mach angle fail");

print("oblique shock relations")
M = 2.0
g = 1.4
local beta = math.rad(44.0)
local theta = math.rad(14.0)
assert(approxEqual(idealgasflow.beta_obl(M, theta, g), beta),
       "Oblique shock, beta from theta fail")
assert(approxEqual(idealgasflow.beta_obl2(M, 2.088, g), beta),
       "Oblique shock, beta from p2_p1 fail")
assert(approxEqual(idealgasflow.theta_obl(M, beta, g), theta),
       "Oblique shock, theta from beta fail")
assert(approxEqual(idealgasflow.M2_obl(M, beta, theta, g), 1.482),
       "Oblique shock, M2 after shock fail")
assert(approxEqual(idealgasflow.T2_T1_obl(M, beta, g), 1.249),
       "Oblique shock, temperature ratio fail")
assert(approxEqual(idealgasflow.p2_p1_obl(M, beta, g), 2.088),
       "Oblique shock, pressure ratio fail")
assert(approxEqual(idealgasflow.r2_r1_obl(M, beta, g), 1.673),
       "Oblique shock, density ratio fail")
assert(approxEqual(idealgasflow.p02_p01_obl(M, beta, g), 0.9608),
       "Oblique shock, total-pressure fail")
assert(approxEqual(idealgasflow.u2_u1_obl(M, beta, g), 0.598),
       "Oblique shock, normal velocity ratio fail")
assert(approxEqual(idealgasflow.V2_V1_obl(M, beta,theta, g),0.828),
       "Oblique shock, absolute velocity ratio fail")

print("Done.")
    `;
    if ( luaL_dostring(L, toStringz(test_code)) != 0 ) {
	writeln("There was a problem interpreting the test code.");
	writeln(to!string(lua_tostring(L, -1)));
    }
    writeln("End luaidealgasflow_demo.");
}

