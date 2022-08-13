/**
 * luagas_model_demo.d
 *
 * Author: Rowan G.
 * Version: Initial cut.
 */

import std.stdio;
import std.conv;
import std.string;

import util.lua;
import gas.luagas_model;

void main()
{
    writeln("Begin demonstration of lua connection to the gas module.");
    auto L = luaL_newstate();
    luaL_openlibs(L);
    lua_pushglobaltable(L);
    registerGasModel(L);

    string test_code = `
-- Initialise an ideal gas model
gm = GasModel:new{'./sample-data/ideal-air-gas-model.lua'}
-- Try out some of the service functions
print("number of species= ", gm:nSpecies())
print("number of modes= ", gm:nModes())
print("speciesName(1)= ", gm:speciesName(0))
print("list of mol masses:")
for i, m in ipairs(gm:molMasses()) do
   print("i, mMass= ", i, m)
end
print("Test thermo evaluations....")
Q = GasState:new{gm}
print("for empty gas state Q.p= ", Q.p)
Q.p = 1.0e5; Q.T = 300.0
print("update based on p-T")
gm:updateThermoFromPT(Q)
print("Q.u= ", Q.u)
print("update based on rho-u")
gm:updateThermoFromRHOU(Q)
print("Q.p= ", Q.p, " Q.T= ", Q.T)
printValues(Q)
    `;
    if ( luaL_dostring(L, toStringz(test_code)) != 0 ) {
        writeln("There was a problem interpreting the test code.");
        writeln(to!string(lua_tostring(L, -1)));
    }
    writeln("Done with luagas_model_demo.");
}
