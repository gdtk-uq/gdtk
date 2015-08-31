/**
 * luaflowstate_demo.d
 * Demonstrate the wrapped FlowState object.
 *
 * Authors: Rowan G. and Peter J.
 * Version: 2015-02-23
 */

import std.stdio;

import util.lua;
import luageom;
import luaflowstate;
import luaglobalconfig;

void main()
{
    writeln("Begin demonstration of Lua connection to FlowState object.");
    auto L = luaL_newstate();
    luaL_openlibs(L);
    registerVector3(L);
    registerGlobalConfig(L);
    registerFlowState(L);
    luaL_dostring(L, `
nsp, nmodes = setGasModel('sample-data/ideal-air-gas-model.lua')
print("GasModel set to ideal air. nsp= ", nsp, " nmodes= ", nmodes)
fs = FlowState:new{p=1.0e5, T=300.0, velx=1000.0, vely=200.0}
fsTab = fs:toTable{}
for k,v in pairs(fsTab) do
    print(k,v)
    if ( k == 'gas' ) then
       for k1,v1 in pairs(v) do
           print(k1, v1)
           if ( k1 == 'massf' ) then
               print("massf[1]= ", v1[1])
           end
       end
    end
end
-- Try to set tke and omega
fs:fromTable{tke=30.0, omega=500.0, velz=3000.0}
fsTab = fs:toTable{}
print("tke= ", fsTab.tke, " omega= ", fsTab.omega, "velz= ", fsTab.velz)

    `);
    writeln("Done with demo.");
}
