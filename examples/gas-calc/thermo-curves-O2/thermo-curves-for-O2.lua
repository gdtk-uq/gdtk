-- A script to output Cp and h for O2 over temperature range 200--20000 K.
--
-- Author: Rowan J. Gollan
-- Date: 2016-12-23
--
-- To run this script:
-- $ prep-gas O2.inp O2-gas-model.lua
-- $ gas-calc thermo-curves-for-O2.lua
--

gasModelFile = 'O2-gas-model.lua'
gmodel = GasModel:new{gasModelFile}

Q = GasState:new{gmodel}
Q.p = 1.0e5 -- Pa
Q.massf = {O2=1.0}

outputFile = 'O2-thermo.dat'
print("Opening file for writing: ", outputFile)
f = assert(io.open(outputFile, "w"))
f:write("#  1:T[K]      2:Cp[J/kg/K]     3:h[J/kg]\n")

Tlow = 200.0
Thigh = 20000.0
dT = 100.0

for T=Tlow,Thigh,dT do
   Q.T = T
   gmodel:updateThermoFromPT(Q)
   Cp = gmodel:Cp(Q)
   h = gmodel:enthalpy(Q)
   f:write(string.format(" %12.6e %12.6e %12.6e\n", T, Cp, h))
end

f:close()
print("File closed. Done.")

