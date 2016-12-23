-- A script to compute the viscosity and thermal conductivity
-- of air (as a mixture of N2 and O2) from 200 -- 20000 K.
--
-- Author: Rowan J. Gollan
-- Date: 2016-12-23
--
-- To run this calculation:
-- $ prep-gas thermally-perfect-N2-O2.inp thermally-perfect-N2-O2.lua
-- $ gas-calc transport-properties-for-air.lua
--

gasModelFile = 'thermally-perfect-N2-O2.lua'
gmodel = GasModel:new{gasModelFile}

Q = GasState:new{gmodel}
Q.p = 1.0e5 -- Pa
Q.massf = {N2=0.78, O2=0.22} -- a good approximation for the composition of air

outputFile = 'trans-props-air.dat'
print("Opening file for writing: ", outputFile)
f = assert(io.open(outputFile, "w"))
f:write("#  1:T[K]      2:mu[Pa.s]      3:k[W/(m.K)]\n")

Tlow = 200.0
Thigh = 20000.0
dT = 100.0

for T=Tlow,Thigh,dT do
   Q.T = T
   gmodel:updateThermoFromPT(Q)
   gmodel:updateTransCoeffs(Q)
   f:write(string.format(" %12.6e %12.6e %12.6e\n", T, Q.mu, Q.k))
end

f:close()
print("File closed. Done.")
