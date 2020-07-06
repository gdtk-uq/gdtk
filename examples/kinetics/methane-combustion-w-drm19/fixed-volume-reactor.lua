-- Author: Rowan J. Gollan
-- Date: 2017-03-23
--
-- A simple fixed-volume reactor.
--
-- This script is designed to run with the gas-calc program
-- that comes as part of the dgd collection. To run this and
-- capture the output in a data file, do:
--
-- > gas-calc fixed-volume-reactor.lua > output.data
--
-- This script is a conversion of Peter Jacobs'
-- fixed_volume_reactor.py that is part of the cfcfd3
-- code collection. That script is, in turn, a cut-down
-- version of Brendan O'Flaherty's master program for
-- testing many versions of the chemical reactor.
--

gmodel = GasModel:new{'drm19-gas-model.lua'}
Q = gmodel:createGasState()
total = 1+2+7.52
molef = {CH4=1/total, O2=2/total, N2=7.52/total}
Q.T = 2000.0 -- K
Q.p = 101.325e3 -- Pa
Q.massf = gmodel:molef2massf(molef)
gmodel:updateThermoFromPT(Q)
nsp = gmodel:nSpecies()

function writeHeader()
   str = "# 1:t  2:T_0  3:p_0  4:T  5:p  6:rho"
   for isp=0,nsp-1 do
      str = str..", "..gmodel:speciesName(isp)
   end
   print(str)
   return
end

function writeData(t)
   str = string.format("%.5e %.5e %.5e %.5e %.5e %.12e ",
		       t, Q.T, Q.p, Q.T, Q.p, Q.rho)
   molef = gmodel:massf2molef(Q)
   for isp=0,nsp-1 do
      str = str..string.format(" %.12e", molef[gmodel:speciesName(isp)])
   end
   print(str)
   return
end

chem = ChemistryUpdate:new{filename='drm19-reaction-scheme.lua', gasmodel=gmodel}

tEnd = 4.0e-4
dtChem = -1.0
dt = tEnd/2000
t = 0.0

writeHeader()
writeData(t)

while t < tEnd do
   dtChem = chem:updateState(Q, dt, dtChem, gmodel)
   gmodel:updateThermoFromRHOU(Q)
   t = t + dt
   writeData(t)
end

print("# Done.")
