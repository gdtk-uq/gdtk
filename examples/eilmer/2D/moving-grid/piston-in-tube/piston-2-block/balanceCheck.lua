-- Author: Fabian Zander
-- Data: 2019-05-03
--
-- This script gives us some checking of the moving mesh processes.
-- We look at the total gas mass and energy at the start and end of the simulation.
--
-- NOTES:
-- We make use of some global data and functions that Eilmer
-- provides use for doing the work. In particular, we use:
-- 1. FlowSolution:
--    Here we can pick up the computed flow solutions so that we can do some
--    data manipulation.
-- 2. GasState:
--    We are using the internally available gas models to compute the internal
--    of the gas.
-- We can not pick up the '0' time step for the initial energy so
-- we are using the first time step, this requires the piston velocity
-- energy to be taken into account at the start as well as the end.

-- Our simulation has two blocks.
-- We will pick up the first,non-initial and the last flow solution.
local nb = 2
local fsolStart = FlowSolution:new{jobName="pit2", dir=".", tindx=1, nBlocks=nb}
local fsolEnd = FlowSolution:new{jobName="pit2", dir=".", tindx=8, nBlocks=nb}
print("fsolStart=", fsolStart)
print("fsolEnd=", fsolEnd)
--
-- Initialise our running sums.
local massStart = 0.0
local energyStart = 0.0
local massEnd = 0.0
local energyEnd = 0.0
-- We need to specify a gas model for our calculations,
-- so that we can compute internal energy for the gas in each cell.
local gmodel = GasModel:new{"ideal-air-gas-model.lua"}
-- Work through the blocks and, for each cell within each block,
-- accumulate the mass and internal-plus-kinetic energy.
for blk = 0, nb-1 do
   for i = 0, fsolStart:get_nic(0)-1 do
      for j = 0, fsolStart:get_njc(0)-1 do
         -- Starting with the first step
         cell0 = fsolStart:get_cell_data{ib=blk,i=i,j=j}
         cell0mass = cell0.volume * 2 * math.pi * cell0.rho
         -- Configure the gas state using the code internal models
         Q0 = GasState:new{gmodel}
         Q0.p = cell0.p; Q0.T = cell0.T
         gmodel:updateThermoFromPT(Q0)
         -- Sum up the internal and kinetic energy of the cell and accumulate.
         energyStart = energyStart + cell0mass * 
            (Q0.u + 0.5 * (cell0['vel.x']^2 + cell0['vel.y']^2))
         -- Calculate the mass
         massStart = massStart + cell0.volume * 2 * math.pi * cell0.rho
         -- Now the last time step
         cell1 = fsolEnd:get_cell_data{ib=blk,i=i,j=j}
         cell1mass = cell1.volume * 2 * math.pi * cell1.rho
         -- Trying to calculate the total internal energy
         Q1 = GasState:new{gmodel}
         Q1.p = cell1.p; Q1.T = cell1.T
         gmodel:updateThermoFromPT(Q1)
         energyEnd = energyEnd + cell1mass *
            (Q1.u  + 0.5 * (cell1['vel.x']^2 + cell1['vel.y']^2))
         -- Calculate the mass
         massEnd = massEnd + cell1.volume * 2 * math.pi * cell1.rho
      end
   end
end
-- These values have to be entered manually at this stage from piston.data
-- I need the 'start' velocity as I can't use flow solution 0, i.e. this is
-- the first time step that has been written, tindx=1
local pStartVel = 5.647227200338181774e+01
local pEndVel = 2.774190962612635190e+02
local pStartEnergy = 0.5*1.0*pStartVel^2
local pEndEnergy = 0.5*1.0*pEndVel^2
-- Print out some results including a calculation of the energy conservation
print('')
print('Gas-mass-early = ', massStart)
print('Gas-mass-at-end = ', massEnd)
print('Gas-energy-early = ', energyStart)
print('Gas-energy-at-end = ', energyEnd)
print('Piston-Energy-at-end = ', pEndEnergy)
print('')
local ee = energyEnd + pEndEnergy - energyStart - pStartEnergy
local ee_rel = ee / (energyStart+pStartEnergy)
print(string.format('Energy-error = %.2f (%.2f%%)', ee, ee_rel*100))
print('')
