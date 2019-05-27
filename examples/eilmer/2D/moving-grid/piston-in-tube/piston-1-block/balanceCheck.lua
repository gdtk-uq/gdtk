-- Author: Fabian Zander
-- Date: 2019-05-03; trimmed to 1 block 2019-05-21
--
-- This script gives us some checking of the moving mesh processes.
-- We will pick up the first,non-initial and the last flow solution.
local fsolStart = FlowSolution:new{jobName="pit1", dir=".", tindx=1, nBlocks=1}
local fsolEnd = FlowSolution:new{jobName="pit1", dir=".", tindx=8, nBlocks=1}
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
for i = 0, fsolStart:get_nic(0)-1 do
   for j = 0, fsolStart:get_njc(0)-1 do
      -- Starting with the first step
      cell0 = fsolStart:get_cell_data{ib=0,i=i,j=j}
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
      cell1 = fsolEnd:get_cell_data{ib=0,i=i,j=j}
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
-- These values have to be entered manually at this stage from piston.data
-- I need the 'start' velocity as I can't use flow solution 0, i.e. this is
-- the first time step that has been written, tindx=1
local pStartVel = 5.631783889448187352e+01
local pEndVel = 2.774202827533125628e+02
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
