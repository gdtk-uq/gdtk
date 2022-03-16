-- DO NOT ACTUALLY USE THIS GAS MODEL
-- This is just a gas model which PITOT3 can use to set 
-- the fill state (and other very low temperature states)
-- for the related cea gas model (cea-mars-c02-n2-with-ions)
-- Chris James (c.james4@uq.edu.au) - 16/03/22

model = "CEAGas"

CEAGas = {
  mixtureName = 'mars-c02-n2-with-ions-room-temperature-only-gas-model',
  speciesList = {'CO2','N2'},
  reactants = {CO2=0.96,N2=0.04},
  inputUnits = "moles",
  withIons = false,
  trace = 1.0e-6
}
