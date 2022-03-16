-- DO NOT ACTUALLY USE THIS GAS MODEL
-- This is just a gas model which PITOT3 can use to set 
-- the fill state (and other very low temperature states)
-- for the related cea gas model (cea-mars-c02-n2-ar-with-ions)
-- Chris James (c.james4@uq.edu.au) - 16/03/22

model = "CEAGas"

CEAGas = {
  mixtureName = 'mars-c02-n2-ar-with-ions-room-temperature-only-gas-model',
  speciesList = {'CO2','N2','Ar'},
  reactants = {CO2=0.958,N2=0.027,Ar=0.015},
  inputUnits = "moles",
  withIons = false,
  trace = 1.0e-6
}
