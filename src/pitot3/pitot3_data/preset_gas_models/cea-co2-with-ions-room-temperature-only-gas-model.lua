-- DO NOT ACTUALLY USE THIS GAS MODEL
-- This is just a gas model which PITOT3 can use to set 
-- the fill state (and other very low temperature states)
-- for the related cea gas model (cea-co2-with-ions)
-- Chris James (c.james4@uq.edu.au) - 15/03/22

model = "CEAGas"

CEAGas = {
  mixtureName = 'co2-with-ions-room-temperature-only-gas-model',
  speciesList = {'CO2'},
  reactants = {CO2=1.0},
  inputUnits = "moles",
  withIons = false,
  trace = 1.0e-6
}
