-- DO NOT ACTUALLY USE THIS GAS MODEL
-- This is just a gas model which PITOT3 can use to set 
-- the fill state (and other very low temperature states)
-- for the related cea gas model (cea-30-h2-15-o2-55-he)
-- Chris James (c.james4@uq.edu.au) - 26/08/23

model = "CEAGas"

CEAGas = {
  mixtureName = 'cea-30-h2-15-o2-55-he',
  speciesList = {"H2", "H", "O2", "O", "He"},
  reactants = {H2=0.3, O2=0.15, He = 0.55},
  inputUnits = "moles",
  withIons = false,
  trace = 1.0e-10
}
