-- CO2 no ions CEA backed gas model
-- This is just the CO2 with ions model with the ions and electrons removed.
-- I had to remove C3O2 as well as it was making PITOT a bit unreliable sometimes...
-- Chris James (c.james4@uq.edu.au) - 15/03/22

model = "CEAGas"

CEAGas = {
  mixtureName = 'co2-no-ions',
  speciesList = {'CO2','C2O','CO','O3','O2','O','C3','C2','C'},
  reactants = {CO2=1.0},
  inputUnits = "moles",
  withIons = false,
  trace = 1.0e-6
}
