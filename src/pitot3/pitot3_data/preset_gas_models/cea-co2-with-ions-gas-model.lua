-- CO2 with ions CEA backed gas model
-- I got the species list from old PITOT code runs with CO2 at various conditions.
-- This was the list below:
-- {'CO2','C2O','CO','O3','O2','O','C2','C','C2+','CO+','C+','C-','O+','O-','e-'}
-- I then went into CEA to check the full set of species and added any other species which did not cause
-- PITOT3 to fail when setting up the test gas state at room temperature. I added C302, C3, and O2+.
-- Chris James (c.james4@uq.edu.au) - 12/06/21

model = "CEAGas"

CEAGas = {
  mixtureName = 'co2-with-ions',
  speciesList = {'CO2','C3O2','C2O','CO','O3','O2','O','C3','C2','C','O2+','C2+','CO+','C+','C-','O+','O-','e-'},
  reactants = {CO2=1.0},
  inputUnits = "moles",
  withIons = true,
  trace = 1.0e-6
}
