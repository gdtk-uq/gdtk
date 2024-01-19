-- updated the trace to 1.0e-10 as running the CEA calculations in massf (which is the GDTk's only way to do it)
-- suppresses ionisation otherwise as electrons are very light so their mass fraction is very small
-- Chris James (c.james4@uq.edu.au) - 19/01/24

model = "CEAGas"

CEAGas = {
  mixtureName = 'he-with-ions',
  speciesList = {"He", "He+","e-"},
  reactants = {He=1.0},
  inputUnits = "moles",
  withIons = true,
  trace = 1.0e-10
}
