-- Ar with ions CEA backed gas model
-- I got the species list from intuition and various old PITOT runs with this composition.
-- Chris James (c.james4@uq.edu.au) - 12/06/21
-- updated the trace to 1.0e-10 as running the CEA calculations in massf (which is the GDTk's only way to do it)
-- suppresses ionisation otherwise as electrons are very light so their mass fraction is very small
-- Chris James (c.james4@uq.edu.au) - 19/01/24

model = "CEAGas"

CEAGas = {
  mixtureName = 'ar-with-ions',
  speciesList = {'Ar','Ar+','e-'},
  reactants = {Ar=1.0},
  inputUnits = "moles",
  withIons = true,
  trace = 1.0e-10
}
