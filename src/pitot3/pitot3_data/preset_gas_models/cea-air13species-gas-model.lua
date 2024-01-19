-- updated the trace to 1.0e-10 as running the CEA calculations in massf (which is the GDTk's only way to do it)
-- suppresses ionisation otherwise as electrons are very light so their mass fraction is very small
-- Chris James (c.james4@uq.edu.au) - 19/01/24

model = "CEAGas"

CEAGas = {
  mixtureName = 'air13species',
  speciesList = {"N2","O2","Ar","N","O","NO","Ar+","NO+","N+","O+","N2+","O2+","e-"},
  reactants = {N2=0.7811, O2=0.2095, Ar=0.0093},
  inputUnits = "moles",
  withIons = true,
  trace = 1.0e-10
}
