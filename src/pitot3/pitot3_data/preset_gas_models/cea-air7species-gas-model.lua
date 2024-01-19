-- updated the trace to 1.0e-20 as running the calculation in massf (which is the GDTk default)
-- suppresses ionisation otherwise as ions are very light on a mass basis.
-- Chris James (c.james4@uq.edu.au) - 19/01/24

model = "CEAGas"

CEAGas = {
  mixtureName = 'air7species',
  speciesList = {"N2","O2","N","O","NO","NO+","e-"},
  reactants = {N2=0.79, O2=0.21},
  inputUnits = "moles",
  withIons = true,
  trace = 1.0e-20
}
