model = "CEAGas"

CEAGas = {
  mixtureName = 'air5species',
  speciesList = {"N2", "O2", "N", "O", "NO"},
  reactants = {N2=0.79, O2=0.21},
  inputUnits = "moles",
  withIons = false,
  trace = 1.0e-6
}
