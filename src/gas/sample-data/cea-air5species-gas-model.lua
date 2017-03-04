model = "CEAGas"

CEAGas = {
  mixtureName = 'air5species',
  reactants = {N2=0.79, O2=0.21},
  speciesList = {"N2", "O2", "N", "O", "NO"},
  inputUnits = "moles",
  outputUnits = "moles",
  withIons = false,
  trace = 1.0e-6
}
