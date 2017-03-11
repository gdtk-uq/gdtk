model = "CEAGas"

CEAGas = {
  mixtureName = 'air7species',
  speciesList = {"N2","O2","N","O","NO","NO+","e-"},
  reactants = {N2=0.79, O2=0.21},
  inputUnits = "moles",
  withIons = true,
  trace = 1.0e-6
}
