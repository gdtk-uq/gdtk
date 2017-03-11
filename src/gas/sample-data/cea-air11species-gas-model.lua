model = "CEAGas"

CEAGas = {
  mixtureName = 'air11species',
  speciesList = {"N2","O2","N","O","NO","NO+","N+","O+","N2+","O2+","e-"},
  reactants = {N2=0.79, O2=0.21},
  inputUnits = "moles",
  withIons = true,
  trace = 1.0e-6
}
