model = "CEAGas"

CEAGas = {
  mixtureName = 'air13species',
  speciesList = {"N2","O2","Ar","N","O","NO","Ar+","NO+","N+","O+","N2+","O2+","e-"},
  reactants = {N2=0.7811, O2=0.2095, Ar=0.0093},
  inputUnits = "moles",
  withIons = true,
  trace = 1.0e-6
}
