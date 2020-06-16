model = "CEAGas"

CEAGas = {
  mixtureName = 'co2-ions',
  speciesList = {'CO2','C2','C','CO','O2','O','C+','CO+','O2+','O+','e-'},
  reactants = {CO2=1.0},
  inputUnits = "moles",
  withIons = true,
  trace = 1.0e-6
}
