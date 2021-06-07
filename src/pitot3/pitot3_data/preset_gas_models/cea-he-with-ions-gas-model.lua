model = "CEAGas"

CEAGas = {
  mixtureName = 'he-with-ions',
  speciesList = {"He", "He+","e-"},
  reactants = {He=1.0},
  inputUnits = "moles",
  withIons = true,
  trace = 1.0e-6
}
