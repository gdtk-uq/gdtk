-- air 23 species CEA backed gas model
-- This is based off the composition when one selects 'Air' in CEA,
-- it is air with both Ar and CO2 (air13species is just air with argon)
-- I got the species list from intuition and some pitot runs, hopefully it is right.
-- Chris James (c.james4@uq.edu.au) - 14/09/21
model = "CEAGas"

CEAGas = {
  mixtureName = 'air23species',
  speciesList = {'N2','N3','O2','NO','N2O','NO2','CO2','CO','CN','C','Ar','N','O','NO+','N2+','O2+','Ar+','C+','N+','N-','O+','O-','e-'},
  reactants = {N2=0.78084, O2=0.20948, Ar=0.009365,CO2=0.000319},
  inputUnits = "moles",
  withIons = true,
  trace = 1.0e-6
}
