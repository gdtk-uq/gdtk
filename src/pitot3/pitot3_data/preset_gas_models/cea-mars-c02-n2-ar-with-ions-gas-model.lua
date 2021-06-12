-- Mars (C02,N2,Ar) with ions CEA backed gas model
-- The old PITOT code had Mars gas defined as just CO2 and N2, which seemed to be a common
-- older composiiton, but this version adds Argon.
-- The composition is based on what was used in Cruden et al. (2016):
-- Radiative Heating During Mars Science Laboratory Entry: Simulation, Ground Test, and Flight
-- Journal of Thermophysics and Heat Transfer
-- I got the species list from my working co2-with-ions species list and species I added after doing
-- old PITOT code runs with this composition at various conditions to see what the N2 and Ar added.
-- Chris James (c.james4@uq.edu.au) - 12/06/21

model = "CEAGas"

CEAGas = {
  mixtureName = 'mars-c02-n2-ar-with-ions',
  speciesList = {'CO2','C3O2','C2O','CO','CN','NO','NO2','O3','O2','O','C3','C2','C','N2','N','Ar','O2+','C2+','CO+','NO+','C+','C-','O+','O-','N+','Ar+','e-'},
  reactants = {CO2=0.958,N2=0.027,Ar=0.015},
  inputUnits = "moles",
  withIons = true,
  trace = 1.0e-6
}
