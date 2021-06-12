-- Mars (C02,N2) with ions CEA backed gas model
-- This is a simpler Mars entry gas composition which is normally used in our lab,
-- and was also used in the old PITOT code which just used CO2 and N2.
-- This composition is used in Cruden et al. (2012):
-- Absolute Radiation Measurement in Venus and Mars Entry Conditions
-- Journal of Spacecraft and Rockets
-- and a lot of work from our lab such as Troy Eichmann's PhD:
-- Radiation Measurements in a Simulated Mars Atmosphere (2012)
-- I got the species list from my working co2-with-ions species list and species I added after doing
-- old PITOT code runs with this composition at various conditions to see what the N2 added.
-- Chris James (c.james4@uq.edu.au) - 12/06/21

model = "CEAGas"

CEAGas = {
  mixtureName = 'mars-c02-n2-with-ions',
  speciesList = {'CO2','C3O2','C2O','CO','CN','NO','NO2','O3','O2','O','C3','C2','C','N2','N','O2+','C2+','CO+','NO+','C+','C-','O+','O-','N+','e-'},
  reactants = {CO2=0.96,N2=0.04},
  inputUnits = "moles",
  withIons = false,
  trace = 1.0e-6
}
