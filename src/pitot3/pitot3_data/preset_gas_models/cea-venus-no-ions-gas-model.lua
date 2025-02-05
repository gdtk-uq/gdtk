-- Venus (C02,N2) WITHOUT ions CEA backed gas model
-- This composition (96.5%CO2/3.5%N2, by volume) is used in Cruden et al. (2012):
-- Absolute Radiation Measurement in Venus and Mars Entry Conditions
-- Journal of Spacecraft and Rockets
-- They don't say it is by volume / moles in the Cruden paper but that seems right from what I know.
-- This atmosphere also lines up with Seiff, A. (1991). 
-- Atmospheres of Earth, Mars, and Venus, as defined by entry probe experiments. 
-- Journal of Spacecraft and Rockets, 28(3), 265-275.
-- Some old xlabs stuff was done with 97%CO2/3%N2, by volume as some work has said that the N2 is 
-- "about 3%" but I think it is good to get it as right as we can.
-- I got the species list from the Mars entry model without Argon which was already in PITOT3
-- which I found by using the co2-with-ions species list and species I added after doing
-- old PITOT code runs with this composition at various conditions to see what the N2 added.
-- Chris James (c.james4@uq.edu.au) - 12/06/25

model = "CEAGas"

CEAGas = {
  mixtureName = 'venus-with-ions',
  speciesList = {'CO2','C3O2','C2O','CO','CN','NO','NO2','O3','O2','O','C3','C2','C','N2','N'},
  reactants = {CO2=0.965,N2=0.035},
  inputUnits = "moles",
  withIons = false,
  trace = 1.0e-10
}
