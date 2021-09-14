-- N2/O2 with ions CEA backed gas model
-- This is to simulate the synthetic Zero air which is currently used for most X2 experiments.
-- Details about the synthetic air can be found here:
-- https://www.coregas.com.au/gases/air/zero-air
-- Composition is 21%O2+-0.5% with an N2 balance.
-- I got the species list from intuition and various old PITOT runs with this composition.
-- Chris James (c.james4@uq.edu.au) - 12/06/21
-- updated 14/09/21 when I realised that I was mising NO+ in the species list... CMJ

model = "CEAGas"

CEAGas = {
  mixtureName = 'n2-o2-with-ions',
  speciesList = {'N2','O2','NO','N2O','NO2','N','O','NO+','N2+','O2+','N+','N-','O+','O-','e-'},
  reactants = {N2=0.79,O2=0.21},
  inputUnits = "moles",
  withIons = true,
  trace = 1.0e-6
}
