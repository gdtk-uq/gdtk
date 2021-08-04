-- giant planet with ions CEA backed gas model with increased neon instead of helium
-- using the Stalker substitution
-- This is based off a condition from the work of Yu Liu:
-- Liu et al (2020) Experimental validation of a test gas substitution for simulating 
-- non‑equilibrium giant planet entry conditions in impulse facilities
-- Experiments in Fluids
-- I kept the species the same as the Uranus entry gas model (but replaced He with Ne)

-- Lewis Barltrop (l.barltrop@uq.net.au) - 28/06/21
-- Adapted from cea-giant-planet-h2-80-ne-gas-model.lua

model = "CEAGas"

CEAGas = {
  mixtureName = 'giant-planet-h2-30-ne',
  speciesList = {"H2","H","Ne","H+","Ne+","e-"},
  reactants = {H2=0.7, Ne=0.3},
  inputUnits = "moles",
  withIons = true,
  trace = 1.0e-6
}
