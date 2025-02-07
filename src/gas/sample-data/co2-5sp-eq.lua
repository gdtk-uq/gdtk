model = "EquilibriumGas"

-- Note that this equilibrium-chemistry model requires
-- a fully-specified thermally-perfect 1T gas model.
-- This may be generated from a list of species with the command:
-- $ prep-gas co2-5sp-1T-input.lua co2-5sp-1T.lua

EquilibriumGas = {
  mixtureName = 'co2-5sp',
  tpGasEqFile = 'co2-5sp-1T.lua',
  reactants = {CO2=1.0},
  inputUnits = "moles",
  trace = 1.0e-6
}
