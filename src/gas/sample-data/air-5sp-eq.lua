model = "EquilibriumGas"

-- Note that this equilibrium-chemistry model requires
-- an fully-specified thermally-perfect 1T gas model.
-- This may be generated from a list of species with the command:
-- $ prep-gas air-5sp-1T-input.lua air-5sp-1T.lua

EquilibriumGas = {
  mixtureName = 'air5species',
  tpGasEqFile = 'air-5sp-1T.lua',
  reactants = {N2=0.79, O2=0.21},
  inputUnits = "moles",
  trace = 1.0e-6
}
