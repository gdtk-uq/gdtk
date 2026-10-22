model = "EquilibriumGas"

-- Note that this equilibrium-chemistry model requires
-- a fully-specified thermally-perfect 1T gas model.
-- This may be generated from a list of species with the command:
-- $ prep-gas n2-2sp-1T.inp n2-2sp-1T.lua

EquilibriumGas = {
  mixtureName = 'n22species',
  tpGasEqFile = 'nitrogen-2sp-1T.lua',
  reactants = {N2=1.0},
  inputUnits = "moles",
  trace = 1.0e-6
}
