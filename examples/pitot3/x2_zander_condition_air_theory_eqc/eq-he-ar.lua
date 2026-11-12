model = "EquilibriumGas"

-- Note that this equilibrium-chemistry model requires
-- a fully-specified thermally-perfect 1T gas model.

EquilibriumGas = {
  mixtureName = 'HeArEq',
  tpGasEqFile = 'tp-he-ar.lua',
  reactants = {He=0.8, Ar=0.2},
  inputUnits = "moles",
  trace = 1.0e-6
}
