model = "EquilibriumGas"

-- Note that this equilibrium-chemistry model requires
-- a fully-specified thermally-perfect 1T gas model.

EquilibriumGas = {
  mixtureName = 'O2N2Ions',
  tpGasEqFile = 'tp-n2-o2-ions.lua',
  reactants = {N2=0.79, O2=0.21},
  inputUnits = "moles",
  trace = 1.0e-6
}
