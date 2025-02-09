model = "TwoTemperatureGas"
-- Revert to including the He species when the collision integral data is resolved.
-- species = {'CO2', 'O2', 'CO', 'C', 'O', 'He'}
-- For the moment, build the simulation without He
species = {'CO2', 'O2', 'CO', 'C', 'O'}
options = {
  ci_database = "wright"
}
