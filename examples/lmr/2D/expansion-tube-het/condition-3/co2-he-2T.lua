model = "TwoTemperatureGas"
species = {'CO2', 'O2', 'CO', 'C', 'O', 'He'}
options = {
   ci_database = "wright",
   -- There are some collision integral pairs that we do not have input data for.
   -- These are mostly the interactions with He.
   -- This below sets a default for missing collision integrals.
   CI_default = "O2:O2"
}
