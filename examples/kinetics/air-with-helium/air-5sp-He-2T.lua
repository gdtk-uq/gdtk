-- A 5 species air model with the addition of Helium
-- The reason for Helium is that it might be a contaminant
-- in the mixture, or perhaps directly injected at a boundary.
--
-- RJG, 2022-06-26
--

model = "TwoTemperatureGas"
species = {'N2', 'O2', 'N', 'O', 'NO', 'He'}

-- There are some collision integral pairs that we do not have input data for.
-- These are mostly the interactions with He.
-- This below sets a default for missing collision integrals.
options = {
   CI_default = "N2:N2"
}
