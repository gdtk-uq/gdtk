-- A 5 species air model with the addition of Helium
-- The reason for Helium is that it might be a contaminant
-- in the mixture, or perhaps directly injected at a boundary.
--
-- RJG, 2022-06-26
--

model = "thermally perfect gas"
species = {'N2', 'O2', 'N', 'O', 'NO', 'He'}
