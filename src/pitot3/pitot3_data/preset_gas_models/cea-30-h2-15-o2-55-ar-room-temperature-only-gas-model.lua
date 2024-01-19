-- cea-300-h2-150-o2-55-ar-gas-model.lua
-- This is me re-creating some of the HYPULSE detonation driver compositions
-- for Liam Heffernan. 
-- This composition is from Table 3 for the M9 conditions in the HYPULSE
-- chapter of Lu and Marren's Advanced Hypersonic Test Facilities textbook:
-- Chue et al. Chapter 3: NASA's HYPULSE Facility at GASL - A Dual Mode,
-- Dual Driver Reflected-Shock/Expansion Tunnel
-- I got the required species from 
-- Jachimowski (1988) An Analytical Study of the Hydrogen-Air Reaction
-- Mechanism With Application to Scramjet Combustion
--  Liam Heffernan (liam.heffernan@uqconnect.edu.au) (13/01/24)

model = "CEAGas"

CEAGas = {
  mixtureName = 'cea-300-h2-150-o2-55-ar',
  speciesList = {"H2", "O2", "Ar"},
  reactants = {H2=0.300, O2=0.150, Ar = 0.55},
  inputUnits = "moles",
  withIons = false,
  trace = 1.0e-10
}
