-- cea-30-h2-15-o2-55-he-gas-model.lua
-- This is me re-creating some of the HYPULSE detonation driver compositions
-- for Liam Heffernan. 
-- This composition is from Table 3 for the M15 and M21 conditions in the HYPULSE
-- chapter of Lu and Marren's Advanced Hypersonic Test Facilities textbook:
-- Chue et al. Chapter 3: NASA's HYPULSE Facility at GASL - A Dual Mode,
-- Dual Driver Reflected-Shock/Expansion Tunnel
-- I got the required species from 
-- Jachimowski (1988) An Analytical Study of the Hydrogen-Air Reaction
-- Mechanism With Application to Scramjet Combustion
--  Chris James (c.james4@uq.edu.au) 25/08/23

model = "CEAGas"

CEAGas = {
  mixtureName = 'cea-30-h2-15-o2-55-he',
  speciesList = {"H2", "H", "O2", "O", "OH", "H2O", "H2O2", "HO2", "He"},
  reactants = {H2=0.3, O2=0.15, He = 0.55},
  inputUnits = "moles",
  withIons = false,
  trace = 1.0e-10
}
