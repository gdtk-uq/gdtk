-- cea-titan-no-ar-with-ions-gas-model.lua
-- Titan with ions CEA backed gas model without Argon added.
--
-- This composition is based on modern (last five years)
-- publications from NASA, namely:
-- Brandis and Cruden (2017)
-- Titan Atmospheric Entry Radiative Heating
-- 47th AIAA Thermophysics Conference
-- Johnston et al. (2019)
-- Features of Afterbody Radiative Heating for Titan Entry
-- AIAA Aviation 2019 Forum
-- from Brandis and Cruden:
-- "Since the entry of the Huygens spacecraft at Titan, it was found that 
-- the CH4 mole fraction was closer to 1.5% with trace amounts of H2 and Ar 
-- (both less than < 1%). An approximate test gas composition of 
-- 1.5% CH4 / 0.5% Ar / 98% N2
-- by mole was used to replicate the Titan atmosphere in Test 61."
-- I got the species list from Johnston et al. (2019)
-- along with a bunch of species which came up in various states during
-- old PITOT runs with this with this composition.
-- IN THIS VERSION THE ARGON HAS BEEN REMOVED AND THE BALANCE REPLACED
-- WITH MORE CH4, AS WAS DONE IN JOHNSTON ET AL. (2019)
-- An interesting thing was that graphite (called 'C(gr') by CEA) was appearing
-- in the freestream, I doubt that that is correct, but the composition was significant.
-- I had to remove it from the species list or I couldn't get the calculations to work with PITOT3,
-- so it isn't there anymore. 
-- I also couldn't fit the final two trace species ('C4H2,butadiyne' and 'C2H2,vinylidene')
-- due to some kind of line length issue? and maybe they were causing a crash, so I didn't add them either
-- Chris James (c.james4@uq.edu.au) - 31/08/21

model = "CEAGas"

CEAGas = {
  mixtureName = 'titan-no-ar-with-ions',
  speciesList = {'CH4', 'CH3', 'CH2', 'CH', 'N2', 'C2', 'H2', 
  	         'C2N2','CN', 'NH', 'HCN', 'N', 'C', 'H',
		  'N2+', 'CN+', 'N+','N-','C+', 'H+', 'e-',
		  'C3','C4','C5', 'CCN', 'HNC', 
  	          'C4N2', 'CNC', 'NCN', 'C2H', 'NH+', 'C2H2,acetylene'},
  reactants = {N2=0.98,CH4=0.02},
  inputUnits = "moles",
  withIons = true,
  trace = 1.0e-6
}
