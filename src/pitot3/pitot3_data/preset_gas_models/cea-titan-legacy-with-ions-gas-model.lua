-- cea-titan-legacy-with-ions-gas-model.lua
-- Legacy Titan with ions CEA backed gas model
--
-- I called this the 'legacy' Titan gas model as it seems
-- to come from estimates of the atmosphere of Titan from before
-- the entry of the Huygens probe.
-- I have also added newer Titan atmosphere models (both with and without Ar)
-- So this model is NOT the model which I would recommend using, 
-- I have just added it for completeness as experiments were performed in
-- X2 and X3 during Hadas Porat and Bianca Capra's PhDs using this composition.
-- (Some experiments were done by Aaron Brandis and Carolyn Jacobs with this composition too,
-- I think, but not so much.)
-- The original reference appears to be from
-- Yelle et al. (1997)
-- Engineering Models for Titan’s Atmosphere, 
-- Huygens Science, Payload and Mission, ESA, SP-1177,
-- I got the set of species from Johnston et al. (2019)
-- Features of Afterbody Radiative Heating for Titan Entry
-- AIAA Aviation 2019 Forum
-- along with a bunch of species which came up in various states during
-- old PITOT runs with this with this composition.
-- An interesting thing was that graphite (called 'C(gr') by CEA) was appearing
-- in the freestream, I doubt that that is correct, but the composition was significant.
-- I had to remove it from the species list or I couldn't get the calculations to work with PITOT3,
-- so it isn't there anymore. 
-- I also couldn't fit the final two trace species ('C4H2,butadiyne' and 'C2H2,vinylidene')
-- due to some kind of line length issue? and maybe they were causing a crash, so I didn't add them either
-- Chris James (c.james4@uq.edu.au) - 31/08/21

model = "CEAGas"

CEAGas = {
  mixtureName = 'titan-legacy-with-ions',
  speciesList = {'CH4', 'CH3', 'CH2', 'CH', 'N2', 'C2', 'H2', 
  	         'C2N2','CN', 'NH', 'HCN', 'N', 'C', 'H',
		  'N2+', 'CN+', 'N+','N-','C+', 'H+', 'e-',
		  'C3','C4','C5', 'CCN', 'HNC', 
  	          'C4N2', 'CNC', 'NCN', 'C2H', 'NH+', 'C2H2,acetylene'},
  reactants = {N2=0.95,CH4=0.05},
  inputUnits = "moles",
  withIons = true,
  trace = 1.0e-6
}
