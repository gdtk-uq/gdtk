-- Saturn entry gas composition including 0.4% CH4 (by moles)
-- which is the gas composition from 
-- Cruden, Brett A., and Augustin C. Tibere-Inglesse. 
-- "Impact of Trace CH4 on Shock Layer Radiation in Outer Planet Entry." 
-- AIAA SCITECH 2024 Forum. 2024.
-- which includes the expected Saturn CH4 fraction from:
-- Moses, Julianne I., et al. 
-- "Photochemistry of Saturn's atmosphere: I. Hydrocarbon chemistry and comparisons with ISO observations." 
-- Icarus 143.2 (2000): 244-298.
-- with the balance then taken from
-- Conrath, B.J., and Gautier, D., 
-- “Saturn Helium Abundance: A Reanalysis of Voyager Measurements,” 
-- Icarus, Vol. 144, No. 1, 2000, pp. 124–134.
-- doi:10.1006/icar.1999.6265
-- which found a helium mass fraction at Saturn from 0.18 to 0.25 
-- (with the rest being H2)
-- from my calculation this is a mole fraction of He from 0.9 to 0.142
-- This is often approximated as somewhere in the middle with a gas composition of 89%H2/11%He (by volume)
-- as was performed in 
-- Palmer et al. (2014) Aeroheating Uncertainties in Uranus and Saturn Entries by the Monte Carlo Method
-- Journal of Spacecraft and Rockets
-- I got the species list from very high temperature (14,000 K) H2/He post-shock calculations using the old PITOT code.
-- which was a calculation that I did to get the related Uranus entry gas composition and then the 
-- C and H chemistry from our Titan chemistry models which include H and C...
-- Chris James (c.james4@uq.edu.au) - 19/01/24

model = "CEAGas"

CEAGas = {
  mixtureName = 'saturn-with-ch4',
  speciesList = {'H2','H','He','CH4', 'CH3', 'CH2', 'CH', 'C2H',
  		'C2', 'C', 'C3','C4','C5',
  		'C2H2,acetylene', 'C4H2,butadiyne', 'C2H2,vinylidene',
  		'H+','He+','C+','e-'},
  reactants = {H2=0.886, He=0.110, CH4=0.004},
  inputUnits = "moles",
  withIons = true,
  trace = 1.0e-10
}

