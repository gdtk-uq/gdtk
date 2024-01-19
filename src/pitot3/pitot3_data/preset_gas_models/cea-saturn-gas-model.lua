-- Saturn entry gas composition based on the outer atmospheric compoisition of Saturn
-- From Conrath, B.J., and Gautier, D., 
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
-- which was a calculation that I did to get the related Uranus entry gas composition
-- Chris James (c.james4@uq.edu.au) - 19/01/24

model = "CEAGas"

CEAGas = {
  mixtureName = 'saturn',
  speciesList = {"H2","H","He","H+","He+","e-"},
  reactants = {H2=0.89, He=0.11},
  inputUnits = "moles",
  withIons = true,
  trace = 1.0e-10
}
