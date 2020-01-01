-- ionizing-argon-reactor.lua
-- Ionizing argon single-temperature reaction to suit Daniel Smith's sims.
--
-- David Gildfind's argon reaction.
-- Peter J. 2019-12-31
--
-- 
Config{
   odeStep = {method='rkf', errTol=1.0e-6},
   tempLimits = {lower=300.0, upper=35000.0},
   tightTempCoupling = true
}

Reaction{
	'Ar + M <=> Ar+ + e- + M',
	fr={'Arrhenius', A=2.50E+34, n=-3.82, C=181700.0},
	label='r22'
}

