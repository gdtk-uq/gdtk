-- Initialise an ideal gas model
gm = GasModel:new{'sample-data/ideal-air-gas-model.lua'}
-- Try out some of the service functions
assert(gm:nSpecies() == 1)
assert(gm:nModes() == 0)
assert(gm:speciesName(0) == 'air')
molMasses = gm:molMasses()
assert(math.abs(molMasses['air'] - 0.02896) < 1.0e-6)
-- Test thermo evaluations....
Q = GasState:new{gm}
assert(Q.massf.air == 1)
Q.p = 1.0e5
Q.T = 300.0
Q.massf = {air=1.0}
gm:updateThermoFromPT(Q)
assert(math.abs(Q.u - 215327.43439227) < 1.0e-6)
gm:updateThermoFromRHOE(Q)
assert(math.abs(Q.T - 300.0) < 1.0e-6)
assert(math.abs(Q.p - 1.0e5) < 1.0e-6)
assert(Q.massf.air == 1)


