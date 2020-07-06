-- fvreactor.lua
-- A simple fixed-volume reactor.
-- PJ & RJG 2018-04-21.

gm = GasModel:new{"nitrogen-2sp.lua"}
chemUpdate = ChemistryUpdate:new{gasmodel=gm, filename="chem.lua"}

gs = GasState:new{gm}
gs.p = 1.0e5 -- Pa
gs.T = 4000.0 -- degree K
molef = {N2=2/3, N=1/3}
gs.massf = gm:molef2massf(molef)
gm:updateThermoFromPT(gs)
conc = gm:massf2conc(gs)

tFinal = 300.0e-6 -- s
t = 0.0
dt = 1.0e-6
dtSuggest = 1.0e-11
print("# Start integration")
f = assert(io.open("fvreactor.data", 'w'))
f:write('# 1:t(s)  2:T(K)  3:p(Pa)  4:massf_N2  5:massf_N  6:conc_N2  7:conc_N\n')
f:write(string.format("%10.3e %10.3f %10.3e %20.12e %20.12e %20.12e %20.12e\n",
                      t, gs.T, gs.p, gs.massf.N2, gs.massf.N, conc.N2, conc.N))
while t <= tFinal do
   dtSuggest = chemUpdate:updateState(gs, dt, dtSuggest, gm)
   t = t + dt
   -- dt = dtSuggest -- uncomment this to get quicker stepping
   gm:updateThermoFromRHOU(gs)
   conc = gm:massf2conc(gs)
   f:write(string.format("%10.3e %10.3f %10.3e %20.12e %20.12e %20.12e %20.12e\n",
                         t, gs.T, gs.p, gs.massf.N2, gs.massf.N, conc.N2, conc.N))
end
f:close()
print("# Done.")
