-- fvreactor.lua
-- A simple fixed-volume reactor.
-- PJ & RJG 2018-04-21.

gm = GasModel:new{"nitrogen-2sp.lua"}
chemUpdate = ChemistryUpdate:new{filename="chem.lua", gasmodel=gm}

Q = GasState:new{gm}
Q.p = 1.0e5 -- Pa
Q.T = 4000.0 -- degree K
molef = {N2=2/3, N=1/3}
Q.massf = gm:molef2massf(molef)
gm:updateThermoFromPT(Q)

tFinal = 200.0e-6 -- s
t = 0.0
dt = 1.0e-6
dtSuggest = 1.0e-11
print("# Start integration")
f = assert(io.open("fvreactor.data", 'w'))
f:write('# 1:t(s)  2:T(K)  3:p(Pa)  4:massf_N2  5:massf_N  6:X_N2  7:X_N\n')
f:write(string.format("%10.3e %10.3f %10.3e %20.12e %20.12e %20.12e %20.12e\n",
                      t, Q.T, Q.p, Q.massf.N2, Q.massf.N, molef.N2, molef.N))
while t <= tFinal do
   dtSuggest = chemUpdate:updateState(Q, dt, dtSuggest, gm)
   t = t + dt
   -- dt = dtSuggest -- uncomment this to get quicker stepping
   gm:updateThermoFromRHOE(Q)
   conc = gm:massf2conc(Q) -- FIX-ME -- Why did I change to concentrations?
   f:write(string.format("%10.3e %10.3f %10.3e %20.12e %20.12e %20.12e %20.12e\n",
                         t, Q.T, Q.p, Q.massf.N2, Q.massf.N, conc.N2, conc.N))
end
f:close()
print("# Done.")
