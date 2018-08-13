-- Authors: Pierre M. and Rowan G.
-- Date: 2018-08-12

gm = GasModel:new{"pseudo-species-nitrogen-3sp.lua"}

Q = GasState:new{gm}
Q.p = 1.0e5 -- Pa
Q.T = 4000.0 -- degree K
molef = {}
molef['N2_eX1SIGGPlus_v0'] = 2/3
molef['N_e2p3Minus4S'] = 1/3
Q.massf = gm:molef2massf(molef)
mfN2 = Q.massf['N2_eX1SIGGPlus_v0']
print("mfN2 = ",mfN2)


-- Now split mass fraction in N2 amonst the two pseudo-species N2(v=0) and N2(v=1)
-- according to Boltzmann populations
eV2J = 1.60219e-19 -- eV2J
k_b = 1.3806e-23
Evib_eV = {
N2_eX1SIGGPlus_v0=-9.754,
N2_eX1SIGGPlus_v1=-9.465,
}
function Qvib_N2X(T)
   -- Compute the vibrational partition function
   qvib=0
   for sp,eV in pairs(Evib_eV) do
      qvib = qvib + math.exp(-eV*eV2J/(k_b*T))
   end
   return qvib
end

qvib = Qvib_N2X(Q.T)
for sp,eV in pairs(Evib_eV) do
   Q.massf[sp] = mfN2*math.exp(-eV*eV2J/(k_b*Q.T))/qvib
   print("Q.massf",sp,Q.massf[sp])
end


gm:updateThermoFromPT(Q)

chemUpdate = PseudoSpeciesKinetics:new{reactionsFile='state-specific-N2-diss.lua', gasModel=gm}

tFinal = 200.0e-6 -- s
t = 0.0
dt = 1.0e-6
dtSuggest = 1.0e-11
print("# Start integration")
f = assert(io.open("state-specific-chem.data", 'w'))
f:write('# 1:t(s)  2:T(K)  3:p(Pa)  4:mf_N_e2p3Minus4S  5:mf_N2_eX1SIGGPlus_v0+1\n')
f:write(string.format("%10.3e %10.3f %10.3e %20.12e %20.12e \n",
                      t, Q.T, Q.p, Q.massf.N_e2p3Minus4S, Q.massf.N2_eX1SIGGPlus_v0+Q.massf.N2_eX1SIGGPlus_v1))
while t <= tFinal do
   dtSuggest = chemUpdate:updateState(Q, dt)
   t = t + dt
   gm:updateThermoFromRHOE(Q)
   conc = gm:massf2conc(Q)
   f:write(string.format("%10.3e %10.3f %10.3e %20.12e %20.12e \n",
                         t, Q.T, Q.p, Q.massf.N_e2p3Minus4S, Q.massf.N2_eX1SIGGPlus_v0+Q.massf.N2_eX1SIGGPlus_v1))
end
f:close()
print("# Done.")

