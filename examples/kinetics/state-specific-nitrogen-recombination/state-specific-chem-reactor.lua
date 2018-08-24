-- Authors: Pierre M. and Rowan G.
-- Date: 2018-08-12

gm = GasModel:new{"pseudo-species-nitogen-49sp.lua"}

Q = GasState:new{gm}
Q.p = 1.0e5 -- Pa
Q.T = 4000.0 -- degree K
molef = {}
mfN2 = 2/3
mfN = 1-mfN2


-- Now split mass fraction in N2 amonst the two pseudo-species N2(v=0) and N2(v=1)
-- according to Boltzmann populations
eV2J = 1.60219e-19 -- eV2J
k_b = 1.3806e-23
Evib_eV = {
N2_eX1SIGGPlus_v0=-9.754,
N2_eX1SIGGPlus_v1=-9.465,
N2_eX1SIGGPlus_v2=-9.180,
N2_eX1SIGGPlus_v3=-8.899,
N2_eX1SIGGPlus_v4=-8.621,
N2_eX1SIGGPlus_v5=-8.346,
N2_eX1SIGGPlus_v6=-8.075,
N2_eX1SIGGPlus_v7=-7.808,
N2_eX1SIGGPlus_v8=-7.544,
N2_eX1SIGGPlus_v9=-7.283,
N2_eX1SIGGPlus_v10=-7.027,
N2_eX1SIGGPlus_v11=-6.773,
N2_eX1SIGGPlus_v12=-6.524,
N2_eX1SIGGPlus_v13=-6.278,
N2_eX1SIGGPlus_v14=-6.035,
N2_eX1SIGGPlus_v15=-5.796,
N2_eX1SIGGPlus_v16=-5.561,
N2_eX1SIGGPlus_v17=-5.329,
N2_eX1SIGGPlus_v18=-5.100,
N2_eX1SIGGPlus_v19=-4.876,
N2_eX1SIGGPlus_v20=-4.654,
N2_eX1SIGGPlus_v21=-4.437,
N2_eX1SIGGPlus_v22=-4.222,
N2_eX1SIGGPlus_v23=-4.012,
N2_eX1SIGGPlus_v24=-3.805,
N2_eX1SIGGPlus_v25=-3.601,
N2_eX1SIGGPlus_v26=-3.401,
N2_eX1SIGGPlus_v27=-3.205,
N2_eX1SIGGPlus_v28=-3.012,
N2_eX1SIGGPlus_v29=-2.823,
N2_eX1SIGGPlus_v30=-2.637,
N2_eX1SIGGPlus_v31=-2.455,
N2_eX1SIGGPlus_v32=-2.276,
N2_eX1SIGGPlus_v33=-2.101,
N2_eX1SIGGPlus_v34=-1.929,
N2_eX1SIGGPlus_v35=-1.761,
N2_eX1SIGGPlus_v36=-1.597,
N2_eX1SIGGPlus_v37=-1.436,
N2_eX1SIGGPlus_v38=-1.278,
N2_eX1SIGGPlus_v39=-1.125,
N2_eX1SIGGPlus_v40=-0.974,
N2_eX1SIGGPlus_v41=-0.828,
N2_eX1SIGGPlus_v42=-0.685,
N2_eX1SIGGPlus_v43=-0.545,
N2_eX1SIGGPlus_v44=-0.409,
N2_eX1SIGGPlus_v45=-0.276,
N2_eX1SIGGPlus_v46=-0.147,
N2_eX1SIGGPlus_v47=-0.022,
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

Q.massf['N_e2p3Minus4S']=mfN

gm:updateThermoFromPT(Q)

chemUpdate = PseudoSpeciesKinetics:new{gasModel=gm}

tFinal = 1.0e-3 -- s
t = 0.0
dt = 1.0e-7
print("# Start integration")
f = assert(io.open("state-specific-chem.data", 'w'))

-- Write header
char='"t(s)","T(K)","p(Pa)"'
for i=0,gm:nSpecies()-1 do
  spName = gm:speciesName(i)
  char = char .. string.format(',"%s"',spName)
end
char = char .. "\n"
f:write(char)
-- Write initial condition
char = string.format("%10.3e, %10.3f, %10.3e", t, Q.T, Q.p)
for i=0,gm:nSpecies()-1 do
  spName = gm:speciesName(i)
  char = char .. string.format(", %20.12e",Q.massf[spName])
end
char = char .. "\n"
f:write(char)
-- Write current condition
while t+dt <= tFinal do
  t = t + dt
  print("t = ",t)
  dtSuggest = chemUpdate:updateState(Q, dt)
  gm:updateThermoFromRHOE(Q)
  conc = gm:massf2conc(Q)
  char = string.format("%10.3e, %10.3f, %10.3e", t, Q.T, Q.p)
  for i=0,gm:nSpecies()-1 do
    spName = gm:speciesName(i)
    char = char .. string.format(", %20.12e ",Q.massf[spName])
  end
  char = char .. "\n"
  f:write(char)
end
f:close()
print("# Done.")
