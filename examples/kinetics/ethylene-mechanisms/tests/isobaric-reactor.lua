-- Author: Rowan J. Gollan
-- Date: 2023-06-11
--
-- This script is used to compute temperture histories of ethylene ignition
-- as a comparison to Fig. 3 in Varatharajan and Williams (2002).
--
-- Reference:
-- Varatharajan and Williams (2002)
-- Ethylene Ignition and Detonation Chemistry, Part 2:
-- Ignition Histories and Reduced Mechanims
-- Journal of Propulsion and Power, 18(2), pp. 352--362
-- 

-- Stoichiometric combustion of ethylene in air
--
--   C2H4 + 3 x (O2 + 3.76 N2) <=> 2 CO2 + 2 H2O + 3 x 3.76 N2
--
--   C2H4 + 3 O2 + 11.28 N2 <=> 2 CO2 + 2 H2O + 11.28 N2
--

specOut = {'C2H4', 'C2H3', 'HO2', 'CH3', 'C2H5', 'C2H4O', 'CO', 'CHO', 'CH2O',
           'H', 'OH', 'O', 'H2O2', 'H2O', 'CH2CHO', 'CH2CO'}

gm = GasModel:new{"vw-2002.gas"}
chem = ChemistryUpdate:new{gasmodel=gm, filename="vw-2002.chem"}

gs = GasState:new{gm}
gs.p = 1.0e5 -- Pa (= 1 bar)
gs.T = 1500 -- K
totalMoles = 1 + 3 + 11.28
molef = {C2H4=1/totalMoles, O2=3/totalMoles, N2=11.28/totalMoles}
gs.massf = gm:molef2massf(molef)
gm:updateThermoFromPT(gs)
conc = gm:massf2conc(gs)
m3tocm3 = (1.0e-2)^3
print("Initial concentration C2H4= ", conc.C2H4*m3tocm3)

tFinal = 1e-3 -- s
t = 0.0
dt = 1.0e-8
dtSuggest = -1.0
print("# Start integration")
f = assert(io.open("isobaric-reactor.data", 'w'))
f:write('# 1:t(s)  2:T(K) ')
offset = 2
for i,sp in ipairs(specOut) do
   f:write(string.format("%d:[%s](mol/cm^3) ", i+offset, sp))
end
f:write("\n")
f:write(string.format("%10.3e %10.3f ", t, gs.T))
for _,sp in ipairs(specOut) do
   f:write(string.format("%12.6e ", conc[sp]*m3tocm3))
end
f:write("\n")

--for i=1,100 do
while t <= tFinal do
   dtSuggest = chem:updateState(gs, dt, dtSuggest, gm)
   t = t + dt
   --dt = dtSuggest
   gm:updateThermoFromPU(gs)
   conc = gm:massf2conc(gs)
   f:write(string.format("%10.3e %10.3f ", t, gs.T))
   for _,sp in ipairs(specOut) do
      f:write(string.format("%12.6e ", conc[sp]*m3tocm3))
   end
   f:write("\n")
end

f:close()
print("# Done.")
   

