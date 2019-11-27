-- brayton.lua
-- Simple Ideal Brayton Cycle using air-standard assumptions.
-- Corresponds to Example 9-5 in the 5th Edition of
-- Cengal and Boles' thermodynamics text.
--
-- To run the script:
-- $ prep-gas ideal-air.inp ideal-air-gas-model.lua
-- $ prep-gas thermal-air.inp thermal-air-gas-model.lua
-- $ gas-calc brayton.lua
--
-- Peter J and Rowan G. 2016-12-19

gasModelFile = "thermal-air-gas-model.lua"
-- gasModelFile = "ideal-air-gas-model.lua" -- Alternative 
gmodel = GasModel:new{gasModelFile}
if gmodel:nSpecies() == 1 then
   print("Ideal air gas model.")
   air_massf = {air=1.0}
else
   print("Thermally-perfect air model")
   air_massf = {N2=0.78, O2=0.22}
end

print("Compute cycle states:")
Q = {} -- We will build up a table of gas states
h = {} -- and enthalpies.
for i=1,4 do
   Q[i] = GasState:new{gmodel}
   Q[i].massf = air_massf
   h[i] = 0.0
end

print("   Start with ambient air")
Q[1].p = 100.0e3; Q[1].T = 300.0
gmodel:updateThermoFromPT(Q[1])
s12 = gmodel:entropy(Q[1])
h[1] = gmodel:enthalpy(Q[1])

print("   Isentropic compression with a pressure ratio of 8")
Q[2].p = 8 * Q[1].p
gmodel:updateThermoFromPS(Q[2], s12)
h[2] = gmodel:enthalpy(Q[2])

print("   Constant pressure heat addition to T=1300K")
Q[3].p = Q[2].p; Q[3].T = 1300.0
gmodel:updateThermoFromPT(Q[3])
h[3] = gmodel:enthalpy(Q[3])
s34 = gmodel:entropy(Q[3])

print("   Isentropic expansion to ambient pressure")
Q[4].p = Q[1].p
gmodel:updateThermoFromPS(Q[4], s34)
h[4] = gmodel:enthalpy(Q[4])

print("")
print("State   Pressure Temperature   Enthalpy")
print("             kPa           K      kJ/kg")
print("---------------------------------------")
for i=1,4 do
   print(string.format(" %4d %10.2f  %10.2f %10.2f",
		       i, Q[i].p/1000, Q[i].T, h[i]/1000))
end
print("---------------------------------------")
print("")
print("Cycle performance:")
work_comp_in = h[2] - h[1]
work_turb_out = h[3] - h[4]
heat_in = h[3] - h[2]
rbw = work_comp_in / work_turb_out
eff = (work_turb_out-work_comp_in) / heat_in
print(string.format("   turbine work out = %.2f kJ/kg", work_turb_out/1000))
print(string.format("   back work ratio = %.3f", rbw))
print(string.format("   thermal_efficiency = %.3f", eff))
