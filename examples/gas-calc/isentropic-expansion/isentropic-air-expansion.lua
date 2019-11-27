-- Step through a steady isentropic expansion,
-- from stagnation condition to sonic condition.
--
-- $ prep-gas ideal-air.inp ideal-air-gas-model.lua
-- $ gas-calc isentropic-air-expansion.lua
--
-- Refresh to align with Python and Ruby port, PJ, 2019-11-21.
--
gmodel = GasModel:new{'ideal-air-gas-model.lua'}
q = GasState:new{gmodel}
q.p = 500e3 -- Pa
q.T = 300.0 -- K
gmodel:updateThermoFromPT(q)
-- Compute enthalpy and entropy at stagnation conditions
h0 = gmodel:enthalpy(q)
s0 = gmodel:entropy(q)
-- Set up for stepping process
dp = 1.0 -- Pa, use 1 Pa as pressure step size
q.p = q.p - dp
mach_tgt = 1.0
-- Begin stepping until Mach = mach_tgt
while true do
   gmodel:updateThermoFromPS(q, s0)
   h1 = gmodel:enthalpy(q)
   v1 = math.sqrt(2*(h0 - h1))
   gmodel:updateSoundSpeed(q)
   m1 = v1/q.a
   if m1 >= mach_tgt then
      print(string.format("Stopping at Mach=%g", m1))
      break
   end
   q.p = q.p - dp
end

print("Gas properties at sonic point are:")
print(string.format("p=%g T=%g", q.p, q.T))


      

