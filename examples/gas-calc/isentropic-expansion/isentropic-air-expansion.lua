gmodel = GasModel:new{'ideal-air-gas-model.lua'}
Q = GasState:new{gmodel}
Q.p = 500e3 -- Pa
Q.T = 300.0 -- K
gmodel:updateThermoFromPT(Q)
-- Compute enthalpy and entropy at stagnation conditions
h0 = gmodel:enthalpy(Q)
s0 = gmodel:entropy(Q)
-- Set up for stepping process
dp = 1.0 -- Pa, use 1 Pa as pressure step size
Q.p = Q.p - dp
M_tgt = 1.0
-- Begin stepping until M = M_tgt
while true do
   gmodel:updateThermoFromPS(Q, s0)
   h1 = gmodel:enthalpy(Q)
   v1 = math.sqrt(2*(h0 - h1))
   gmodel:updateSoundSpeed(Q)
   M1 = v1/Q.a
   if M1 >= M_tgt then
      print("Stopping at M= ", M1)
      break
   end
   Q.p = Q.p - dp
end

print("Gas state at sonic conditions are:")
print("p= ", Q.p)
print("T= ", Q.T)

      

