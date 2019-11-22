# Step through a steady isentropic expansion,
# from stagnation condition to sonic condition.
#
# $ prep-gas ideal-air.inp ideal-air-gas-model.lua
# $ python3 isentropic-air-expansion.py
#
# Python port, PJ, 2019-11-21
# 
import math
from eilmer.gas import GasModel, GasState

gmodel = GasModel('ideal-air-gas-model.lua')
q = GasState(gmodel)
q.p = 500e3 # Pa
q.T = 300.0 # K
q.update_thermo_from_pT()
# Compute enthalpy and entropy at stagnation conditions
h0 = q.enthalpy
s0 = q.entropy
# Set up for stepping process
dp = 1.0 # Pa, use 1 Pa as pressure step size
q.p = q.p - dp
mach_tgt = 1.0
# Begin stepping until Mach = mach_tgt
while True:
    q.update_thermo_from_ps(s0)
    h1 = q.enthalpy
    v1 = math.sqrt(2*(h0 - h1))
    q.update_sound_speed()
    m1 = v1/q.a
    if m1 >= mach_tgt:
        print("Stopping at Mach=%g" % m1)
        break
    q.p = q.p - dp

print("Gas properties at sonic point are:")
print("p=%g T=%g" % (q.p, q.T))
