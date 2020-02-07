# osher-riemann-test.py
#
# $ prep-gas ideal-air.inp ideal-air-gas-model.lua
# $ python3 osher-riemann-test.py
#
# PJ, 2020-02-07
# 
import math
def approxEqual(a, b):
    result = math.isclose(a, b, rel_tol=1.0e-2, abs_tol=1.0e-5)
    # print("a=",a, "b=",b, "rel=",(a-b)/b, "abs=",a-b, "result=",result) 
    return result
from eilmer.gas import GasModel, GasState, GasFlow

# Set up something like Sod's shock tube.
gmodel = GasModel('ideal-air-gas-model.lua')
state1 = GasState(gmodel)
state1.p = 100.0e3 # Pa
state1.T = 300.0 # K
state1.update_thermo_from_pT()
state1.update_sound_speed()
print("state1: %s" % state1)
state4 = GasState(gmodel)
state4.p = 125.0e3 # Pa
state4.T = 300.0 # K
state4.update_thermo_from_pT()
state4.update_sound_speed()
print("state4: %s" % state4)

# Intermediate states
state2 = GasState(gmodel)
state3 = GasState(gmodel)
state0 = GasState(gmodel)

print("Solve Riemann problem")
flow = GasFlow(gmodel)
pstar, wstar, wL, wR, velX0 = \
    flow.osher_riemann(state4, state1, 0.0, 0.0, state3, state2, state0)
print("pstar=%g wstar=%g wL=%g wR=%g, velX0=%g" % (pstar, wstar, wL, wR, velX0))
print("state2: %s" % state2)
# assert approxEqual(vg, 2052.1), "vg number after shock fail"
# assert approxEqual(state2.p, 7.314e6), "p2 number after shock fail"
# assert approxEqual(state2.T, 2630.0), "T2 number after shock fail"
print("state3: %s" % state3)
# assert approxEqual(vg, 2052.1), "vg number after shock fail"
# assert approxEqual(state2.p, 7.314e6), "p2 number after shock fail"
# assert approxEqual(state2.T, 2630.0), "T2 number after shock fail"
print("state0: %s" % state0)
# assert approxEqual(vg, 2052.1), "vg number after shock fail"
# assert approxEqual(state2.p, 7.314e6), "p2 number after shock fail"
# assert approxEqual(state2.T, 2630.0), "T2 number after shock fail"
