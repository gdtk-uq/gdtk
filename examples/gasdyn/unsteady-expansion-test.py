# unsteady-expansion-test.py
#
# $ cp ${DGD_REPO}/src/gas/sample-data/cea-air13species-gas-model.lua .
# $ python3 unsteady-expansion-test.py
#
# PJ, 2019-12-01
# 
import math
def approxEqual(a, b):
    result = math.isclose(a, b, rel_tol=1.0e-2, abs_tol=1.0e-5)
    # print("a=",a, "b=",b, "rel=",(a-b)/b, "abs=",a-b, "result=",result) 
    return result
from eilmer.gas import GasModel, GasState, GasFlow

print("Unsteady expansion.")
gmodel = GasModel('cea-air13species-gas-model.lua')
state1 = GasState(gmodel)
state1.p = 100.0e3 # Pa
state1.T = 320.0 # K  ideal air, not high T
state1.update_thermo_from_pT()
state1.update_sound_speed()
print("  state1: %s" % state1)
v1 = 0.0
jplus = v1 + 2*state1.a/(1.4-1)
print("  v1=%g jplus=%g" % (v1,jplus))

print("Finite wave process along a cplus characteristic, stepping in pressure.")
state2 = GasState(gmodel)
flow = GasFlow(gmodel)
v2 = flow.finite_wave_dp(state1, v1, "cplus", 60.0e3, state2, 500)
print("  v2=%g" % v2)
print("  state2: %s" % state2)
print("  ideal v2=%g" % (jplus - 2*state2.a/(1.4-1)))
assert approxEqual(v2, 126.2), "velocity after finite_wave_dp fail"
assert approxEqual(state2.p, 60.0e3), "pressure after finite_wave_dp fail"
assert approxEqual(state2.T, 276.5), "temperature after finite_wave_dp fail"

print("Finite wave process along a cplus characteristic, stepping in velocity.")
v2 = flow.finite_wave_dv(state1, v1, "cplus", 125.0, state2)
print("  v2=%g" % v2)
print("  state2: %s" % state2)
print("  ideal v2=%g" % (jplus - 2*state2.a/(1.4-1)))
assert approxEqual(v2, 125.0), "velocity after finite_wave_dv fail"
assert approxEqual(state2.p, 60.3e3), "pressure after finite_wave_dv fail"
assert approxEqual(state2.T, 276.9), "temperature after finite_wave_dv fail"
