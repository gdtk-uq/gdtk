# normal-shock-test.py
#
# $ cp ${DGD_REPO}/src/gas/sample-data/cea-air13species-gas-model.lua .
# $ python3 normal-shock-test.py
#
# PJ, 2019-11-28
# 
import math
def approxEqual(a, b):
    result = math.isclose(a, b, rel_tol=1.0e-2, abs_tol=1.0e-5)
    # print("a=",a, "b=",b, "rel=",(a-b)/b, "abs=",a-b, "result=",result) 
    return result
from eilmer.gas import GasModel, GasState, GasFlow

gmodel = GasModel('cea-air13species-gas-model.lua')
state1 = GasState(gmodel)
state1.p = 125.0e3 # Pa
state1.T = 300.0 # K
state1.update_thermo_from_pT()
state1.update_sound_speed()
print("state1: %s" % state1)

print("normal shock, given shock speed")
vs = 2414.0
print("vs=%g" % vs)
state2 = GasState(gmodel)
flow = GasFlow(gmodel)
v2, vg = flow.normal_shock(state1, vs, state2)
print("v2=%g vg=%g" % (v2, vg))
print("state2: %s" % state2)
assert approxEqual(v2, 361.9), "v2 number after shock fail"
assert approxEqual(vg, 2052.1), "vg number after shock fail"
assert approxEqual(state2.p, 7.314e6), "p2 number after shock fail"
assert approxEqual(state2.T, 2630.0), "T2 number after shock fail"

print("normal shock, given pressure ratio")
p2p1 = 58.516
print("p2p1=%g" % p2p1)
vs, v2, vg = flow.normal_shock_p2p1(state1, p2p1, state2)
print("vs=%g v2=%g vg=%g" % (vs, v2, vg))
print("state2: %s" % state2)
assert approxEqual(vs, 2414.0), "vs number after p2p1 shock fail"
assert approxEqual(v2, 361.9), "v2 number after p2p1 shock fail"
assert approxEqual(vg, 2052.1), "vg number after p2p1 shock fail"

print("reflected shock")
state5 = GasState(gmodel)
vr_b = flow.reflected_shock(state2, vg, state5)
print("vr_b=%g" % vr_b)
print("state5: %s" % state5)
assert approxEqual(vr_b, 573.9), "vr_b number after reflected shock fail"
assert approxEqual(state5.p, 59.47e6), "p5 number after reflected shock fail"
assert approxEqual(state5.T, 4551.8), "T5 number after reflected shock fail"
