# t4-reflected-shock-tunnel.py
# Compute the state-to-state processes
# for a particular shot of the the T4 shock tunnel.
#
# $ cp ${DGD_REPO}/src/gas/sample-data/cea-air13species-gas-model.lua .
# $ python3 t4-reflected-shock-tunnel.py
#
# PJ, 2019-11-30
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
print("Initial test gas:")
print("  state1: %s" % state1)

print("normal shock, given shock speed")
vs = 2414.0
print("  vs=%g" % vs)
state2 = GasState(gmodel)
flow = GasFlow(gmodel)
v2, vg = flow.normal_shock(state1, vs, state2)
print("  v2=%g vg=%g" % (v2, vg))
print("  state2: %s" % state2)
assert approxEqual(v2, 361.9), "v2 number after shock fail"
assert approxEqual(vg, 2052.1), "vg number after shock fail"
assert approxEqual(state2.p, 7.314e6), "p2 number after shock fail"
assert approxEqual(state2.T, 2630.0), "T2 number after shock fail"

print("normal shock, given pressure ratio")
p2p1 = 58.516
print("  p2p1=%g" % p2p1)
vs, v2, vg = flow.normal_shock_p2p1(state1, p2p1, state2)
print("  vs=%g v2=%g vg=%g" % (vs, v2, vg))
print("  state2: %s" % state2)
assert approxEqual(vs, 2414.0), "vs number after p2p1 shock fail"
assert approxEqual(v2, 361.9), "v2 number after p2p1 shock fail"
assert approxEqual(vg, 2052.1), "vg number after p2p1 shock fail"

print("reflected shock")
state5 = GasState(gmodel)
vr_b = flow.reflected_shock(state2, vg, state5)
print("  vr_b=%g" % vr_b)
print("  state5: %s" % state5)
assert approxEqual(vr_b, 573.9), "vr_b number after reflected shock fail"
assert approxEqual(state5.p, 59.47e6), "p5 number after reflected shock fail"
assert approxEqual(state5.T, 4551.8), "T5 number after reflected shock fail"

print("Expand from stagnation (with ratio of pressure to match observation)")
state5s = GasState(gmodel)
v5s = flow.expand_from_stagnation(state5, 34.37/59.47, state5s)
print("  v5s=%g Mach=%g" % (v5s, v5s/state5s.a))
print("  state5s: %s" % state5s)
print("  (h5s-h1)=%g" % (state5s.enthalpy - state1.enthalpy))
assert approxEqual(v5s, 1184.7), "v5s number after expand_from_stagnation fail"
assert approxEqual(state5s.p, 34.37e6), "p5s number after expand_from_stagnation fail"
assert approxEqual(state5s.T, 4161.8), "T5s number after expand_from_stagnation fail"

print("Expand to throat condition (Mach 1.0001)")
state6 = GasState(gmodel)
v6 = flow.expand_to_mach(state5s, 1.0001, state6)
print("  v6=%g Mach=%g" % (v6, v6/state6.a))
print("  state6: %s" % state6)
print("  ceaSavedData=%s" % state6.ceaSavedData)
assert approxEqual(v6, 1155.8), "v6 number after expand_to_mach fail"
assert approxEqual(v6/state6.a, 1.0001), "mach number after expand_to_mach fail"
assert approxEqual(state6.p, 19.32e6), "p6 number after expand_to_mach fail"
assert approxEqual(state6.T, 3788.0), "T6 number after expand_to_mach fail"

print("Something like Mach 4 nozzle")
state7 = GasState(gmodel)
v7 = flow.steady_flow_with_area_change(state6, v6, 27.0, state7)
print("  v7=%g Mach=%g" % (v7, v7/state7.a))
print("  state7: %s" % state7)
assert approxEqual(v7, 2950.7), "v7 number after steady_flow_with_area_change fail"
assert approxEqual(v7/state7.a, 4.238), "mach number after steady_flow_with_area_change fail"
assert approxEqual(state7.p, 93.598e3), "p7 number after steady_flow_with_area_change fail"
assert approxEqual(state7.T, 1282.5), "T7 number after steady_flow_with_area_change fail"

print("Total condition")
state8 = GasState(gmodel)
flow.total_condition(state7, v7, state8)
print("  state8: %s" % state8)

print("Pitot condition")
state9 = GasState(gmodel)
flow.pitot_condition(state7, v7, state9)
print("  state9: %s" % state9)
assert approxEqual(state9.p/state8.p, 0.06253), "pitot/total pressure ratio fail" 
