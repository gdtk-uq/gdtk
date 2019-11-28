# normal-shock-test.py
#
# $ cp ${DGD_REPO}/src/gas/sample-data/cea-air13species-gas-model.lua .
# $ python3 normal-shock-test.py
#
# PJ, 2019-11-28
# 
import math
from eilmer.gas import GasModel, GasState, GasFlow

gmodel = GasModel('cea-air13species-gas-model.lua')
state1 = GasState(gmodel)
state1.p = 100.0e3 # Pa
state1.T = 300.0 # K
state1.update_thermo_from_pT()
state1.update_sound_speed()
print("state1: %s" % state1)
print("normal shock, given shock speed")
vs = 2414.0
state2 = GasState(gmodel)
flow = GasFlow(gmodel)
v2, vg = flow.normal_shock(state1, vs, state2)
print("v2=%g vg=%g" % (v2, vg))
print("state2: %s" % state2)
