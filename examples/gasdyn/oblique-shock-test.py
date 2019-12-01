# oblique-shock-test.py
#
# $ cp ${DGD_REPO}/src/gas/sample-data/cea-air13species-gas-model.lua .
# $ python3 oblique-shock-test.py
#
# PJ, 2019-12-01
# 
import math
from eilmer.gas import GasModel, GasState, GasFlow

m1 = 1.5
print("Oblique-shock demo for m1=%g" % m1)

gmodel = GasModel('cea-air13species-gas-model.lua')
state1 = GasState(gmodel)
state1.p = 100.0e3 # Pa
state1.T = 300.0 # K ideal air, not high T
state1.update_thermo_from_pT()
state1.update_sound_speed()
print("state1: %s" % state1)
v1 = m1*state1.a
print("v1=%g" % v1)

beta = 45.0 * math.pi/180.0
print("  given beta(degrees)=%g" % (beta*180/math.pi))
state2 = GasState(gmodel)
flow = GasFlow(gmodel)
theta, v2 = flow.theta_oblique(state1, v1, beta, state2)
print("  theta=%g degrees" % (theta*180/math.pi))
print("  v2=%g" % (v2))
print("  state2: %s" % state2)

print("Oblique shock angle from deflection.")
beta2 = flow.beta_oblique(state1, v1, theta)
print("  beta2(degrees)=%g" % (beta2*180/math.pi))
