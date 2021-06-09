# conical-shock-test.py
#
# $ cp ${DGD_REPO}/src/gas/sample-data/cea-air13species-gas-model.lua .
# $ python3 conical-shock-test.py
#
# PJ, 2019-12-01
#
import math
def approxEqual(a, b):
    result = math.isclose(a, b, rel_tol=1.0e-2, abs_tol=1.0e-5)
    # print("a=",a, "b=",b, "rel=",(a-b)/b, "abs=",a-b, "result=",result)
    return result
from eilmer.gas import GasModel, GasState, GasFlow

m1 = 1.5
print("Conical-shock demo for m1=%g" % m1)

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
state_c = GasState(gmodel)
flow = GasFlow(gmodel)
theta_c, v_c = flow.theta_cone(state1, v1, beta, state_c)
print("  theta_c=%g degrees" % (theta_c*180/math.pi))
print("  v_c=%g" % (v_c))
print("  state_c: %s" % state_c)

print("Conical shock angle from deflection.")
beta2 = flow.beta_cone(state1, v1, theta_c)
print("  beta2(degrees)=%g" % (beta2*180/math.pi))
assert approxEqual(beta, beta2), "conical shock wave angle fail"
