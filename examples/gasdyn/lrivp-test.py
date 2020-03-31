# lrivp-test.py
#
# $ prep-gas ideal-air.inp ideal-air-gas-model.lua
# $ python3 lrivp-test.py
#
# PJ, 2020-03-30 adapted from osher_riemann-test.py
# 
import math
def approxEqual(a, b):
    result = math.isclose(a, b, rel_tol=1.0e-2, abs_tol=1.0e-5)
    print("a=",a, "b=",b, "rel=",(a-b)/b, "abs=",a-b, "result=",result) 
    return result
from eilmer.gas import GasModel, GasState, GasFlow
from eilmer.zero_solvers import secant

# Set up something like Sod's shock tube.
gmodel = GasModel('ideal-air-gas-model.lua')
state1 = GasState(gmodel)
state1.p = 10.0e3 # Pa
state1.T = 278.8 # K
state1.update_thermo_from_pT()
print("state1: %s" % state1)
state4 = GasState(gmodel)
state4.p = 100.0e3 # Pa
state4.T = 348.4 # K
state4.update_thermo_from_pT()
print("state4: %s" % state4)

print("Solve Lagrangian flavour of Riemann problem")
flow = GasFlow(gmodel)
pstar, wstar = flow.lrivp(state4, state1, 0.0, 0.0)
print("pstar=%g wstar=%g" % (pstar, wstar))

print("Solve piston-at-left problem (given contact-surface speed)")
pstar = flow.piston_at_left(state1, 0.0, 293.4)
print("pstar=%g" % pstar)

print("Solve piston-at-right problem (given contact-surface speed)")
pstar = flow.piston_at_right(state4, 0.0, 293.4)
print("pstar=%g" % pstar)

print("Solve again using the state-to-state functions")
# as used in the classic shock tube analysis script.
states = []
for i in range(5): states.append(GasState(gmodel))
states[0].p = 1.0e4; states[0].T = 278.8 # spare, to be used later
states[1].p = 1.0e4; states[1].T = 278.8 # driven tube, initial
states[2].p = 1.0e4; states[2].T = 278.8 # intermediate, post-shock
states[3].p = 1.0e4; states[3].T = 278.8 # intermediate, post-expansion
states[4].p = 1.0e5; states[4].T = 348.4 # driver tube, initial
for gs in states: gs.update_thermo_from_pT()
#
def error_in_velocity(p3p4):
    """
    Compute the velocity mismatch for a given pressure ratio across the expansion.
    Across the expansion, we get a test-gas velocity, v3g.
    """
    p3 = p3p4*states[4].p
    v3g = flow.finite_wave_dp(states[4], 0.0, 'cplus', p3, states[3])
    # Across the contact surface.
    p2 = p3
    print("current guess for p3 and p2=", p2)
    v1s, v2, v2g = flow.normal_shock_p2p1(states[1], p2/states[1].p, states[2])
    return (v3g - v2g)/v3g
#
p3p4 = secant(error_in_velocity, 0.105, 0.11, 1.0e-3)
print("From secant solve: p3/p4=", p3p4)
print("Expanded driver gas:")
p3 = p3p4*states[4].p
v3g = flow.finite_wave_dp(states[4], 0.0, 'cplus', p3, states[3])
print("v3g=", v3g)
print("state 3: %s" % states[3])
print("Shock-processed test gas:")
v1s, v2, v2g = flow.normal_shock_p2p1(states[1], p3/states[1].p, states[2])
print("v1s=", v1s, "v2g=", v2g)
print("state 2: %s" % states[2])
assert (abs(v2g - v3g)/v3g < 1.0e-3), "mismatch in velocities"

# Check for differences between the solution methods.
# Given the approximate nature of the Osher-type solution
# we will tolerate a couple of percent differences.
assert (abs(v2g - wstar)/v2g < 2.0e-2), "mismatch in contact-surface velocities"
assert (abs(states[2].p - pstar)/states[2].p < 2.0e-2), "mismatch in pressure"
