#! /usr/bin/env python3
# analytic_he_air_eq.py
#
# PJ, 2020-06-14, adapted from analytic_he_n2.py script.
#     2023-01-27, Tamara's nominal condition.

from gdtk.gas import GasModel, GasState, GasFlow
from gdtk.numeric.zero_solvers import secant

print("Compute the flow conditions expected in shock tube.")
#
print("shock-tube fill conditions with helium driving eq-air")
gm_he = GasModel('ideal-helium-gas-model.lua')
gm_air = GasModel('cea-air5species-gas-model.lua')
flow_air = GasFlow(gm_air)
flow_he = GasFlow(gm_he)
states = []
states.append(GasState(gm_he))  # [0] to interpolate inside expansion
states.append(GasState(gm_air))  # [1] driven
states.append(GasState(gm_air))  # [2]
states.append(GasState(gm_he))  # [3]
states.append(GasState(gm_he))  # [4] driver
states[0].p = 0.80e6; states[0].T = 800.0 # spare, to be used later
states[1].p = 290.0;  states[1].T = 293.0  # driven tube, initial
states[2].p = 290.0;  states[2].T = 293.0  # intermediate, post-shock
states[3].p = 0.80e6; states[3].T = 800.0 # intermediate, post-expansion
states[4].p = 0.80e6; states[4].T = 800.0 # driver tube, initial
for gs in states: gs.update_thermo_from_pT()

print("state 4: %s" % states[4])
print("state 1: %s" % states[1])
print("state 1 cea data: %s", states[1].ceaSavedData)
#
# For the unsteady expansion of the driver gas, regulation of the amount
# of expansion is determined by the shock-processed test gas.
# Across the contact surface between these gases, the pressure and velocity
# have to match so we set up some trials of various pressures and check
# that velocities match.
def error_in_velocity(p3p4):
    """
    Compute the velocity mismatch for a given pressure ratio across the expansion.
    Across the expansion, we get a test-gas velocity, v3g.
    """
    p3 = p3p4*states[4].p
    v3g = flow_he.finite_wave_dp(states[4], 0.0, 'cplus', p3, states[3])
    # Across the contact surface.
    p2 = p3
    print("current guess for p3 and p2=", p2)
    v1s, v2, v2g = flow_air.normal_shock_p2p1(states[1], p2/states[1].p, states[2])
    return (v3g - v2g)/v3g
#
p3p4 = secant(error_in_velocity, 0.105, 0.11, 1.0e-3)
print("From secant solve: p3/p4=", p3p4)
print("Expanded driver gas:")
p3 = p3p4*states[4].p
v3g = flow_he.finite_wave_dp(states[4], 0.0, 'cplus', p3, states[3])
print("v3g=", v3g)
print("state 3: %s" % states[3])
print("Shock-processed test gas:")
v1s, v2, v2g = flow_air.normal_shock_p2p1(states[1], p3/states[1].p, states[2])
print("v1s=", v1s, "v2g=", v2g)
print("state 2: %s" % states[2])
print("state 2 cea data: %s", states[2].ceaSavedData)
assert (abs(v2g - v3g)/v3g < 1.0e-3), "mismatch in velocities"
#
# Make a record for plotting against the simulation data.
# We reconstruct the expected data along a tube -3.0 <= x <= 8.0
# at t=2.8ms, where the diaphragm is at x=0.0.
x_centre = 0.0 # metres
t = 2.80e-3 # seconds
f = open('analytic_eq.data', 'w')
f.write('# 1:x(m)  2:rho(kg/m**3) 3:p(Pa) 4:T(K) 5:V(m/s)\n')
print('Left end')
x = x_centre - states[4].a * t - 1.0 # A bit to the left of the unsteady expansion.
f.write('%g %g %g %g %g\n' % (x, states[4].rho, states[4].p, states[4].T, 0.0))
print('Upstream head of the unsteady expansion.')
x = x_centre - states[4].a * t
f.write('%g %g %g %g %g\n' % (x, states[4].rho, states[4].p, states[4].T, 0.0))
print('The unsteady expansion in n steps.')
n = 100
dp = (states[3].p - states[4].p) / n
states[0].copy_values(states[4])
s = gm_he.entropy(states[4])
v = 0.0
p = states[0].p
for i in range(n):
    rhoa = states[0].rho * states[0].a
    dv = -dp / rhoa
    v = v + dv
    p = p + dp
    states[0].p = p
    states[0].update_thermo_from_ps(s)
    x = x_centre + t * (v - states[0].a)
    f.write('%g %g %g %g %g\n' % (x, states[0].rho, states[0].p, states[0].T, v))
print('Downstream tail of expansion.')
x = x_centre + t * (v3g - states[3].a)
f.write('%g %g %g %g %g\n' % (x, states[3].rho, states[3].p, states[3].T, v3g))
print('Contact surface.')
x = x_centre + t * v3g
f.write('%g %g %g %g %g\n' % (x, states[3].rho, states[3].p, states[3].T, v3g))
x = x_centre + t * v2g
f.write('%g %g %g %g %g\n' % (x, states[2].rho, states[2].p, states[2].T, v2g))
print('Shock front')
x = x_centre + t * v1s
f.write('%g %g %g %g %g\n' % (x, states[2].rho, states[2].p, states[2].T, v2g))
f.write('%g %g %g %g %g\n' % (x, states[1].rho, states[1].p, states[1].T, 0.0))
print('Right end')
x = x_centre + t * v1s + 1.0 # A bit beyond the shock location.
f.write('%g %g %g %g %g\n' % (x, states[1].rho, states[1].p, states[1].T, 0.0))
f.close()
