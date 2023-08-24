# sod-driving-left.py
#
# Run with commands like:
# $ prep-gas ideal-air.inp ideal-air-gas-model.lua
# $ python3 sod-driving-left.py
#
# PJ, 2021-08-09 adapted from classic-shock-tube.py

from gdtk.gas import GasModel, GasState, GasFlow
from gdtk.numeric.zero_solvers import secant

print("Compute the flow conditions expected in the Sod shock tube.")
#
print("shock-tube fill conditions with air driving air")
gm = GasModel('ideal-air-gas-model.lua')
flow = GasFlow(gm)
states = []
for i in range(5): states.append(GasState(gm))
states[0].p = 1.0e4; states[0].T = 278.65 # spare, to be used later
states[1].p = 1.0e4; states[1].T = 278.65 # driven tube, initial
states[2].p = 1.0e4; states[2].T = 278.65 # intermediate, post-shock
states[3].p = 1.0e4; states[3].T = 278.65 # intermediate, post-expansion
states[4].p = 1.0e5; states[4].T = 348.31 # driver tube, initial
for gs in states: gs.update_thermo_from_pT()

print("state 1: %s" % states[1])
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
#
# Make a record for plotting against the Eilmer3 simulation data.
# We reconstruct the expected data along a tube 0.0 <= x <= 1.0
# at t=100us, where the diaphragm is at x=0.5.
x_centre = 0.5 # metres
t = 600.0e-6 # seconds
f = open('analytic.data', 'w')
f.write('# 1:x(m)  2:rho(kg/m**3) 3:p(Pa) 4:T(K) 5:V(m/s)\n')
print('Left end')
x = 1.0
f.write('%g %g %g %g %g\n' % (x, states[4].rho, states[4].p, states[4].T, 0.0))
print('Upstream head of the unsteady expansion.')
x = x_centre + states[4].a * t
f.write('%g %g %g %g %g\n' % (x, states[4].rho, states[4].p, states[4].T, 0.0))
print('The unsteady expansion in n steps.')
n = 100
dp = (states[3].p - states[4].p) / n
states[0].copy_values(states[4])
s = gm.entropy(states[4])
v = 0.0
p = states[4].p
for i in range(n):
    rhoa = states[0].rho * states[0].a
    dv = -dp / rhoa
    v = v + dv
    p = p + dp
    states[0].p = p
    states[0].update_thermo_from_ps(s)
    x = x_centre - t * (v - states[0].a)
    f.write('%g %g %g %g %g\n' % (x, states[0].rho, states[0].p, states[0].T, v))
print('Downstream tail of expansion.')
x = x_centre - t * (v3g - states[3].a)
f.write('%g %g %g %g %g\n' % (x, states[3].rho, states[3].p, states[3].T, v3g))
print('Contact surface.')
x = x_centre - t * v3g
f.write('%g %g %g %g %g\n' % (x, states[3].rho, states[3].p, states[3].T, v3g))
x = x_centre - t * v2g  # should not have moved
f.write('%g %g %g %g %g\n' % (x, states[2].rho, states[2].p, states[2].T, v2g))
print('Shock front')
x = x_centre - t * v1s  # should not have moved
f.write('%g %g %g %g %g\n' % (x, states[2].rho, states[2].p, states[2].T, v2g))
f.write('%g %g %g %g %g\n' % (x, states[1].rho, states[1].p, states[1].T, 0.0))
print('Right end')
x = 0.0
f.write('%g %g %g %g %g\n' % (x, states[1].rho, states[1].p, states[1].T, 0.0))
f.close()
