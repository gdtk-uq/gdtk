# prepare-transient-inflow.py
#
# Run with commands like:
# $ prep-gas ideal-air.inp ideal-air.gas
# $ python3 prepare-transient-inflow.py
#
# PJ, 2024-01-27 adapted from classic-shock-tube.py

from gdtk.gas import GasModel, GasState, GasFlow
from gdtk.numeric.zero_solvers import secant

print("Generate the inflow conditions.")
#
print("shock-tube fill conditions with air driving air")
print("expected to provide a supersonic post-shock flow")
gm = GasModel('ideal-air.gas')
flow = GasFlow(gm)
states = []
for i in range(5): states.append(GasState(gm))
states[0].p = 1.0e4; states[0].T = 278.8 # spare, to be used later
states[1].p = 1.0e3; states[1].T = 278.8 # driven tube, initial
states[2].p = 1.0e4; states[2].T = 278.8 # intermediate, post-shock
states[3].p = 1.0e4; states[3].T = 278.8 # intermediate, post-expansion
states[4].p = 1.0e5; states[4].T = 348.4 # driver tube, initial
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
# Write the transient-inflow data file.
f = open('transient-inflow.data', 'w')
f.write('time p T vel.x vel.y\n')
f.write('%g %g %g %g %g\n' % (0.0, states[1].p, states[1].T, 0.0, 0.0))
f.write('%g %g %g %g %g\n' % (100.0e-6, states[1].p, states[1].T, 0.0, 0.0))
f.write('%g %g %g %g %g\n' % (101.0e-6, states[2].p, states[2].T, v2g, 0.0))
f.close()
