# east-shock-tube.py
# Given initial fill conditions and incident shock speed,
# compute incident-shock state (2) the refleted-shock state (5)
# and the slightly-supersonic condition at the nozzle throat.
#
# Run with command like:
# $ python3 east-shock-tube.py > east-shock-tube.transcript
#
# PJ, 2021-11-19
#     2024-04-05: Update and add throat calculation.

from gdtk.gas import GasModel, GasState, GasFlow
from gdtk.numeric.zero_solvers import secant

print("Compute the flow conditions expected in the EAST shock tube.")
#
print("Shock-tube fill conditions with nitrogen")
gm = GasModel('cea-n2-gas-model.lua')
flow = GasFlow(gm)
state1 = GasState(gm)
state1.p = 150.0/760.0*101.325e3
state1.T = 300.0
state1.update_thermo_from_pT()
print("state 1: %s" % state1)
#
print("Incident shock")
v1s = 2.6e3
state2 = GasState(gm)
v2, v2g = flow.normal_shock(state1, v1s, state2)
mf2 = {'N2':state2.ceaSavedData['massf']['N2'],
       'N':state2.ceaSavedData['massf']['N']}
print("state 2: %s" % state2)
print("  massf: %s" % mf2)
print("  v2:", v2)
print("  v2g:", v2g)
#
print("Reflected shock")
state5 = GasState(gm)
vr_b = flow.reflected_shock(state2, v2g, state5)
mf5 = {'N2':state5.ceaSavedData['massf']['N2'],
       'N':state5.ceaSavedData['massf']['N']}
print("state5: %s" % state5)
print("  massf: %s" % mf5)
print("  vr_b=%g" % vr_b)
#
print("Expand to throat condition (Mach 1.001)")
state6 = GasState(gm)
v6 = flow.expand_to_mach(state5, 1.001, state6)
mf6 = {'N2':state6.ceaSavedData['massf']['N2'],
       'N':state6.ceaSavedData['massf']['N']}
print("  v6=%g Mach=%g" % (v6, v6/state6.a))
print("  state6: %s" % state6)
print("  massf: %s" % mf6)
