# flux-calculators.py
#
# $ prep-gas ideal-air.inp ideal-air-gas-model.lua
# $ python3 flux-calculators.py
#
# PJ, 2021-04-28
#
import math
def approxEqual(a, b):
    result = math.isclose(a, b, rel_tol=1.0e-2, abs_tol=1.0e-5)
    print("a=",a, "b=",b, "rel=",(a-b)/b, "abs=",a-b, "result=",result)
    return result
from eilmer.gas import GasModel, GasState, GasFlow

# Set up something like Sod's shock tube.
gmodel = GasModel('ideal-air-gas-model.lua')
state1 = GasState(gmodel)
state1.p = 90.0e3 # Pa
state1.T = 278.8 # K
state1.update_thermo_from_pT()
print("state1: %s" % state1)
state4 = GasState(gmodel)
state4.p = 100.0e3 # Pa
state4.T = 348.4 # K
state4.update_thermo_from_pT()
print("state4: %s" % state4)

print("Solve Riemann problem approximately")
flow = GasFlow(gmodel)

print("Osher flux calculator.")
F_mass, F_x_momentum, F_energy = \
    flow.osher_flux(state4, state1, 0.0, 0.0)
print("F_mass=%g F_x_momentum=%g F_energy=%g" % (F_mass, F_x_momentum, F_energy))

print("Roe flux calculator.")
F_mass, F_x_momentum, F_energy = \
    flow.roe_flux(state4, state1, 0.0, 0.0)
print("F_mass=%g F_x_momentum=%g F_energy=%g" % (F_mass, F_x_momentum, F_energy))
