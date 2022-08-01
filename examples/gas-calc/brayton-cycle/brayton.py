# brayton.py
# Simple Ideal Brayton Cycle using air-standard assumptions.
# Corresponds to Example 9-5 in the 5th Edition of
# Cengal and Boles' thermodynamics text.
#
# To run the script:
# $ prep-gas ideal-air.inp ideal-air-gas-model.lua
# $ prep-gas thermal-air.inp thermal-air-gas-model.lua
# $ python3 brayton.rb
#
# Peter J and Rowan G. 2019-11-21
from gdtk.gas import GasModel, GasState

gasModelFile = "thermal-air-gas-model.lua"
# gasModelFile = "ideal-air-gas-model.lua" # Alternative
gmodel = GasModel(gasModelFile)
if gmodel.n_species == 1:
    print("Ideal air gas model.")
    air_massf = {"air":1.0}
else:
    print("Thermally-perfect air model.")
    air_massf = {"N2":0.78, "O2":0.22}

print("Compute cycle states:")
gs = [] # We will build up a list of gas states
h = [] # and enthalpies.
# Note that we want to use indices consistent with the Lua script,
# so we set up 5 elements but ignore the one with 0 index.
for i in range(5):
    gs.append(GasState(gmodel))
    h.append(0.0)
for i in range(1,5):
    gs[i].massf = air_massf

print("   Start with ambient air")
gs[1].p = 100.0e3; gs[1].T = 300.0
gs[1].update_thermo_from_pT()
s12 = gs[1].entropy
h[1] = gs[1].enthalpy

print("   Isentropic compression with a pressure ratio of 8")
gs[2].p = 8 * gs[1].p
gs[2].update_thermo_from_ps(s12)
h[2] = gs[2].enthalpy

print("   Constant pressure heat addition to T=1300K")
gs[3].p = gs[2].p; gs[3].T = 1300.0
gs[3].update_thermo_from_pT()
h[3] = gs[3].enthalpy
s34 = gs[3].entropy

print("   Isentropic expansion to ambient pressure")
gs[4].p = gs[1].p
gs[4].update_thermo_from_ps(s34)
h[4] = gs[4].enthalpy

print("")
print("State   Pressure Temperature   Enthalpy")
print("             kPa           K      kJ/kg")
print("---------------------------------------")
for i in range(1,5):
    print(" %4d %10.2f  %10.2f %10.2f" %
	  (i, gs[i].p/1000, gs[i].T, h[i]/1000))
print("---------------------------------------")
print("")
print("Cycle performance:")
work_comp_in = h[2] - h[1]
work_turb_out = h[3] - h[4]
heat_in = h[3] - h[2]
rbw = work_comp_in / work_turb_out
eff = (work_turb_out-work_comp_in) / heat_in
print("   turbine work out = %.2f kJ/kg" % (work_turb_out/1000))
print("   back work ratio = %.3f" % (rbw))
print("   thermal_efficiency = %.3f" % (eff))
