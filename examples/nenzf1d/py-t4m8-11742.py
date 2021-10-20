"""
Copy of job t4m8-11742.yaml, executed using the python interface.

@author: Nick Gibbons
"""

from nenzf1d import Nenzf1d
from numpy import array, zeros

# Conditions matching T4 shot 11742
T1 = 300.0;
molef = array([1.0, 0.0, 0.0])

species =  ['N2', 'O2', 'N', 'O', 'NO']
molef = {'N2': 0.79, 'O2': 0.21}
gas_model_1 = "air-5sp-eq.lua"
gas_model_2 = "air-5sp-1T.lua"
reactions = "air-5sp-1T-reactions.lua"

T1 = 300        # K
p1 = 192.0e3    # Pa
Vs = 2268.0     # m/s
pe = 45.8e6     # Pa
ar = 271.16     # Mach 8 nozzle
pp_ps = 7.01e-3 # Check to make sure this is being used...
C = 0.96        # pPitot/(rho*v^2)

# Define the expanding part of the nozzle as a schedule of diameters with position.
# Values are sampled from M8_COORDS.txt file.
xi = array([0.0000, 5.007e-3, 1.038e-2, 1.998e-2, 5.084e-2, 0.10097, 0.20272, 0.40123, 0.60631, 0.80419, 1.110])
di = array([0.0164, 0.01676, 0.01840, 0.02330, 0.04332, 0.07457, 0.12397, 0.18691,0.22705, 0.25263, 0.27006])

molefa = zeros(len(species))
for k,v in molef.items(): molefa[species.index(k)] = v
molef = molefa

sim = Nenzf1d(species, molef, gas_model_1, gas_model_2, reactions, T1, p1, Vs, xi, di, pe=pe, ar=ar, pp_ps=pp_ps, C=C)
sim.run()
sim.print_exit_condition()
