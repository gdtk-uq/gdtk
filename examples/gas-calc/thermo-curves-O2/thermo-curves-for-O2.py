# A script to output Cp and h for O2 over temperature range 200--20000 K.
#
# Author: Peter J. and Rowan J. Gollan
# Date: 2019-11-21
#
# To run this script:
# $ prep-gas O2.inp O2-gas-model.lua
# $ python3 thermo-curves-for-O2.py
#
from gasmodule import GasModel, GasState

gasModelFile = 'O2-gas-model.lua'
gmodel = GasModel(gasModelFile)

q = GasState(gmodel)
q.p = 1.0e5 # Pa
q.massf = {"O2":1.0}

outputFile = 'O2-thermo.dat'
print("Opening file for writing: %s" % outputFile)
f = open(outputFile, "w")
f.write("#  1:T[K]      2:Cp[J/kg/K]     3:h[J/kg]\n")

lowT = 200.0
dT = 100.0

for i in range(199):
    q.T = dT*i + lowT
    q.update_thermo_from_pT()
    f.write(" %12.6e %12.6e %12.6e\n" % (q.T, q.Cp, q.enthalpy))

f.close()
print("File closed. Done.")
