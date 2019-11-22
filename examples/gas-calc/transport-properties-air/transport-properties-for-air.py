# A script to compute the viscosity and thermal conductivity
# of air (as a mixture of N2 and O2) from 200 -- 20000 K.
#
# Author: Peter J. and Rowan J. Gollan
# Date: 2019-11-21
#
# To run this script:
# $ prep-gas thermally-perfect-N2-O2.inp thermally-perfect-N2-O2.lua
# $ python3 transport-properties-for-air.py
#
from eilmer.gas import GasModel, GasState

gasModelFile = 'thermally-perfect-N2-O2.lua'
gmodel = GasModel(gasModelFile)

q = GasState(gmodel)
q.p = 1.0e5 # Pa
q.massf = {"N2":0.78, "O2":0.22} # approximation for the composition of air

outputFile = 'trans-props-air.dat'
print("Opening file for writing: %s" % outputFile)
f = open(outputFile, "w")
f.write("#  1:T[K]      2:mu[Pa.s]      3:k[W/(m.K)]\n")

lowT = 200.0
dT = 100.0

for i in range(199):
    q.T = dT*i + lowT
    q.update_thermo_from_pT()
    q.update_trans_coeffs()
    f.write(" %12.6e %12.6e %12.6e\n" % (q.T, q.mu, q.k))

f.close()
print("File closed. Done.")
