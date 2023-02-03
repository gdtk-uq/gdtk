# A script to compute the effective Lewis numbers of high
# temperature air using the binary diffusion routines.
#
# Author: Nick Gibbons
# Date: 2022-02-03
#
# To run this calculation:
# $ prep-gas air-5sp-1T.inp air-5sp-1T.lua
# $ python3 lewis_numbers.py

import math
from gdtk.gas import GasModel, GasState

gmodel = GasModel('air-5sp-1T.lua')
gs = GasState(gmodel)
gs.p = 500e3 # Pa
gs.T = 300.0 # K
gs.p = 1.0e5 # Pa
gs.T = 3000.0
gs.massf = {"N2":0.70, "O2":0.10, "N":0.05, "O":0.05, "NO":0.1}

gmodel.update_thermo_from_pT(gs)
gmodel.update_trans_coeffs(gs)
D = gmodel.binary_diffusion_coefficients(gs)
molef = gmodel.massf2molef(gs.massf)
cp = gmodel.Cp(gs)

print("Problem Description: ")
print("    T: {}".format(gs.T))
print("    p: {}".format(gs.p))
print("    cp: {}".format(cp))
print("    rho: {}".format(gs.rho))
print("    k: {}".format(gs.k))
print("    molef: {}".format(molef))

# Take the binary diffusion coefficients and compute a single 
# average diffusion coefficient for each species.
# Code from mass_diffusion.d
Dav = [0.0]*gmodel.n_species
for i,Di in enumerate(D):
    sum = 0.0
    for j,Dij in enumerate(Di):
        if i!=j:
            molefj = molef[j]
            sum = sum + molefj/Dij
    molefi = molef[i]
    Dav[i] = (1.0 - molefi)/sum

print("---------------------------------------")
print("Species Specific Lewis Numbers: ")
for i,Davi in enumerate(Dav):
    Lei = gs.rho*Davi*cp/gs.k
    print("    {}: {:.6f}".format(gmodel.species_names[i], Lei))

print("Done.")
