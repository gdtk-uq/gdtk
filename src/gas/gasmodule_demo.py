#! /usr/bin/env python3
# gasmodule_demo.py
#
# Try out the GasModel and GasState CFFI wrappers.
#
# You will need to "make install" the loadable library
# and then invoke this script from a location where
# the gas-model Lua file is visible, as given by the path below.
#
# PJ 2019-07-24 direct use of FFI
#    2019-07-25 using Pythonic wrapper
#
from gasmodule import GasModel, GasState

gmodel = GasModel("sample-data/ideal-air-gas-model.lua")
print("gmodel=", gmodel)
print("n_species=", gmodel.n_species)

Q = GasState(gmodel)
print("freshly minted Q=", Q)
Q.rho = 1.1
print("after setting density Q.rho=", Q.rho)
