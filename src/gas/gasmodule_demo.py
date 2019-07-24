#! /usr/bin/env python3
# gasmodule_demo.py
#
# Try out the GasModel and GasState CFFI wrappers.
#
# You will need to "make install" the loadable library
# and then invoke this script from a location where
# the gas-model Lua file is visible, as given by the path below.
#
# PJ 2019-07-24
#
import gasmodule

gmodel = gasmodule.so.gas_model_new(b"sample-data/ideal-air-gas-model.lua")
print("gmodel=", gmodel)
print("n_species=", gasmodule.so.gas_model_n_species(gmodel))

Q = gasmodule.so.gas_state_new(gmodel)
print("Q=", Q)
flag = gasmodule.so.gas_state_set_scalar_field(Q, b"rho", 1.1)
myrho = gasmodule.so.gas_state_get_scalar_field(Q, b"rho")
print("Q.rho=", myrho)
