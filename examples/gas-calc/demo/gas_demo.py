#! /usr/bin/env python3
# gas_demo.py
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
from gdtk.gas import GasModel, GasState

gmodel = GasModel("ideal-air-gas-model.lua")
print("gmodel=", gmodel)
print("n_species=", gmodel.n_species, ", n_modes=", gmodel.n_modes)
print("species_names=", gmodel.species_names)
print("mol_masses=", gmodel.mol_masses)

gs = GasState(gmodel)
print("freshly minted gs=", gs)
gs.rho=1.1; gs.p=1.0e5; gs.T=300.0; gs.u=1.44e6; gs.massf=[1.0]
print("after setting some values")
print("  gs.rho=%g p=%g T=%g u=%g massf=%s a=%g k=%g mu=%g" %
      (gs.rho, gs.p, gs.T, gs.u, gs.massf, gs.a, gs.k, gs.mu))
gmodel.update_thermo_from_pT(gs) # the way that we do the update in D
gmodel.update_sound_speed(gs) # not necessary from Python but is available
gmodel.update_trans_coeffs(gs)
print("after update thermo from pT")
print("  gs.rho=%g p=%g T=%g u=%g massf=%s a=%g k=%g mu=%g" %
      (gs.rho, gs.p, gs.T, gs.u, gs.massf, gs.a, gs.k, gs.mu))
gs.p=3000.0; gs.T=99.0; gs.massf={'air':1.0}
gs.update_thermo_from_rhou() # update another way
gs.update_trans_coeffs()
print("after update thermo from rhou")
print("  gs.rho=%g p=%g T=%g u=%g massf=%s a=%g k=%g mu=%g" %
      (gs.rho, gs.p, gs.T, gs.u, gs.massf, gs.a, gs.k, gs.mu))

print("Some derived properties")
print("gs.Cv=%g Cp=%g R=%g enthalpy=%g entropy=%g molecular_mass=%g" %
      (gs.Cv, gs.Cp, gs.R, gs.enthalpy, gs.entropy, gs.molecular_mass))
