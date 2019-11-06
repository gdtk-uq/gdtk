#! /usr/bin/env ruby
# gasmodule_demo.rb
#
# Try out the GasModel and GasState CFFI wrappers.
#
# You will need to "make install" the loadable library
# and then invoke this script from a location where
# the gas-model Lua file is visible, as given by the path below.
#
# PJ 2019-07-24 direct use of FFI
#    2019-07-27 using Ruby wrapper
#
$LOAD_PATH << '~/dgdinst/lib'
require 'gasmodule'

gmodel = GasModel.new("sample-data/ideal-air-gas-model.lua")
puts "gmodel= #{gmodel}"
puts "n_species= #{gmodel.n_species}, n_modes= #{gmodel.n_modes}"
puts "species_names= #{gmodel.species_names}"
puts "mol_masses= #{gmodel.mol_masses}"

Q = GasState.new(gmodel)
puts "freshly minted Q= #{Q}"
Q.rho=1.1; Q.p=1.0e5; Q.T=300.0; Q.u=1.44e6; Q.massf=[1.0]
puts "after setting some values Q.rho= %g, p=%g, T=%g, u=%g, massf=%s" %
      [Q.rho, Q.p, Q.T, Q.u, Q.massf]
gmodel.update_thermo_from_pT(Q) # the way that we do the update in D
puts "after update thermo from pT Q.rho=%g, p=%g, T=%g, u=%g, massf=%s" %
      [Q.rho, Q.p, Q.T, Q.u, Q.massf]
Q.p = 3000.0; Q.T=99.0; Q.massf={'air'=>1.0}
Q.update_thermo_from_rhou() # update another way
puts "after update thermo from rhou Q.rho=%g, p=%g, T=%g, u=%g, massf=%s" %
      [Q.rho, Q.p, Q.T, Q.u, Q.massf]
