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

Q = GasState.new(gmodel)
puts "freshly minted Q= #{Q}"
Q.rho = 1.1; Q.p = 1.0e5; Q.T = 300.0; Q.u = 1.44e6
puts "after setting some values Q.rho=#{Q.rho}, p=#{Q.p}, T=#{Q.T}, u=#{Q.u}"
gmodel.update_thermo_from_pT(Q) # the way that we do the update in D
puts "after update thermo from pT Q.rho=#{Q.rho}, p=#{Q.p}, T=#{Q.T}, u=#{Q.u}"
Q.p = 3000.0; Q.T=99.0
Q.update_thermo_from_rhou() # update another way
puts "after update thermo from rhou Q.rho=#{Q.rho}, p=#{Q.p}, T=#{Q.T}, u=#{Q.u}"
