#! /usr/bin/env ruby
# gas_demo.rb
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
require 'eilmer/gas'

gmodel = GasModel.new("ideal-air-gas-model.lua")
puts "gmodel= #{gmodel}"
puts "n_species= #{gmodel.n_species}, n_modes= #{gmodel.n_modes}"
puts "species_names= #{gmodel.species_names}"
puts "mol_masses= #{gmodel.mol_masses}"

q = GasState.new(gmodel)
puts "freshly minted q= #{q}"
q.rho=1.1; q.p=1.0e5; q.T=300.0; q.u=1.44e6; q.massf=[1.0]
puts "after setting some values"
puts "  q.rho=%g p=%g T=%g u=%g massf=%s a=%g k=%g mu=%g" %
      [q.rho, q.p, q.T, q.u, q.massf, q.a, q.k, q.mu]
gmodel.update_thermo_from_pT(q) # the way that we do the update in D
gmodel.update_sound_speed(q)
gmodel.update_trans_coeffs(q)
puts "after update thermo from pT"
puts "  q.rho=%g p=%g T=%g u=%g massf=%s a=%g k=%g mu=%g" %
      [q.rho, q.p, q.T, q.u, q.massf, q.a, q.k, q.mu]
q.p = 3000.0; q.T=99.0; q.massf={'air'=>1.0}
q.update_thermo_from_rhou() # update another way
q.update_sound_speed()
q.update_trans_coeffs()
puts "after update thermo from rhou"
puts "  q.rho=%g p=%g T=%g u=%g massf=%s a=%g k=%g mu=%g" %
      [q.rho, q.p, q.T, q.u, q.massf, q.a, q.k, q.mu]

puts "Some derived properties"
puts "q.Cv=%g q.Cp=%g q.R=%g q.enthalpy=%g q.entropy=%g q.molecular_mass=%g" %
     [q.Cv, q.Cp, q.R, q.enthalpy, q.entropy, q.molecular_mass]
