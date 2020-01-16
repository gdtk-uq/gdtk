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

gs = GasState.new(gmodel)
puts "freshly minted gs= #{gs}"
gs.rho=1.1; gs.p=1.0e5; gs.T=300.0; gs.u=1.44e6; gs.massf=[1.0]
puts "after setting some values"
puts "  gs.rho=%g p=%g T=%g u=%g massf=%s a=%g k=%g mu=%g" %
      [gs.rho, gs.p, gs.T, gs.u, gs.massf, gs.a, gs.k, gs.mu]
gmodel.update_thermo_from_pT(gs) # the way that we do the update in D
gmodel.update_sound_speed(gs) # not necessary from Ruby but is available
gmodel.update_trans_coeffs(gs)
puts "after update thermo from pT"
puts "  gs.rho=%g p=%g T=%g u=%g massf=%s a=%g k=%g mu=%g" %
      [gs.rho, gs.p, gs.T, gs.u, gs.massf, gs.a, gs.k, gs.mu]
gs.p = 3000.0; gs.T=99.0; gs.massf={'air'=>1.0}
gs.update_thermo_from_rhou() # update another way
gs.update_trans_coeffs()
puts "after update thermo from rhou"
puts "  gs.rho=%g p=%g T=%g u=%g massf=%s a=%g k=%g mu=%g" %
      [gs.rho, gs.p, gs.T, gs.u, gs.massf, gs.a, gs.k, gs.mu]

puts "Some derived properties"
puts "gs.Cv=%g Cp=%g R=%g enthalpy=%g entropy=%g molecular_mass=%g" %
     [gs.Cv, gs.Cp, gs.R, gs.enthalpy, gs.entropy, gs.molecular_mass]
