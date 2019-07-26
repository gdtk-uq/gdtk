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
puts "n_species= #{gmodel.n_species}"

Q = GasState.new(gmodel)
puts "freshly minted Q= #{Q}"
Q.rho = 1.1
puts "after setting density Q.rho= #{Q.rho}"
