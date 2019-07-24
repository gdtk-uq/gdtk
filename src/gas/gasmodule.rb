# gasmodule.rb
# A Ruby wrapper for GasModel and GasState classes.
# PJ 2019-07-24: start of experiment.
#
require 'fiddle'
require 'fiddle/import'

module Gasmodule
  extend Fiddle::Importer
  dlload 'libgasmodule.so'
  extern 'int cwrap_gas_module_init()'
end

puts "Initialize the gas module."
Gasmodule.cwrap_gas_module_init()
