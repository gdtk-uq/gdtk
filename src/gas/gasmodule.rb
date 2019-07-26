# gasmodule.rb
# A Ruby wrapper for GasModel and GasState classes.
# PJ 2019-07-24: start of experiment.
#    2019-07-27: add Ruby wrapping code
#
require 'fiddle'
require 'fiddle/import'

module Gasmodule
  extend Fiddle::Importer
  dlload 'libgasmodule.so'
  extern 'int cwrap_gas_module_init()'

  extern 'int gas_model_new(char* file_name)'
  extern 'int gas_model_n_species(int gm_i)'

  extern 'int gas_state_new(int gm_i)'
  extern 'int gas_state_set_scalar_field(int gs_i, char* field_name, double value)'
  extern 'double gas_state_get_scalar_field(int gs_i, char* field_name)'
end

puts "Initialize the gas module."
Gasmodule.cwrap_gas_module_init()

# Service classes that wrap the C-API in a nice Ruby API...

class GasModel
  include Gasmodule
  attr_reader :id
  
  def initialize(file_name)
    @file_name = file_name
    @id = Gasmodule.gas_model_new(file_name)
  end

  def to_s()
    "GasModel(file=\"#{@file_name}\", id=#{@id})"
  end
  
  def n_species()
    Gasmodule.gas_model_n_species(@id)
  end
end


class GasState
  include Gasmodule
  
  def initialize(gmodel)
    @gmodel = gmodel
    @id = Gasmodule.gas_state_new(gmodel.id)
  end

  def to_s()
    text = "GasState(rho=#{self.rho}"
    text << ", id=#{@id}, gmodel.id=#{@gmodel.id})"
  end
    
  def rho()
    Gasmodule.gas_state_get_scalar_field(@id, "rho")
  end
  def rho=(value)
    flag = Gasmodule.gas_state_set_scalar_field(@id, "rho", value)
    if flag < 0 then raise "Oops, could not set density." end
  end
end
