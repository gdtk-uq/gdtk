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
  extern 'int gas_model_n_modes(int gm_i)'

  extern 'int gas_state_new(int gm_i)'
  extern 'int gas_state_set_scalar_field(int gs_i, char* field_name, double value)'
  extern 'double gas_state_get_scalar_field(int gs_i, char* field_name)'

  extern 'int gas_model_gas_state_update_thermo_from_pT(int gm_i, int gs_i)'
  extern 'int gas_model_gas_state_update_thermo_from_rhou(int gm_i, int gs_i)'
  extern 'int gas_model_gas_state_update_thermo_from_rhoT(int gm_i, int gs_i)'
  extern 'int gas_model_gas_state_update_thermo_from_rhop(int gm_i, int gs_i)'
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
  def n_modes()
    Gasmodule.gas_model_n_modes(@id)
  end

  def update_thermo_from_pT(gstate)
    flag = Gasmodule.gas_model_gas_state_update_thermo_from_pT(@id, gstate.id)
    if flag < 0 then raise "Oops, could not update thermo from p,T." end
  end
  def update_thermo_from_rhou(gstate)
    flag = Gasmodule.gas_model_gas_state_update_thermo_from_rhou(@id, gstate.id)
    if flag < 0 then raise "Oops, could not update thermo from rho,u." end
  end
  def update_thermo_from_rhoT(gstate)
    flag = Gasmodule.gas_model_gas_state_update_thermo_from_rhoT(@id, gstate.id)
    if flag < 0 then raise "Oops, could not update thermo from rho,T." end
  end
  def update_thermo_from_rhop(gstate)
    flag = Gasmodule.gas_model_gas_state_update_thermo_from_rhop(@id, gstate.id)
    if flag < 0 then raise "Oops, could not update thermo from rho,p." end
  end
end


class GasState
  include Gasmodule
  attr_reader :id
  
  def initialize(gmodel)
    @gmodel = gmodel
    @id = Gasmodule.gas_state_new(gmodel.id)
  end

  def to_s()
    text = "GasState(rho=#{self.rho}"
    text << ", p=#{self.p}"
    text << ", T=#{self.T}"
    text << ", u=#{self.u}"
    text << ", id=#{@id}, gmodel.id=#{@gmodel.id})"
  end
    
  def rho()
    Gasmodule.gas_state_get_scalar_field(@id, "rho")
  end
  def rho=(value)
    flag = Gasmodule.gas_state_set_scalar_field(@id, "rho", value)
    if flag < 0 then raise "Oops, could not set density." end
  end
    
  def p()
    Gasmodule.gas_state_get_scalar_field(@id, "p")
  end
  def p=(value)
    flag = Gasmodule.gas_state_set_scalar_field(@id, "p", value)
    if flag < 0 then raise "Oops, could not set pressure." end
  end
    
  def T()
    Gasmodule.gas_state_get_scalar_field(@id, "T")
  end
  def T=(value)
    flag = Gasmodule.gas_state_set_scalar_field(@id, "T", value)
    if flag < 0 then raise "Oops, could not set temperature." end
  end
    
  def u()
    Gasmodule.gas_state_get_scalar_field(@id, "u")
  end
  def u=(value)
    flag = Gasmodule.gas_state_set_scalar_field(@id, "u", value)
    if flag < 0 then raise "Oops, could not set internal energy." end
  end

  def update_thermo_from_pT()
    @gmodel.update_thermo_from_pT(self)
  end
  def update_thermo_from_rhou()
    @gmodel.update_thermo_from_rhou(self)
  end
  def update_thermo_from_rhoT()
    @gmodel.update_thermo_from_rhoT(self)
  end
  def update_thermo_from_rhop()
    @gmodel.update_thermo_from_rhop(self)
  end
end
