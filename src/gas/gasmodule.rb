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
  extern 'int gas_model_species_name(int gm_i, int isp, char* name, int n)'
  extern 'int gas_model_mol_masses(int gm_i, double* mm)'

  extern 'int gas_state_new(int gm_i)'
  extern 'int gas_state_set_scalar_field(int gs_i, char* field_name, double value)'
  extern 'int gas_state_get_scalar_field(int gs_i, char* field_name, double* value)'
  extern 'int gas_state_set_array_field(int gs_i, char* field_name, double* values, int n)'
  extern 'int gas_state_get_array_field(int gs_i, char* field_name, double* values, int n)'

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
  attr_reader :species_names
  
  def initialize(file_name)
    @file_name = file_name
    @id = Gasmodule.gas_model_new(file_name)
    nsp = Gasmodule.gas_model_n_species(@id)
    @species_names = []
    buf = Fiddle::Pointer.malloc(32)
    nsp.times do |i|
      Gasmodule.gas_model_species_name(@id, i, buf, 32)
      @species_names << buf.to_s
    end
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
  def mol_masses()
    nsp = Gasmodule.gas_model_n_species(@id)
    mm = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE*nsp)
    Gasmodule.gas_model_mol_masses(@id, mm)
    return mm[0, mm.size].unpack("d*")
  end

  def update_thermo_from_pT(gstate)
    flag = Gasmodule.gas_model_gas_state_update_thermo_from_pT(@id, gstate.id)
    if flag < 0 then raise "could not update thermo from p,T." end
  end
  def update_thermo_from_rhou(gstate)
    flag = Gasmodule.gas_model_gas_state_update_thermo_from_rhou(@id, gstate.id)
    if flag < 0 then raise "could not update thermo from rho,u." end
  end
  def update_thermo_from_rhoT(gstate)
    flag = Gasmodule.gas_model_gas_state_update_thermo_from_rhoT(@id, gstate.id)
    if flag < 0 then raise "could not update thermo from rho,T." end
  end
  def update_thermo_from_rhop(gstate)
    flag = Gasmodule.gas_model_gas_state_update_thermo_from_rhop(@id, gstate.id)
    if flag < 0 then raise "could not update thermo from rho,p." end
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
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gasmodule.gas_state_get_scalar_field(@id, "rho", valuep)
    if flag < 0 then raise "could not get density." end
    return valuep[0, valuep.size].unpack("d")[0]
  end
  def rho=(value)
    flag = Gasmodule.gas_state_set_scalar_field(@id, "rho", value)
    if flag < 0 then raise "could not set density." end
  end
    
  def p()
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gasmodule.gas_state_get_scalar_field(@id, "p", valuep)
    if flag < 0 then raise "could not get pressure." end
    return valuep[0, valuep.size].unpack("d")[0]
  end
  def p=(value)
    flag = Gasmodule.gas_state_set_scalar_field(@id, "p", value)
    if flag < 0 then raise "could not set pressure." end
  end
    
  def T()
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gasmodule.gas_state_get_scalar_field(@id, "T", valuep)
    if flag < 0 then raise "could not get temperature." end
    return valuep[0, valuep.size].unpack("d")[0]
  end
  def T=(value)
    flag = Gasmodule.gas_state_set_scalar_field(@id, "T", value)
    if flag < 0 then raise "could not set temperature." end
  end
    
  def u()
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gasmodule.gas_state_get_scalar_field(@id, "u", valuep)
    if flag < 0 then raise "could not get internal energy." end
    return valuep[0, valuep.size].unpack("d")[0]
  end
  def u=(value)
    flag = Gasmodule.gas_state_set_scalar_field(@id, "u", value)
    if flag < 0 then raise "could not set internal energy." end
  end

  def massf()
    nsp = @gmodel.n_species
    mf = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE*nsp)
    Gasmodule.gas_state_get_array_field(@id, "massf", mf, nsp)
    return mf[0, mf.size].unpack("d*")
  end
  def massf=(mf_given)
    nsp = @gmodel.n_species
    if mf_given.class == [].class then
      mf_array = mf_given
    elsif mf_given.class == {}.class then
      mf_array = []
      @gmodel.species_names.each do |name|
        if mf_given.has_key?(name) then
          mf_array << mf_given[name]
        else
          mf_array << 0.0
        end
      end
    end
    mf_sum = 0.0; mf_array.each do |mfi| mf_sum += mfi end
    if (mf_sum - 1.0).abs > 1.0e-6 then
      raise "mass fractions do not sum to 1."
    end
    # mf = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE*nsp)
    mf = mf_array.pack("d*")
    Gasmodule.gas_state_set_array_field(@id, "massf", mf, nsp)
    return mf_array
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
