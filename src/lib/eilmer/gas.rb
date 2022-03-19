# gas.rb
# A Ruby wrapper for GasModel, GasState and ThermochemicalReactor classes.
# It also wraps the simple 1D flow calculations such as shocks and expansions.
#
# PJ 2019-07-24: start of experiment.
#    2019-07-27: add Ruby wrapping code
#
require 'fiddle'
require 'fiddle/import'

module Gas
  extend Fiddle::Importer
  dlload 'libgas.so'
  extern 'int cwrap_gas_init()'

  extern 'int gas_model_new(char* file_name)'
  extern 'int gas_model_type_str(int gm_i, char* dest_str, int n)'
  extern 'int gas_model_n_species(int gm_i)'
  extern 'int gas_model_n_modes(int gm_i)'
  extern 'int gas_model_species_name(int gm_i, int isp, char* name, int n)'
  extern 'int gas_model_mol_masses(int gm_i, double* mm)'

  extern 'int gas_state_new(int gm_i)'
  extern 'int gas_state_set_scalar_field(int gs_i, char* field_name, double value)'
  extern 'int gas_state_get_scalar_field(int gs_i, char* field_name, double* value)'
  extern 'int gas_state_set_array_field(int gs_i, char* field_name, double* values, int n)'
  extern 'int gas_state_get_array_field(int gs_i, char* field_name, double* values, int n)'
  extern 'int gas_state_get_ceaSavedData_field(int gs_i, char* field_name, double* value)'
  extern 'int gas_state_get_ceaSavedData_massf(int gs_i, char* species_name, double* value)'
  extern 'int gas_state_get_ceaSavedData_species_names(int gs_i, char* dest_str, int n)'
  extern 'int gas_state_copy_values(int gs_to_i, int gs_from_i)'

  extern 'int gas_model_gas_state_update_thermo_from_pT(int gm_i, int gs_i)'
  extern 'int gas_model_gas_state_update_thermo_from_rhou(int gm_i, int gs_i)'
  extern 'int gas_model_gas_state_update_thermo_from_rhoT(int gm_i, int gs_i)'
  extern 'int gas_model_gas_state_update_thermo_from_rhop(int gm_i, int gs_i)'
  extern 'int gas_model_gas_state_update_thermo_from_ps(int gm_i, int gs_i, double s)'
  extern 'int gas_model_gas_state_update_thermo_from_hs(int gm_i, int gs_i, double h, double s)'
  extern 'int gas_model_gas_state_update_sound_speed(int gm_i, int gs_i)'
  extern 'int gas_model_gas_state_update_trans_coeffs(int gm_i, int gs_i)'

  extern 'int gas_model_gas_state_Cv(int gm_i, int gs_i, double* result)'
  extern 'int gas_model_gas_state_Cp(int gm_i, int gs_i, double* result)'
  extern 'int gas_model_gas_state_dpdrho_const_T(int gm_i, int gs_i, double* result)'
  extern 'int gas_model_gas_state_R(int gm_i, int gs_i, double* result)'
  extern 'int gas_model_gas_state_gamma(int gm_i, int gs_i, double* result)'
  extern 'int gas_model_gas_state_Prandtl(int gm_i, int gs_i, double* result)'
  extern 'int gas_model_gas_state_internal_energy(int gm_i, int gs_i, double* result)'
  extern 'int gas_model_gas_state_enthalpy(int gm_i, int gs_i, double* result)'
  extern 'int gas_model_gas_state_entropy(int gm_i, int gs_i, double* result)'
  extern 'int gas_model_gas_state_molecular_mass(int gm_i, int gs_i, double* result)'

  extern 'int gas_model_gas_state_enthalpy_isp(int gm_i, int gs_i, int isp, double* result)'
  extern 'int gas_model_gas_state_entropy_isp(int gm_i, int gs_i, int isp, double* result)'
  extern 'int gas_model_gas_state_gibbs_free_energy_isp(int gm_i, int gs_i, int isp, double* result)'

  extern 'int gas_model_massf2molef(int gm_i, double* massf, double* molef)'
  extern 'int gas_model_molef2massf(int gm_i, double* molef, double* massf)'
  extern 'int gas_model_gas_state_get_molef(int gm_i, int gs_i, double* molef)'
  extern 'int gas_model_gas_state_get_conc(int gm_i, int gs_i, double* conc)'

  extern 'int thermochemical_reactor_new(int gm_i, char* filename1, char* filename2)'
  extern 'int thermochemical_reactor_gas_state_update(int cr_i, int gs_i, double t_interval,
                                                      double* dt_suggest)'

  extern 'int gasflow_shock_ideal(int state1_id, double Vs, int state2_id, int gm_id,
                                  double* results)'
  extern 'int gasflow_normal_shock(int state1_id, double Vs, int state2_id, int gm_id,
                                   double* results, double rho_tol, double T_tol)'
  extern 'int gasflow_normal_shock_1(int state1_id, double Vs, int state2_id, int gm_id,
                                     double* results, double p_tol, double T_tol)'
  extern 'int gasflow_normal_shock_p2p1(int state1_id, double p2p1, int state2_id, int gm_id,
                                        double* results)'
  extern 'int gasflow_reflected_shock(int state2_id, double vg, int state5_id, int gm_id,
                                      double* results)'

  extern 'int gasflow_expand_from_stagnation(int state0_id, double p_over_p0, int state1_id,
                                             int gm_id, double* results)'
  extern 'int gasflow_expand_to_mach(int state0_id, double mach, int state1_id,
                                     int gm_id, double* results)'
  extern 'int gasflow_total_condition(int state1_id, double v1, int state0_id, int gm_id)'
  extern 'int gasflow_pitot_condition(int state1_id, double v1, int state2pitot_id, int gm_id)'
  extern 'int gasflow_steady_flow_with_area_change(int state1_id, double v1, double a2_over_a1,
                                                   int state2_id, int gm_id, double tol, double p2p1_min,
                                                   double* results)'

  extern 'int gasflow_finite_wave_dp(int state1_id, double v1, char* characteristic, double p2,
                                     int state2_id, int gm_id, int steps, double* results)'
  extern 'int gasflow_finite_wave_dv(int state1_id, double v1, char* characteristic, double v2_target,
                                     int state2_id, int gm_id, int steps, double Tmin, double* results)'

  extern 'int gasflow_osher_riemann(int stateL_id, int stateR_id, double velL, double velR,
                                    int stateLstar_id, int stateRstar_id,
                                    int stateX0_id, int gm_id, double* results)'
  extern 'int gasflow_osher_flux(int stateL_id, int stateR_id, double velL, double velR,
                                 int gm_id, double* results)'
  extern 'int gasflow_roe_flux(int stateL_id, int stateR_id, double velL, double velR,
                               int gm_id, double* results)'

  extern 'int gasflow_lrivp(int stateL_id, int stateR_id, double velL, double velR,
                            int gmL_id, int gmR_id, double* wstar, double* pstar)'
  extern 'int gasflow_piston_at_left(int stateR_id, double velR, int gm_id,
                                     double wstar, double* pstar)'
  extern 'int gasflow_piston_at_right(int stateL_id, double velL, int gm_id,
                                      double wstar, double* pstar)'

  extern 'int gasflow_theta_oblique(int state1_id, double v1, double beta,
                                    int state2_id, int gm_id, double* results)'
  extern 'int gasflow_beta_oblique(int state1_id, double v1, double theta,
                                   int gm_id, double* results)'

  extern 'int gasflow_theta_cone(int state1_id, double v1, double beta,
                                 int state_c_id, int gm_id, double dtheta, double* results)'
  extern 'int gasflow_beta_cone(int state1_id, double v1, double theta,
                                int gm_id, double dtheta, double* results)'
  PC_P_atm = 101.325e3
end

Gas.cwrap_gas_init()

# Service classes that wrap the C-API in a nice Ruby API...

class GasModel
  include Gas
  attr_reader :id
  attr_reader :species_names

  def initialize(file_name)
    @file_name = file_name
    @id = Gas.gas_model_new(file_name)
    nsp = Gas.gas_model_n_species(@id)
    @species_names = []
    buf = Fiddle::Pointer.malloc(32)
    nsp.times do |i|
      Gas.gas_model_species_name(@id, i, buf, 32)
      @species_names << buf.to_s
    end
  end

  def to_s()
    "GasModel(file=\"#{@file_name}\", id=#{@id}, species=#{@species_names})"
  end

  def type_str()
    buf = Fiddle::Pointer.malloc(32)
    Gas.gas_model_type_str(@id, buf, 32)
    return buf.to_s
  end
  def n_species()
    Gas.gas_model_n_species(@id)
  end
  def n_modes()
    Gas.gas_model_n_modes(@id)
  end
  def mol_masses()
    nsp = Gas.gas_model_n_species(@id)
    mm = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE*nsp)
    Gas.gas_model_mol_masses(@id, mm)
    return mm[0, mm.size].unpack("d*")
  end

  def update_thermo_from_pT(gstate)
    flag = Gas.gas_model_gas_state_update_thermo_from_pT(@id, gstate.id)
    if flag < 0 then raise "could not update thermo from p,T." end
    update_sound_speed(gstate)
  end
  def update_thermo_from_rhou(gstate)
    flag = Gas.gas_model_gas_state_update_thermo_from_rhou(@id, gstate.id)
    if flag < 0 then raise "could not update thermo from rho,u." end
    update_sound_speed(gstate)
  end
  def update_thermo_from_rhoT(gstate)
    flag = Gas.gas_model_gas_state_update_thermo_from_rhoT(@id, gstate.id)
    if flag < 0 then raise "could not update thermo from rho,T." end
    update_sound_speed(gstate)
  end
  def update_thermo_from_rhop(gstate)
    flag = Gas.gas_model_gas_state_update_thermo_from_rhop(@id, gstate.id)
    if flag < 0 then raise "could not update thermo from rho,p." end
    update_sound_speed(gstate)
  end
  def update_thermo_from_ps(gstate, s)
    flag = Gas.gas_model_gas_state_update_thermo_from_ps(@id, gstate.id, s)
    if flag < 0 then raise "could not update thermo from p,s." end
    update_sound_speed(gstate)
  end
  def update_thermo_from_hs(gstate, h, s)
    flag = Gas.gas_model_gas_state_update_thermo_from_hs(@id, gstate.id, h, s)
    if flag < 0 then raise "could not update thermo from h,s." end
    update_sound_speed(gstate)
  end
  def update_sound_speed(gstate)
    flag = Gas.gas_model_gas_state_update_sound_speed(@id, gstate.id)
    if flag < 0 then raise "could not update sound speed." end
  end
  def update_trans_coeffs(gstate)
    flag = Gas.gas_model_gas_state_update_trans_coeffs(@id, gstate.id)
    if flag < 0 then raise "could not update transport coefficients." end
  end

  def Cv(gstate)
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gas.gas_model_gas_state_Cv(@id, gstate.id, valuep)
    if flag < 0 then raise "could not compute Cv." end
    return valuep[0, valuep.size].unpack("d")[0]
  end
  def Cp(gstate)
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gas.gas_model_gas_state_Cp(@id, gstate.id, valuep)
    if flag < 0 then raise "could not compute Cp." end
    return valuep[0, valuep.size].unpack("d")[0]
  end
  def dpdrho_const_T(gstate)
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gas.gas_model_gas_state_dpdrho_const_T(@id, gstate.id, valuep)
    if flag < 0 then raise "could not compute dpdrho_const_T." end
    return valuep[0, valuep.size].unpack("d")[0]
  end
  def R(gstate)
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gas.gas_model_gas_state_R(@id, gstate.id, valuep)
    if flag < 0 then raise "could not compute R." end
    return valuep[0, valuep.size].unpack("d")[0]
  end
  def gamma(gstate)
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gas.gas_model_gas_state_gamma(@id, gstate.id, valuep)
    if flag < 0 then raise "could not compute gamma." end
    return valuep[0, valuep.size].unpack("d")[0]
  end
  def Prandtl(gstate)
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gas.gas_model_gas_state_Prandtl(@id, gstate.id, valuep)
    if flag < 0 then raise "could not compute Prandtl." end
    return valuep[0, valuep.size].unpack("d")[0]
  end
  def internal_energy(gstate)
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gas.gas_model_gas_state_internal_energy(@id, gstate.id, valuep)
    if flag < 0 then raise "could not compute internal energy." end
    return valuep[0, valuep.size].unpack("d")[0]
  end
  def enthalpy(gstate)
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gas.gas_model_gas_state_enthalpy(@id, gstate.id, valuep)
    if flag < 0 then raise "could not compute enthalpy." end
    return valuep[0, valuep.size].unpack("d")[0]
  end
  def entropy(gstate)
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gas.gas_model_gas_state_entropy(@id, gstate.id, valuep)
    if flag < 0 then raise "could not compute entropy." end
    return valuep[0, valuep.size].unpack("d")[0]
  end
  def molecular_mass(gstate)
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gas.gas_model_gas_state_molecular_mass(@id, gstate.id, valuep)
    if flag < 0 then raise "could not compute molecular mass." end
    return valuep[0, valuep.size].unpack("d")[0]
  end

  def enthalpy_isp(gstate, isp)
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gas.gas_model_gas_state_enthalpy_isp(@id, gstate.id, isp, valuep)
    if flag < 0 then raise "could not compute enthalpy for species." end
    return valuep[0, valuep.size].unpack("d")[0]
  end
  def entropy_isp(gstate, isp)
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gas.gas_model_gas_state_entropy_isp(@id, gstate.id, isp, valuep)
    if flag < 0 then raise "could not compute entropy for species." end
    return valuep[0, valuep.size].unpack("d")[0]
  end
  def gibbs_free_energy_isp(gstate, isp)
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gas.gas_model_gas_state_gibbs_free_energy_isp(@id, gstate.id, isp, valuep)
    if flag < 0 then raise "could not compute gibbs free energy for species." end
    return valuep[0, valuep.size].unpack("d")[0]
  end

  def massf2molef(massf_given)
    nsp = Gas.gas_model_n_species(@id)
    if massf_given.class == [].class then
      massf_array = massf_given
    elsif massf_given.class == {}.class then
      massf_array = []
      @species_names.each do |name|
        if massf_given.has_key?(name) then
          massf_array << massf_given[name]
        else
          massf_array << 0.0
        end
      end
    end
    mf_sum = 0.0; massf_array.each do |mfi| mf_sum += mfi end
    if (mf_sum - 1.0).abs > 1.0e-6 then
      raise "mass fractions do not sum to 1."
    end
    my_massf = massf_array.pack("d*")
    my_molef = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE*nsp)
    Gas.gas_model_massf2molef(@id, my_massf, my_molef)
    return my_molef[0, my_molef.size].unpack("d*")
  end

  def molef2massf(molef_given)
    nsp = Gas.gas_model_n_species(@id)
    if molef_given.class == [].class then
      molef_array = molef_given
    elsif molef_given.class == {}.class then
      molef_array = []
      @species_names.each do |name|
        if molef_given.has_key?(name) then
          molef_array << molef_given[name]
        else
          molef_array << 0.0
        end
      end
    end
    mf_sum = 0.0; molef_array.each do |mfi| mf_sum += mfi end
    if (mf_sum - 1.0).abs > 1.0e-6 then
      raise "mole fractions do not sum to 1."
    end
    my_molef = molef_array.pack("d*")
    my_massf = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE*nsp)
    Gas.gas_model_molef2massf(@id, my_molef, my_massf)
    return my_massf[0, my_massf.size].unpack("d*")
  end
end


class GasState
  include Gas
  attr_reader :id
  attr_reader :gmodel

  def initialize(gmodel)
    @gmodel = gmodel
    @id = Gas.gas_state_new(gmodel.id)
  end

  def to_s()
    text = "GasState(rho=#{self.rho}"
    text << ", p=#{self.p}"
    text << ", T=#{self.T}"
    text << ", u=#{self.u}"
    if @gmodel.n_modes > 0 then
      text << ", T_modes=#{self.T_modes}"
      text << ", u_modes=#{self.u_modes}"
    end
    text << ", a=#{self.a}"
    if @gmodel.n_species > 1 then
      text << ", massf=#{self.massf}"
    end
    text << ", id=#{@id}, gmodel.id=#{@gmodel.id})"
  end

  def rho()
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gas.gas_state_get_scalar_field(@id, "rho", valuep)
    if flag < 0 then raise "could not get density." end
    return valuep[0, valuep.size].unpack("d")[0]
  end
  def rho=(value)
    flag = Gas.gas_state_set_scalar_field(@id, "rho", value)
    if flag < 0 then raise "could not set density." end
  end

  def p()
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gas.gas_state_get_scalar_field(@id, "p", valuep)
    if flag < 0 then raise "could not get pressure." end
    return valuep[0, valuep.size].unpack("d")[0]
  end
  def p=(value)
    flag = Gas.gas_state_set_scalar_field(@id, "p", value)
    if flag < 0 then raise "could not set pressure." end
  end

  def T()
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gas.gas_state_get_scalar_field(@id, "T", valuep)
    if flag < 0 then raise "could not get temperature." end
    return valuep[0, valuep.size].unpack("d")[0]
  end
  def T=(value)
    flag = Gas.gas_state_set_scalar_field(@id, "T", value)
    if flag < 0 then raise "could not set temperature." end
  end

  def u()
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gas.gas_state_get_scalar_field(@id, "u", valuep)
    if flag < 0 then raise "could not get internal energy." end
    return valuep[0, valuep.size].unpack("d")[0]
  end
  def u=(value)
    flag = Gas.gas_state_set_scalar_field(@id, "u", value)
    if flag < 0 then raise "could not set internal energy." end
  end

  def a()
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gas.gas_state_get_scalar_field(@id, "a", valuep)
    if flag < 0 then raise "could not get sound speed." end
    return valuep[0, valuep.size].unpack("d")[0]
  end

  def k()
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gas.gas_state_get_scalar_field(@id, "k", valuep)
    if flag < 0 then raise "could not get conductivity." end
    return valuep[0, valuep.size].unpack("d")[0]
  end

  def mu()
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    flag = Gas.gas_state_get_scalar_field(@id, "mu", valuep)
    if flag < 0 then raise "could not get viscosity." end
    return valuep[0, valuep.size].unpack("d")[0]
  end

  def massf()
    nsp = @gmodel.n_species
    mf = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE*nsp)
    flag = Gas.gas_state_get_array_field(@id, "massf", mf, nsp)
    if flag < 0 then raise "could not get mass-fractions." end
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
    mf = mf_array.pack("d*")
    flag = Gas.gas_state_set_array_field(@id, "massf", mf, nsp)
    if flag < 0 then raise "could not set mass-fractions." end
    return mf_array
  end

  def molef()
    nsp = @gmodel.n_species
    mf = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE*nsp)
    flag = Gas.gas_model_gas_state_get_molef(@gmodel.id, @id, mf)
    if flag < 0 then raise "could not get mole-fractions." end
    return mf[0, mf.size].unpack("d*")
  end
  def molef=(molef_given)
    nsp = @gmodel.n_species
    mf_array = @gmodel.molef2massf(molef_given)
    mf = mf_array.pack("d*")
    flag = Gas.gas_state_set_array_field(@id, "massf", mf, nsp)
    if flag < 0 then raise "could not set mass-fractions from mole-fractions." end
    # At this point, we may not have the mole-fractions as an array.
    return nil
  end

  def conc()
    nsp = @gmodel.n_species
    mc = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE*nsp)
    flag = Gas.gas_model_gas_state_get_conc(@gmodel.id, @id, mc)
    if flag < 0 then raise "could not get concentrations." end
    return mc[0, mc.size].unpack("d*")
  end

  def u_modes()
    n = @gmodel.n_modes
    um = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE*n)
    flag = Gas.gas_state_get_array_field(@id, "u_modes", um, n)
    if flag < 0 then raise "could not get u_modes." end
    return um[0, um.size].unpack("d*")
  end
  def u_modes=(um_given)
    n = @gmodel.n_modes
    if um_given.class == [].class then
      raise "u_modes needs to be supplied as an array."
    end
    um = um_given.pack("d*")
    flag = Gas.gas_state_set_array_field(@id, "u_modes", um, n)
    if flag < 0 then raise "could not set u_modes." end
    return um_given
  end

  def T_modes()
    n = @gmodel.n_modes
    myTm = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE*n)
    flag = Gas.gas_state_get_array_field(@id, "T_modes", myTm, n)
    if flag < 0 then raise "could not get T_modes." end
    return myTm[0, myTm.size].unpack("d*")
  end
  def T_modes=(myTm_given)
    n = @gmodel.n_modes
    if myTm_given.class == [].class then
      raise "T_modes needs to be supplied as an array."
    end
    myTm = myTm_given.pack("d*")
    flag = Gas.gas_state_set_array_field(@id, "T_modes", myTm, n)
    if flag < 0 then raise "could not set T_modes." end
    return myTm_given
  end

  def k_modes()
    n = @gmodel.n_modes
    km = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE*n)
    flag = Gas.gas_state_get_array_field(@id, "k_modes", km, n)
    if flag < 0 then raise "could not get k_modes." end
    return km[0, km.size].unpack("d*")
  end

  def ceaSavedData()
    my_data = {}
    scalar_fields = ["p", "rho", "u", "h", "T", "a",
                     "Mmass", "Rgas", "gamma", "Cp",
                     "s", "k", "mu"]
    valuep = Fiddle::Pointer.malloc(Fiddle::SIZEOF_DOUBLE)
    scalar_fields.each do |name|
      flag = Gas.gas_state_get_ceaSavedData_field(@id, name, valuep)
      if flag < 0 then raise "could not get ceaSavedData field %s."%[name] end
      my_data[name] = valuep[0, valuep.size].unpack("d")[0]
    end
    buf = Fiddle::Pointer.malloc(1024)
    flag = Gas.gas_state_get_ceaSavedData_species_names(@id, buf, 1024)
    if flag < 0 then raise "could not get ceaSavedData species names." end
    cea_species_names = buf.to_s.split("\t")
    massf_data = {}
    cea_species_names.each do |name|
      flag = Gas.gas_state_get_ceaSavedData_massf(@id, name, valuep)
      if flag < 0 then raise "could not get ceaSavedData massf[%s]."%[name] end
      massf_data[name] = valuep[0, valuep.size].unpack("d")[0]
    end
    my_data["massf"] = massf_data
    return my_data
  end

  def copy_values(gstate)
    flag = Gas.gas_state_copy_values(self.id, gstate.id)
    if flag < 0 then raise "could not copy values" end
    return nil
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
  def update_thermo_from_ps(s)
    @gmodel.update_thermo_from_ps(self, s)
  end
  def update_thermo_from_hs(h, s)
    @gmodel.update_thermo_from_hs(self, h, s)
  end
  def update_sound_speed()
    @gmodel.update_sound_speed(self)
  end
  def update_trans_coeffs()
    @gmodel.update_trans_coeffs(self)
  end

  def Cv()
    @gmodel.Cv(self)
  end
  def Cp()
    @gmodel.Cp(self)
  end
  def dpdrho_const_T()
    @gmodel.dpdrho_const_T(self)
  end
  def R()
    @gmodel.R(self)
  end
  def gamma()
    @gmodel.gamma(self)
  end
  def Prandtl()
    @gmodel.Prandtl(self)
  end
  def internal_energy()
    @gmodel.internal_energy(self)
  end
  def enthalpy()
    @gmodel.enthalpy(self)
  end
  def entropy()
    @gmodel.entropy(self)
  end
  def molecular_mass()
    @gmodel.molecular_mass(self)
  end

  def enthalpy_isp(isp)
    @gmodel.enthalpy(self, isp)
  end
  def entropy_isp(isp)
    @gmodel.entropy(self, isp)
  end
  def gibbs_free_energy_isp(isp)
    @gmodel.gibbs_free_energy(self, isp)
  end
end


class ThermochemicalReactor
  include Gas
  attr_reader :id

  def initialize(gmodel, filename1, filename2="")
    @filename1 = filename1
    @filename2 = filename2
    @gmodel = gmodel
    @id = Gas.thermochemical_reactor_new(gmodel.id, filename1, filename2)
  end

  def to_s()
    text = "ThermochemicalReactor(id=#{@id}, gmodel.id=#{@gmodel.id}"
    text << " file1=#{@filename1}, file2=#{@filename2})"
  end

  def update_state(gstate, t_interval, dt_suggest)
    dt_suggestp = [dt_suggest].pack("d*")
    flag = Gas.thermochemical_reactor_gas_state_update(@id, gstate.id,
                                                       t_interval, dt_suggestp)
    if flag < 0 then raise "could not update state." end
    return dt_suggestp[0, dt_suggestp.size].unpack("d")[0]
  end
end

#------------------------------------------------------------------------------------

class GasFlow
  include Gas

  def initialize(gmodel)
    @gmodel = gmodel
  end

  def to_s()
    text = "GasFlow(gmodel.id=#{@gmodel.id})"
  end

  def ideal_shock(state1, vs, state2)
    my_results = [0.0, 0.0].pack("d*")
    flag = Gas.gasflow_shock_ideal(state1.id, vs, state2.id, @gmodel.id, my_results)
    if flag < 0 then raise "failed to compute ideal shock jump." end
    return my_results[0, my_results.size].unpack("dd") # [v2, vg]
  end

  def normal_shock(state1, vs, state2, rho_tol=1.0e-6, t_tol=0.1)
    my_results = [0.0, 0.0].pack("d*")
    flag = Gas.gasflow_normal_shock(state1.id, vs, state2.id, @gmodel.id, my_results,
                                    rho_tol, t_tol)
    if flag < 0 then raise "failed to compute normal shock jump from shock speed." end
    return my_results[0, my_results.size].unpack("dd") # [v2, vg]
  end

  def normal_shock_1(state1, vs, state2, p_tol=0.5, t_tol=0.1)
    my_results = [0.0, 0.0].pack("d*")
    flag = Gas.gasflow_normal_shock_1(state1.id, vs, state2.id, @gmodel.id, my_results,
                                      p_tol, t_tol)
    if flag < 0 then raise "failed to compute normal shock jump from shock speed." end
    return my_results[0, my_results.size].unpack("dd") # [v2, vg]
  end

  def normal_shock_p2p1(state1, p2p1, state2)
    my_results = [0.0, 0.0, 0.0].pack("d*")
    flag = Gas.gasflow_normal_shock_p2p1(state1.id, p2p1, state2.id, @gmodel.id, my_results)
    if flag < 0 then raise "failed to compute normal shock jump from pressure ratio." end
    return my_results[0, my_results.size].unpack("ddd") # [vs, v2, vg]
  end

  def reflected_shock(state2, vg, state5)
    my_results = [0.0].pack("d")
    flag = Gas.gasflow_reflected_shock(state2.id, vg, state5.id, @gmodel.id, my_results)
    if flag < 0 then raise "failed to compute reflected shock." end
    return my_results[0, my_results.size].unpack("d")[0] # vr
  end

  def expand_from_stagnation(state0, p_over_p0, state1)
    my_results = [0.0].pack("d")
    flag = Gas.gasflow_expand_from_stagnation(state0.id, p_over_p0, state1.id,
                                              @gmodel.id, my_results)
    if flag < 0 then raise "failed to compute expansion from stagnation." end
    return my_results[0, my_results.size].unpack("d")[0] # v
  end

  def expand_to_mach(state0, mach, state1)
    my_results = [0.0].pack("d")
    flag = Gas.gasflow_expand_to_mach(state0.id, mach, state1.id,
                                      @gmodel.id, my_results)
    if flag < 0 then raise "failed to compute expansion to mach number." end
    return my_results[0, my_results.size].unpack("d")[0] # v
  end

  def total_condition(state1, v1, state0)
    flag = Gas.gasflow_total_condition(state1.id, v1, state0.id, @gmodel.id)
    if flag < 0 then raise "failed to compute total condition." end
    return nil
  end

  def pitot_condition(state1, v1, state2pitot)
    flag = Gas.gasflow_pitot_condition(state1.id, v1, state2pitot.id, @gmodel.id)
    if flag < 0 then raise "failed to compute pitot condition." end
    return nil
  end

  def steady_flow_with_area_change(state1, v1, area2_over_area1, state2,
                                   tol=1.0e-4, p2p1_min=0.0001)
    my_results = [0.0].pack("d")
    flag = Gas.gasflow_steady_flow_with_area_change(state1.id, v1, area2_over_area1,
                                                    state2.id, @gmodel.id, tol,
                                                    p2p1_min, my_results)
    if flag < 0 then raise "failed to compute steady flow with area change." end
    return my_results[0, my_results.size].unpack("d")[0] # v2
  end

  def finite_wave_dp(state1, v1, characteristic, p2, state2, steps=100)
    my_results = [0.0].pack("d")
    flag = Gas.gasflow_finite_wave_dp(state1.id, v1, characteristic, p2,
                                      state2.id, @gmodel.id, steps,
                                      my_results)
    if flag < 0 then raise "failed to compute (unsteady) finite wave dp." end
    return my_results[0, my_results.size].unpack("d")[0] # v2
  end

  def finite_wave_dv(state1, v1, characteristic, v2_target, state2, steps=100,
                     t_min=200.0)
    my_results = [0.0].pack("d")
    flag = Gas.gasflow_finite_wave_dv(state1.id, v1, characteristic, v2_target,
                                      state2.id, @gmodel.id, steps, t_min,
                                      my_results)
    if flag < 0 then raise "failed to compute (unsteady) finite wave dv." end
    return my_results[0, my_results.size].unpack("d")[0] # v2
  end

  def osher_riemann(stateL, stateR, velL, velR, stateLstar, stateRstar, stateX0)
    my_results = [0.0, 0.0, 0.0, 0.0, 0.0].pack("d*")
    flag = Gas.gasflow_osher_riemann(stateL.id, stateR.id, velL, velR,
                                     stateLstar.id, stateRstar.id,
                                     stateX0.id, @gmodel.id, my_results)
    if flag < 0 then raise "failed to compute solution to Osher Riemann problem." end
    return my_results[0, my_results.size].unpack("ddddd") # [pstar, wstar, wL, wR, velX0]
  end

  def osher_flux(stateL, stateR, velL, velR)
    my_results = [0.0, 0.0, 0.0].pack("d*")
    flag = Gas.gasflow_osher_flux(stateL.id, stateR.id, velL, velR, @gmodel.id, my_results)
    if flag < 0 then raise "failed to compute Osher flux." end
    return my_results[0, my_results.size].unpack("ddd") # [f_mass, f_x_momentum, f_energy]
  end

  def roe_flux(stateL, stateR, velL, velR)
    my_results = [0.0, 0.0, 0.0].pack("d*")
    flag = Gas.gasflow_roe_flux(stateL.id, stateR.id, velL, velR, @gmodel.id, my_results)
    if flag < 0 then raise "failed to compute Roe flux." end
    return my_results[0, my_results.size].unpack("ddd") # [f_mass, f_x_momentum, f_energy]
  end

  def lrivp(stateL, stateR, velL, velR)
    my_pstar = [0.0].pack("d")
    my_wstar = [0.0].pack("d")
    flag = Gas.gasflow_lrivp(stateL.id, stateR.id, velL, velR,
                             stateL.gmodel.id, stateR.gmodel.id,
                             my_wstar, my_pstar)
    if flag < 0 then raise "failed to compute solution for lrivp." end
    pstar = my_pstar[0, my_pstar.size].unpack("d")[0]
    wstar = my_wstar[0, my_wstar.size].unpack("d")[0]
    return [pstar, wstar]
  end

  def piston_at_left(stateR, velR, wstar)
    my_pstar = [0.0].pack("d")
    flag = Gas.gasflow_piston_at_left(stateR.id, velR, stateR.gmodel.id,
                                      wstar, my_pstar)
    if flag < 0 then raise "failed to compute solution for piston_at_left." end
    return my_pstar[0, my_pstar.size].unpack("d")[0]
  end

  def piston_at_right(stateL, velL, wstar)
    my_pstar = [0.0].pack("d")
    flag = Gas.gasflow_piston_at_right(stateL.id, velL, stateL.gmodel.id,
                                      wstar, my_pstar)
    if flag < 0 then raise "failed to compute solution for piston_at_right." end
    return my_pstar[0, my_pstar.size].unpack("d")[0]
  end

  def theta_oblique(state1, v1, beta, state2)
    my_results = [0.0, 0.0].pack("dd")
    flag = Gas.gasflow_theta_oblique(state1.id, v1, beta, state2.id, @gmodel.id,
                                     my_results)
    if flag < 0 then raise "failed to compute theta oblique." end
    return my_results[0, my_results.size].unpack("dd") # [theta, v2]
  end

  def beta_oblique(state1, v1, theta)
    my_results = [0.0].pack("d")
    flag = Gas.gasflow_beta_oblique(state1.id, v1, theta, @gmodel.id,
                                    my_results)
    if flag < 0 then raise "failed to compute beta oblique." end
    return my_results[0, my_results.size].unpack("d")[0] # beta
  end

  def theta_cone(state1, v1, beta, state_c, dtheta=-0.01*Math::PI/180.0)
    my_results = [0.0, 0.0].pack("dd")
    flag = Gas.gasflow_theta_cone(state1.id, v1, beta, state_c.id, @gmodel.id,
                                  dtheta, my_results)
    if flag < 0 then raise "failed to compute theta cone." end
    return my_results[0, my_results.size].unpack("dd") # [theta_c, v_c]
  end

  def beta_cone(state1, v1, theta_c, dtheta=-0.01*Math::PI/180.0)
    my_results = [0.0].pack("d")
    flag = Gas.gasflow_beta_cone(state1.id, v1, theta_c, @gmodel.id,
                                 dtheta, my_results)
    if flag < 0 then raise "failed to compute beta cone." end
    return my_results[0, my_results.size].unpack("d")[0] # beta
  end
end
