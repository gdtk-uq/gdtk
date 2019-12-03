# gas.py
# A Python wrapper for GasModel and GasState classes.
# PJ 2019-07-24: start of experiment with FFI.
#    2019-07-25: added Python wrapper
#
from cffi import FFI

ffi = FFI()
ffi.cdef("""
    int cwrap_gas_init();

    int gas_model_new(char* file_name);
    int gas_model_n_species(int gm_i);
    int gas_model_n_modes(int gm_i);
    int gas_model_species_name(int gm_i, int isp, char* name, int n);
    int gas_model_mol_masses(int gm_i, double* mm);

    int gas_state_new(int gm_i);
    int gas_state_set_scalar_field(int gs_i, char* field_name, double value);
    int gas_state_get_scalar_field(int gs_i, char* field_name, double* value);
    int gas_state_set_array_field(int gs_i, char* field_name, double* values, int n);
    int gas_state_get_array_field(int gs_i, char* field_name, double* values, int n);
    int gas_state_copy_values(int gs_to_i, int gs_from_i);

    int gas_model_gas_state_update_thermo_from_pT(int gm_i, int gs_i);
    int gas_model_gas_state_update_thermo_from_rhou(int gm_i, int gs_i);
    int gas_model_gas_state_update_thermo_from_rhoT(int gm_i, int gs_i);
    int gas_model_gas_state_update_thermo_from_rhop(int gm_i, int gs_i);
    int gas_model_gas_state_update_thermo_from_ps(int gm_i, int gs_i, double s);
    int gas_model_gas_state_update_thermo_from_hs(int gm_i, int gs_i, double h, double s);
    int gas_model_gas_state_update_sound_speed(int gm_i, int gs_i);
    int gas_model_gas_state_update_trans_coeffs(int gm_i, int gs_i);

    int gas_model_gas_state_Cv(int gm_i, int gs_i, double* result);
    int gas_model_gas_state_Cp(int gm_i, int gs_i, double* result);
    int gas_model_gas_state_dpdrho_const_T(int gm_i, int gs_i, double* result);
    int gas_model_gas_state_R(int gm_i, int gs_i, double* result);
    int gas_model_gas_state_internal_energy(int gm_i, int gs_i, double* result);
    int gas_model_gas_state_enthalpy(int gm_i, int gs_i, double* result);
    int gas_model_gas_state_entropy(int gm_i, int gs_i, double* result);
    int gas_model_gas_state_molecular_mass(int gm_i, int gs_i, double* result);

    int gas_model_gas_state_enthalpy_isp(int gm_i, int gs_i, int isp, double* result);
    int gas_model_gas_state_entropy_isp(int gm_i, int gs_i, int isp, double* result);
    int gas_model_gas_state_gibbs_free_energy_isp(int gm_i, int gs_i, int isp, double* result);

    int gas_model_massf2molef(int gm_i, double* massf, double* molef);
    int gas_model_molef2massf(int gm_i, double* molef, double* massf);
    int gas_model_gas_state_get_molef(int gm_i, int gs_i, double* molef);
    int gas_model_gas_state_get_conc(int gm_i, int gs_i, double* conc);

    int thermochemical_reactor_new(int gm_i, char* filename1, char* filename2);
    int thermochemical_reactor_gas_state_update(int cr_i, int gs_i, double t_interval,
                                                double* dt_suggest);

    int gasflow_shock_ideal(int state1_id, double vs, int state2_id, int gm_id,
                            double* results);
    int gasflow_normal_shock(int state1_id, double vs, int state2_id, int gm_id,
                             double* results, double rho_tol, double T_tol);
    int gasflow_normal_shock_p2p1(int state1_id, double p2p1, int state2_id, int gm_id,
                                  double* results);
    int gasflow_reflected_shock(int state2_id, double vg, int state5_id, int gm_id,
                                double* results);

    int gasflow_expand_from_stagnation(int state0_id, double p_over_p0, int state1_id,
                                       int gm_id, double* results);
    int gasflow_expand_to_mach(int state0_id, double mach, int state1_id,
                               int gm_id, double* results);
    int gasflow_total_condition(int state1_id, double v1, int state0_id, int gm_id);
    int gasflow_pitot_condition(int state1_id, double v1, int state2pitot_id, int gm_id);
    int gasflow_steady_flow_with_area_change(int state1_id, double v1, double a2_over_a1,
                                             int state2_id, int gm_id, double tol,
                                             double* results);

    int gasflow_finite_wave_dp(int state1_id, double v1, char* characteristic, double p2,
                               int state2_id, int gm_id, int steps, double* results);
    int gasflow_finite_wave_dv(int state1_id, double v1, char* characteristic, double v2_target,
                               int state2_id, int gm_id, int steps, double Tmin, double* results);

    int gasflow_theta_oblique(int state1_id, double v1, double beta,
                              int state2_id, int gm_id, double* results);
    int gasflow_beta_oblique(int state1_id, double v1, double theta,
                             int gm_id, double* results);

    int gasflow_theta_cone(int state1_id, double v1, double beta,
                           int state_c_id, int gm_id, double* results);
    int gasflow_beta_cone(int state1_id, double v1, double theta,
                          int gm_id, double* results);
""")
so = ffi.dlopen("libgas.so")
so.cwrap_gas_init()

# Service classes that wrap the C-API in a nice Pythonic API...

class GasModel(object):
    def __init__(self, file_name):
        self.file_name = file_name
        self.id = so.gas_model_new(bytes(self.file_name, 'utf-8'))
        self.species_names = []
        buf = ffi.new("char[]", b'\000'*32)
        for i in range(self.n_species):
            so.gas_model_species_name(self.id, i, buf, 32)
            self.species_names.append(ffi.string(buf).decode('utf-8'))
        return

    def __str__(self):
        text = 'GasModel(file="%s", id=%d, species=%s)' % \
            (self.file_name, self.id, self.species_names)
        return text
    
    @property
    def n_species(self):
        return so.gas_model_n_species(self.id)
    
    @property
    def n_modes(self):
        return so.gas_model_n_modes(self.id)

    @property
    def mol_masses(self):
        mm = ffi.new("double[]", [0.0]*self.n_species)
        so.gas_model_mol_masses(self.id, mm)
        return [mm[i] for i in range(self.n_species)]
        
    def update_thermo_from_pT(self, gstate):
        flag = so.gas_model_gas_state_update_thermo_from_pT(self.id, gstate.id)
        if flag < 0: raise Exception("could not update thermo from p,T.")
        return
    def update_thermo_from_rhou(self, gstate):
        flag = so.gas_model_gas_state_update_thermo_from_rhou(self.id, gstate.id)
        if flag < 0: raise Exception("could not update thermo from rho,u.")
        return
    def update_thermo_from_rhoT(self, gstate):
        flag = so.gas_model_gas_state_update_thermo_from_rhoT(self.id, gstate.id)
        if flag < 0: raise Exception("could not update thermo from rho,T.")
        return
    def update_thermo_from_rhop(self, gstate):
        flag = so.gas_model_gas_state_update_thermo_from_rhop(self.id, gstate.id)
        if flag < 0: raise Exception("could not update thermo from rho,p.")
        return
    def update_thermo_from_ps(self, gstate, s):
        flag = so.gas_model_gas_state_update_thermo_from_ps(self.id, gstate.id, s)
        if flag < 0: raise Exception("could not update thermo from p,s.")
        return
    def update_thermo_from_hs(self, gstate, h, s):
        flag = so.gas_model_gas_state_update_thermo_from_hs(self.id, gstate.id, h, s)
        if flag < 0: raise Exception("could not update thermo from h,s.")
        return
    def update_sound_speed(self, gstate):
        flag = so.gas_model_gas_state_update_sound_speed(self.id, gstate.id)
        if flag < 0: raise Exception("could not update sound speed.")
        return
    def update_trans_coeffs(self, gstate):
        flag = so.gas_model_gas_state_update_trans_coeffs(self.id, gstate.id)
        if flag < 0: raise Exception("could not update transport coefficients.")
        return

    def Cv(self, gstate):
        valuep = ffi.new("double *")
        flag = so.gas_model_gas_state_Cv(self.id, gstate.id, valuep)
        if flag < 0: raise Exception("could not compute Cv.")
        return valuep[0]
    def Cp(self, gstate):
        valuep = ffi.new("double *")
        flag = so.gas_model_gas_state_Cp(self.id, gstate.id, valuep)
        if flag < 0: raise Exception("could not compute Cp.")
        return valuep[0]
    def dpdrho_const_T(self, gstate):
        valuep = ffi.new("double *")
        flag = so.gas_model_gas_state_dpdrho_const_T(self.id, gstate.id, valuep)
        if flag < 0: raise Exception("could not compute dpdrho_const_T.")
        return valuep[0]
    def R(self, gstate):
        valuep = ffi.new("double *")
        flag = so.gas_model_gas_state_R(self.id, gstate.id, valuep)
        if flag < 0: raise Exception("could not compute R.")
        return valuep[0]
    def internal_energy(self, gstate):
        valuep = ffi.new("double *")
        flag = so.gas_model_gas_state_internal_energy(self.id, gstate.id, valuep)
        if flag < 0: raise Exception("could not compute internal energy.")
        return valuep[0]
    def enthalpy(self, gstate):
        valuep = ffi.new("double *")
        flag = so.gas_model_gas_state_enthalpy(self.id, gstate.id, valuep)
        if flag < 0: raise Exception("could not compute enthalpy.")
        return valuep[0]
    def entropy(self, gstate):
        valuep = ffi.new("double *")
        flag = so.gas_model_gas_state_entropy(self.id, gstate.id, valuep)
        if flag < 0: raise Exception("could not compute entropy.")
        return valuep[0]
    def molecular_mass(self, gstate):
        valuep = ffi.new("double *")
        flag = so.gas_model_gas_state_molecular_mass(self.id, gstate.id, valuep)
        if flag < 0: raise Exception("could not compute molecular mass.")
        return valuep[0]

    def enthalpy_isp(self, gstate, isp):
        valuep = ffi.new("double *")
        flag = so.gas_model_gas_state_enthalpy_isp(self.id, gstate.id, isp, valuep)
        if flag < 0: raise Exception("could not compute enthalpy for species.")
        return valuep[0]
    def entropy_isp(self, gstate, isp):
        valuep = ffi.new("double *")
        flag = so.gas_model_gas_state_entropy_isp(self.id, gstate.id, isp, valuep)
        if flag < 0: raise Exception("could not compute entropy for species.")
        return valuep[0]
    def gibbs_free_energy_isp(self, gstate, isp):
        valuep = ffi.new("double *")
        flag = so.gas_model_gas_state_gibbs_free_energy_isp(self.id, gstate.id, isp, valuep)
        if flag < 0: raise Exception("could not compute gibbs free energy for species.")
        return valuep[0]

    def massf2molef(self, massf_given):
        nsp = self.n_species
        if type(massf_given) == type([]):
            massf_list = massf_given.copy()
            assert len(massf_list) == self.n_species, "incorrect massf list length"
        elif type(massf_given) == type({}):
            massf_list = []
            for name in self.species_names:
                if name in massf_given.keys():
                    massf_list.append(massf_given[name])
                else:
                    massf_list.append(0.0)
        if abs(sum(massf_list) - 1.0) > 1.0e-6:
            raise Exception("mass fractions do not sum to 1.")
        my_massf = ffi.new("double[]", massf_list)
        my_molef = ffi.new("double[]", [0.0]*self.n_species)
        so.gas_model_massf2molef(self.id, my_massf, my_molef)
        return [my_molef[i] for i in range(self.n_species)]

    def molef2massf(self, molef_given):
        nsp = self.n_species
        if type(molef_given) == type([]):
            molef_list = molef_given.copy()
            assert len(molef_list) == self.n_species, "incorrect molef list length"
        elif type(molef_given) == type({}):
            molef_list = []
            for name in self.species_names:
                if name in molef_given.keys():
                    molef_list.append(molef_given[name])
                else:
                    molef_list.append(0.0)
        if abs(sum(molef_list) - 1.0) > 1.0e-6:
            raise Exception("mole fractions do not sum to 1.")
        my_molef = ffi.new("double[]", molef_list)
        my_massf = ffi.new("double[]", [0.0]*self.n_species)
        so.gas_model_molef2massf(self.id, my_molef, my_massf)
        return [my_massf[i] for i in range(self.n_species)]

    
class GasState(object):
    def __init__(self, gmodel):
        self.gmodel = gmodel
        self.id = so.gas_state_new(self.gmodel.id)

    def __str__(self):
        text = 'GasState(rho=%g' % self.rho
        text += ', p=%g' % self.p
        text += ', T=%g' % self.T
        text += ', u=%g' % self.u
        text += ', a=%g' % self.a
        text += ', massf=%s' % str(self.massf)
        text += ', id=%d, gmodel.id=%d)' % (self.id, self.gmodel.id)
        return text
    
    @property
    def rho(self):
        valuep = ffi.new("double *")
        flag = so.gas_state_get_scalar_field(self.id, b"rho", valuep)
        if flag < 0: raise Exception("could not get density.")
        return valuep[0]
    @rho.setter
    def rho(self, value):
        flag = so.gas_state_set_scalar_field(self.id, b"rho", value)
        if flag < 0: raise Exception("could not set density.")
        return
    
    @property
    def p(self):
        valuep = ffi.new("double *")
        flag = so.gas_state_get_scalar_field(self.id, b"p", valuep)
        if flag < 0: raise Exception("could not get pressure.")
        return valuep[0]
    @p.setter
    def p(self, value):
        flag = so.gas_state_set_scalar_field(self.id, b"p", value)
        if flag < 0: raise Exception("could not set pressure.")
        return
    
    @property
    def T(self):
        valuep = ffi.new("double *")
        flag = so.gas_state_get_scalar_field(self.id, b"T", valuep)
        if flag < 0: raise Exception("could not get temperature.")
        return valuep[0]
    @T.setter
    def T(self, value):
        flag = so.gas_state_set_scalar_field(self.id, b"T", value)
        if flag < 0: raise Exception("could not set temperature.")
        return
    
    @property
    def u(self):
        valuep = ffi.new("double *")
        flag = so.gas_state_get_scalar_field(self.id, b"u", valuep)
        if flag < 0: raise Exception("could not get internal energy.")
        return valuep[0]
    @u.setter
    def u(self, value):
        flag = so.gas_state_set_scalar_field(self.id, b"u", value)
        if flag < 0: raise Exception("could not set internal energy.")
        return
    
    @property
    def a(self):
        valuep = ffi.new("double *")
        flag = so.gas_state_get_scalar_field(self.id, b"a", valuep)
        if flag < 0: raise Exception("could not get sound speed.")
        return valuep[0]
    
    @property
    def k(self):
        valuep = ffi.new("double *")
        flag = so.gas_state_get_scalar_field(self.id, b"k", valuep)
        if flag < 0: raise Exception("could not get conductivity.")
        return valuep[0]
    
    @property
    def mu(self):
        valuep = ffi.new("double *")
        flag = so.gas_state_get_scalar_field(self.id, b"mu", valuep)
        if flag < 0: raise Exception("could not get viscosity.")
        return valuep[0]

    @property
    def massf(self):
        nsp = self.gmodel.n_species
        mf = ffi.new("double[]", [0.0]*nsp)
        flag = so.gas_state_get_array_field(self.id, b"massf", mf, nsp)
        if flag < 0: raise Exception("could not get mass-fractions.")
        return [mf[i] for i in range(nsp)]
    @massf.setter
    def massf(self, mf_given):
        nsp = self.gmodel.n_species
        if type(mf_given) == type([]):
            mf_list = mf_given.copy()
        elif type(mf_given) == type({}):
            mf_list = []
            for name in self.gmodel.species_names:
                if name in mf_given.keys():
                    mf_list.append(mf_given[name])
                else:
                    mf_list.append(0.0)
        if abs(sum(mf_list) - 1.0) > 1.0e-6:
            raise Exception("mass fractions do not sum to 1.")
        mf = ffi.new("double[]", mf_list)
        flag = so.gas_state_set_array_field(self.id, b"massf", mf, nsp)
        if flag < 0: raise Exception("could not set mass-fractions.")
        return mf_list

    @property
    def molef(self):
        nsp = self.gmodel.n_species
        mf = ffi.new("double[]", [0.0]*nsp)
        flag = so.gas_model_gas_state_get_molef(self.gmodel.id, self.id, mf)
        if flag < 0: raise Exception("could not get mole-fractions.")
        return [mf[i] for i in range(nsp)]
    @molef.setter
    def molef(self, molef_given):
        nsp = self.gmodel.n_species
        mf_list = self.gmodel.molef2massf(molef_given)
        mf = ffi.new("double[]", mf_list)
        flag = so.gas_state_set_array_field(self.id, b"massf", mf, nsp)
        if flag < 0: raise Exception("could not set mass-fractions from mole-fractions.")
        # At this point, we may not have the mole-fractions as a list.
        return None

    @property
    def conc(self):
        nsp = self.gmodel.n_species
        myconc = ffi.new("double[]", [0.0]*nsp)
        flag = so.gas_model_gas_state_get_conc(self.gmodel.id, self.id, myconc)
        if flag < 0: raise Exception("could not get concentrations.")
        return [myconc[i] for i in range(nsp)]

    @property
    def u_modes(self):
        n = self.gmodel.n_modes
        if n == 0: return []
        um = ffi.new("double[]", [0.0]*n)
        flag = so.gas_state_get_array_field(self.id, b"u_modes", um, nsp)
        if flag < 0: raise Exception("could not get u_modes.")
        return [um[i] for i in range(n)]
    @u_modes.setter
    def u_modes(self, um_given):
        n = self.gmodel.n_modes
        if n == 0: return []
        if type(um_given) != type([]):
            raise Exception("u_modes needs to be supplied as a list.")
        um = ffi.new("double[]", um_given)
        flag = so.gas_state_set_array_field(self.id, b"u_modes", um, n)
        if flag < 0: raise Exception("could not set u_modes.")
        return um_given

    @property
    def T_modes(self):
        n = self.gmodel.n_modes
        if n == 0: return []
        Tm = ffi.new("double[]", [0.0]*n)
        flag = so.gas_state_get_array_field(self.id, b"T_modes", Tm, nsp)
        if flag < 0: raise Exception("could not get T_modes.")
        return [um[i] for i in range(n)]
    @T_modes.setter
    def T_modes(self, Tm_given):
        n = self.gmodel.n_modes
        if n == 0: return []
        if type(Tm_given) != type([]):
            raise Exception("T_modes needs to be supplied as a list.")
        Tm = ffi.new("double[]", Tm_given)
        flag = so.gas_state_set_array_field(self.id, b"T_modes", Tm, n)
        if flag < 0: raise Exception("could not set T_modes.")
        return Tm_given

    @property
    def k_modes(self):
        n = self.gmodel.k_modes
        if n == 0: return []
        km = ffi.new("double[]", [0.0]*n)
        flag = so.gas_state_get_array_field(self.id, b"k_modes", km, nsp)
        if flag < 0: raise Exception("could not get k_modes.")
        return [km[i] for i in range(n)]

    def copy_values(self, gstate):
        flag = so.gas_state_copy_values(self.id, gstate.id)
        if flag < 0: raise Exception("could not copy values.")
        return
            
    def update_thermo_from_pT(self):
        self.gmodel.update_thermo_from_pT(self)
        return
    def update_thermo_from_rhou(self):
        self.gmodel.update_thermo_from_rhou(self)
        return
    def update_thermo_from_rhoT(self):
        self.gmodel.update_thermo_from_rhoT(self)
        return
    def update_thermo_from_rhop(self):
        self.gmodel.update_thermo_from_rhop(self)
        return
    def update_thermo_from_ps(self, s):
        self.gmodel.update_thermo_from_ps(self, s)
        return
    def update_thermo_from_hs(self, h, s):
        self.gmodel.update_thermo_from_hs(self, h, s)
        return
    def update_sound_speed(self):
        self.gmodel.update_sound_speed(self)
        return
    def update_trans_coeffs(self):
        self.gmodel.update_trans_coeffs(self)
        return

    @property
    def Cv(self):
        return self.gmodel.Cv(self)
    @property
    def Cp(self):
        return self.gmodel.Cp(self)
    @property
    def dpdrho_const_T(self):
        return self.gmodel.dpdrho_const_T(self)
    @property
    def R(self):
        return self.gmodel.R(self)
    @property
    def internal_energy(self):
        return self.gmodel.internal_energy(self)
    @property
    def enthalpy(self):
        return self.gmodel.enthalpy(self)
    @property
    def entropy(self):
        return self.gmodel.entropy(self)
    @property
    def molecular_mass(self):
        return self.gmodel.molecular_mass(self)

    def enthalpy_isp(self, isp):
        return self.gmodel.enthalpy_isp(self, isp)
    def entropy_isp(self, isp):
        return self.gmodel.entropy_isp(self, isp)
    def gibbs_free_energy_isp(self, isp):
        return self.gmodel.gibbs_free_energy_isp(self, isp)


class ThermochemicalReactor(object):
    def __init__(self, gmodel, filename1, filename2=""):
        self.filename1 = filename1
        self.filename2 = filename2
        self.gmodel = gmodel
        self.id = so.thermochemical_reactor_new(self.gmodel.id,
                                                bytes(self.filename1, 'utf-8'),
                                                bytes(self.filename2, 'utf-8'))

    def __str__(self):
        text = 'ThermochemicalReactor(id=%d, gmodel.id=%d, file1="%s", file2="%s")' % \
            (self.id, self.gmodel.id, self.filename1, self.filename2)
        return text

    def update_state(self, gstate, t_interval, dt_suggest):
        dt_suggestp = ffi.new("double *")
        dt_suggestp[0] = dt_suggest
        flag = so.thermochemical_reactor_gas_state_update(self.id, gstate.id,
                                                          t_interval, dt_suggestp)
        if flag < 0: raise Exception("could not update state.")
        return dt_suggestp[0]

#------------------------------------------------------------------------------------

class GasFlow(object):
    def __init__(self, gmodel):
        self.gmodel = gmodel

    def __str__(self):
        return "GasFlow(gmodel.id=%d)" % self.gmodel.id

    def ideal_shock(self, state1, vs, state2):
        my_results = ffi.new("double[]", [0.0]*2)
        flag = so.gasflow_shock_ideal(state1.id, vs, state2.id, self.gmodel.id, my_results)
        if flag < 0: raise Exception("failed to compute ideal shock jump.")
        v2 = my_results[0]
        vg = my_results[1]
        return [v2, vg]

    def normal_shock(self, state1, vs, state2, rho_tol=1.0e-6, T_tol=0.1):
        my_results = ffi.new("double[]", [0.0]*2)
        flag = so.gasflow_normal_shock(state1.id, vs, state2.id, self.gmodel.id, my_results,
                                       rho_tol, T_tol)
        if flag < 0: raise Exception("failed to compute normal shock jump.")
        v2 = my_results[0]
        vg = my_results[1]
        return [v2, vg]

    def normal_shock_p2p1(self, state1, p2p1, state2):
        my_results = ffi.new("double[]", [0.0]*3)
        flag = so.gasflow_normal_shock_p2p1(state1.id, p2p1, state2.id, self.gmodel.id,
                                            my_results)
        if flag < 0: raise Exception("failed to compute normal shock jump from p2p1.")
        vs = my_results[0]
        v2 = my_results[1]
        vg = my_results[2]
        return [vs, v2, vg]

    def reflected_shock(self, state2, vg, state5):
        my_results = ffi.new("double[]", [0.0])
        flag = so.gasflow_reflected_shock(state2.id, vg, state5.id, self.gmodel.id,
                                          my_results)
        if flag < 0: raise Exception("failed to compute reflected shock.")
        vr = my_results[0]
        return vr

    def expand_from_stagnation(self, state0, p_over_p0, state1):
        my_results = ffi.new("double[]", [0.0])
        flag = so.gasflow_expand_from_stagnation(state0.id, p_over_p0, state1.id,
                                                 self.gmodel.id, my_results)
        if flag < 0: raise Exception("failed to compute expansion from stagnation.")
        v = my_results[0]
        return v

    def expand_to_mach(self, state0, mach, state1):
        my_results = ffi.new("double[]", [0.0])
        flag = so.gasflow_expand_to_mach(state0.id, mach, state1.id,
                                         self.gmodel.id, my_results)
        if flag < 0: raise Exception("failed to compute expansion to mach number.")
        v = my_results[0]
        return v

    def total_condition(self, state1, v1, state0):
        flag = so.gasflow_total_condition(state1.id, v1, state0.id, self.gmodel.id)
        if flag < 0: raise Exception("failed to compute total condition.")
        return

    def pitot_condition(self, state1, v1, state2pitot):
        flag = so.gasflow_pitot_condition(state1.id, v1, state2pitot.id, self.gmodel.id)
        if flag < 0: raise Exception("failed to compute pitot condition.")
        return

    def steady_flow_with_area_change(self, state1, v1, area2_over_area1, state2,
                                     tol=1.0e-4):
        my_results = ffi.new("double[]", [0.0])
        flag = so.gasflow_steady_flow_with_area_change(state1.id, v1, area2_over_area1,
                                                       state2.id, self.gmodel.id, tol,
                                                       my_results)
        if flag < 0: raise Exception("failed to compute steady flow with area change.")
        v2 = my_results[0]
        return v2

    def finite_wave_dp(self, state1, v1, characteristic, p2, state2, steps=100):
        char_name = bytes(characteristic, 'utf-8')
        my_results = ffi.new("double[]", [0.0])
        flag = so.gasflow_finite_wave_dp(state1.id, v1, char_name, p2,
                                         state2.id, self.gmodel.id, steps,
                                         my_results)
        if flag < 0: raise Exception("failed to compute (unsteady) finite wave dp.")
        v2 = my_results[0]
        return v2

    def finite_wave_dv(self, state1, v1, characteristic, v2_target, state2,
                       steps=100, t_min=200.0):
        char_name = bytes(characteristic, 'utf-8')
        my_results = ffi.new("double[]", [0.0])
        flag = so.gasflow_finite_wave_dv(state1.id, v1, char_name, v2_target,
                                         state2.id, self.gmodel.id, steps, t_min,
                                         my_results)
        if flag < 0: raise Exception("failed to compute (unsteady) finite wave dv.")
        v2 = my_results[0]
        return v2

    def theta_oblique(self, state1, v1, beta, state2):
        my_results = ffi.new("double[]", [0.0, 0.0])
        flag = so.gasflow_theta_oblique(state1.id, v1, beta,
                                        state2.id, self.gmodel.id, my_results)
        if flag < 0: raise Exception("failed to compute theta oblique.")
        theta = my_results[0]
        v2 = my_results[1]
        return theta, v2

    def beta_oblique(self, state1, v1, theta):
        my_results = ffi.new("double[]", [0.0])
        flag = so.gasflow_beta_oblique(state1.id, v1, theta,
                                       self.gmodel.id, my_results)
        if flag < 0: raise Exception("failed to compute beta oblique.")
        beta = my_results[0]
        return beta

    def theta_cone(self, state1, v1, beta, state_c):
        my_results = ffi.new("double[]", [0.0, 0.0])
        flag = so.gasflow_theta_cone(state1.id, v1, beta,
                                     state_c.id, self.gmodel.id, my_results)
        if flag < 0: raise Exception("failed to compute theta cone.")
        theta_c = my_results[0]
        v2_c = my_results[1]
        return theta_c, v2_c

    def beta_cone(self, state1, v1, theta):
        my_results = ffi.new("double[]", [0.0])
        flag = so.gasflow_beta_cone(state1.id, v1, theta,
                                    self.gmodel.id, my_results)
        if flag < 0: raise Exception("failed to compute beta cone.")
        beta = my_results[0]
        return beta
