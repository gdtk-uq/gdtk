# gas.py
# A Python wrapper for GasModel, GasState and ThermochemicalReactor classes.
# It also wraps the simple 1D flow calculations such as shocks and expansions.
#
# PJ 2019-07-24: start of experiment with FFI.
#    2019-07-25: added Python wrapper
#    2023-06-01: added PyGasState with shadow attributes
#
PC_P_atm = 101.325e3

from cffi import FFI
import math

ffi = FFI()
ffi.cdef("""
    int cwrap_gas_init();

    int gas_model_new(char* file_name);
    int gas_model_type_str(int gm_i, char* dest_str, int n);
    int gas_model_n_species(int gm_i);
    int gas_model_n_modes(int gm_i);
    int gas_model_species_name(int gm_i, int isp, char* name, int n);
    int gas_model_mol_masses(int gm_i, double* mm);

    int gas_state_new(int gm_i);
    int gas_state_set_scalar_field(int gs_i, char* field_name, double value);
    int gas_state_get_scalar_field(int gs_i, char* field_name, double* value);
    int gas_state_get_thermo_scalars(int gs_i, double* values);
    int gas_state_set_array_field(int gs_i, char* field_name, double* values, int n);
    int gas_state_get_array_field(int gs_i, char* field_name, double* values, int n);
    int gas_state_get_ceaSavedData_field(int gs_i, char* field_name, double* value);
    int gas_state_get_ceaSavedData_massf(int gs_i, char* species_name, double* value);
    int gas_state_get_ceaSavedData_species_names(int gs_i, char* dest_str, int n);
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
    int gas_model_gas_state_gamma(int gm_i, int gs_i, double* result);
    int gas_model_gas_state_Prandtl(int gm_i, int gs_i, double* result);
    int gas_model_gas_state_internal_energy(int gm_i, int gs_i, double* result);
    int gas_model_gas_state_enthalpy(int gm_i, int gs_i, double* result);
    int gas_model_gas_state_entropy(int gm_i, int gs_i, double* result);
    int gas_model_gas_state_molecular_mass(int gm_i, int gs_i, double* result);
    int gas_model_gas_state_binary_diffusion_coefficients(int gm_i, int gs_i, double* dij);

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
    int thermochemical_reactor_eval_source_terms(int cr_i, int gm_i, int gs_i, int nsp, int nmodes, double* source);

    int reaction_mechanism_new(int gm_i, char* filename);
    int reaction_mechanism_n_reactions(int rm_i);
    int reaction_mechanism_tickrates(int rm_i, int gm_i, int gs_i, double* forwardrates, double* backwardrates);

    int gasflow_shock_ideal(int state1_id, double vs, int state2_id, int gm_id,
                            double* results);
    int gasflow_normal_shock(int state1_id, double vs, int state2_id, int gm_id,
                             double* results, double rho_tol, double T_tol);
    int gasflow_normal_shock_1(int state1_id, double vs, int state2_id, int gm_id,
                               double* results, double p_tol, double T_tol);
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
                                             int state2_id, int gm_id, double tol, double p2p1_min,
                                             double* results);

    int gasflow_finite_wave_dp(int state1_id, double v1, char* characteristic, double p2,
                               int state2_id, int gm_id, int steps, double* results);
    int gasflow_finite_wave_dv(int state1_id, double v1, char* characteristic, double v2_target,
                               int state2_id, int gm_id, int steps, double Tmin, double* results);

    int gasflow_osher_riemann(int stateL_id, int stateR_id, double velL, double velR,
                              int stateLstar_id, int stateRstar_id,
                              int stateX0_id, int gm_id, double* results);
    int gasflow_osher_flux(int stateL_id, int stateR_id, double velL, double velR,
                           int gm_id, double* results);
    int gasflow_roe_flux(int stateL_id, int stateR_id, double velL, double velR,
                         int gm_id, double* results);

    int gasflow_lrivp(int stateL_id, int stateR_id, double velL, double velR,
                      int gmL_id, int gmR_id, double* wstar, double* pstar);
    int gasflow_piston_at_left(int stateR_id, double velR, int gm_id,
                               double wstar, double* pstar);
    int gasflow_piston_at_right(int stateL_id, double velL, int gm_id,
                                double wstar, double* pstar);

    int gasflow_theta_oblique(int state1_id, double v1, double beta,
                              int state2_id, int gm_id, double* results);
    int gasflow_beta_oblique(int state1_id, double v1, double theta,
                             int gm_id, double* results);

    int gasflow_theta_cone(int state1_id, double v1, double beta,
                           int state_c_id, int gm_id, double dtheta, double* results);
    int gasflow_beta_cone(int state1_id, double v1, double theta,
                          int gm_id, double dtheta, double* results);
""")
so = ffi.dlopen("libgas.so")
so.cwrap_gas_init()

# -----------------------------------------------------------------------------------
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
        # Allocate ffi buffers here, so that they can be reused.
        self._valuep = ffi.new("double *")
        self._massf = ffi.new("double[]", [0.0]*self.n_species)
        self._molef = ffi.new("double[]", [0.0]*self.n_species)
        self._modes = ffi.new("double[]", [0.0]*self.n_modes)
        return

    def __copy__(self):
        gm_new = GasModel(self.file_name)
        return gm_new

    def __deepcopy__(self, memo):
        gm_new = GasModel(self.file_name)
        return gm_new

    def __str__(self):
        text = 'GasModel(file="%s", id=%d, species=%s)' % \
            (self.file_name, self.id, self.species_names)
        return text

    @property
    def type_str(self):
        buf = ffi.new("char[]", b'\000'*32)
        so.gas_model_type_str(self.id, buf, 32)
        return ffi.string(buf).decode('utf-8')

    @property
    def n_species(self):
        return so.gas_model_n_species(self.id)

    @property
    def n_modes(self):
        return so.gas_model_n_modes(self.id)

    @property
    def mol_masses(self):
        so.gas_model_mol_masses(self.id, self._massf)
        return [self._massf[i] for i in range(self.n_species)]

    def update_thermo_from_pT(self, gstate):
        flag = so.gas_model_gas_state_update_thermo_from_pT(self.id, gstate.id)
        if flag < 0: raise Exception("could not update thermo from p,T.")
        self.update_sound_speed(gstate)
        return
    def update_thermo_from_rhou(self, gstate):
        flag = so.gas_model_gas_state_update_thermo_from_rhou(self.id, gstate.id)
        if flag < 0: raise Exception("could not update thermo from rho,u.")
        self.update_sound_speed(gstate)
        return
    def update_thermo_from_rhoT(self, gstate):
        flag = so.gas_model_gas_state_update_thermo_from_rhoT(self.id, gstate.id)
        if flag < 0: raise Exception("could not update thermo from rho,T.")
        self.update_sound_speed(gstate)
        return
    def update_thermo_from_rhop(self, gstate):
        flag = so.gas_model_gas_state_update_thermo_from_rhop(self.id, gstate.id)
        if flag < 0: raise Exception("could not update thermo from rho,p.")
        self.update_sound_speed(gstate)
        return
    def update_thermo_from_ps(self, gstate, s):
        flag = so.gas_model_gas_state_update_thermo_from_ps(self.id, gstate.id, s)
        if flag < 0: raise Exception("could not update thermo from p,s.")
        self.update_sound_speed(gstate)
        return
    def update_thermo_from_hs(self, gstate, h, s):
        flag = so.gas_model_gas_state_update_thermo_from_hs(self.id, gstate.id, h, s)
        if flag < 0: raise Exception("could not update thermo from h,s.")
        self.update_sound_speed(gstate)
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
        flag = so.gas_model_gas_state_Cv(self.id, gstate.id, self._valuep)
        if flag < 0: raise Exception("could not compute Cv.")
        return self._valuep[0]
    def Cp(self, gstate):
        flag = so.gas_model_gas_state_Cp(self.id, gstate.id, self._valuep)
        if flag < 0: raise Exception("could not compute Cp.")
        return self._valuep[0]
    def dpdrho_const_T(self, gstate):
        flag = so.gas_model_gas_state_dpdrho_const_T(self.id, gstate.id, self._valuep)
        if flag < 0: raise Exception("could not compute dpdrho_const_T.")
        return self._valuep[0]
    def R(self, gstate):
        flag = so.gas_model_gas_state_R(self.id, gstate.id, self._valuep)
        if flag < 0: raise Exception("could not compute R.")
        return self._valuep[0]
    def gamma(self, gstate):
        flag = so.gas_model_gas_state_gamma(self.id, gstate.id, self._valuep)
        if flag < 0: raise Exception("could not compute gamma.")
        return self._valuep[0]
    def Prandtl(self, gstate):
        flag = so.gas_model_gas_state_Prandtl(self.id, gstate.id, self._valuep)
        if flag < 0: raise Exception("could not compute Prandtl.")
        return self._valuep[0]
    def internal_energy(self, gstate):
        flag = so.gas_model_gas_state_internal_energy(self.id, gstate.id, self._valuep)
        if flag < 0: raise Exception("could not compute internal energy.")
        return self._valuep[0]
    def enthalpy(self, gstate):
        flag = so.gas_model_gas_state_enthalpy(self.id, gstate.id, self._valuep)
        if flag < 0: raise Exception("could not compute enthalpy.")
        return self._valuep[0]
    def entropy(self, gstate):
        flag = so.gas_model_gas_state_entropy(self.id, gstate.id, self._valuep)
        if flag < 0: raise Exception("could not compute entropy.")
        return self._valuep[0]
    def molecular_mass(self, gstate):
        flag = so.gas_model_gas_state_molecular_mass(self.id, gstate.id, self._valuep)
        if flag < 0: raise Exception("could not compute molecular mass.")
        return self._valuep[0]
    def binary_diffusion_coefficients(self, gstate):
        nsp = self.n_species
        Dij = ffi.new("double[]", [0.0]*nsp*nsp)
        flag = so.gas_model_gas_state_binary_diffusion_coefficients(self.id, gstate.id, Dij)
        if flag < 0: raise Exception("could not compute binary diffusion coefficients.")
        return [[Dij[i*self.n_species + j] for j in range(self.n_species)] for i in range(self.n_species)]


    def enthalpy_isp(self, gstate, isp):
        flag = so.gas_model_gas_state_enthalpy_isp(self.id, gstate.id, isp, self._valuep)
        if flag < 0: raise Exception("could not compute enthalpy for species.")
        return self._valuep[0]
    def entropy_isp(self, gstate, isp):
        flag = so.gas_model_gas_state_entropy_isp(self.id, gstate.id, isp, self._valuep)
        if flag < 0: raise Exception("could not compute entropy for species.")
        return self._valuep[0]
    def gibbs_free_energy_isp(self, gstate, isp):
        flag = so.gas_model_gas_state_gibbs_free_energy_isp(self.id, gstate.id, isp, self._valuep)
        if flag < 0: raise Exception("could not compute gibbs free energy for species.")
        return self._valuep[0]

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
        for i in range(self.nspecies): self._massf[i] = massf_list[i]
        so.gas_model_massf2molef(self.id, self._massf, self._molef)
        return [self._molef[i] for i in range(self.n_species)]

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
        for i in range(self.n_species): self._molef[i] = molef_list[i]
        so.gas_model_molef2massf(self.id, self._molef, self._massf)
        return [self._massf[i] for i in range(self.n_species)]


# -----------------------------------------------------------------------------------
class GasState(object):
    """
    A GasState that has most of its data in the Dlang domain.

    This is the basic wrapper class but you might like to consider
    using objects of the PyGasState class if you want to access
    the data attributes a lot from the Python domain.
    """
    def __init__(self, gmodel):
        self.gmodel = gmodel
        self.id = so.gas_state_new(self.gmodel.id)
        # Allocate ffi buffers here, so that they can be reused.
        self._valuep = ffi.new("double *")
        self._mf = ffi.new("double[]", [0.0]*gmodel.n_species)
        self._modes = ffi.new("double[]", [0.0]*gmodel.n_modes)
        return

    def __copy__(self):
        gs_new = GasState(self.gmodel)
        gs_new.p = self.p
        gs_new.T = self.T
        if self.gmodel.n_modes > 0:
            gs_new.T_modes = self.T_modes
        gs_new.massf = self.massf
        if (self.p > 0.0) and (self.T > 0.0):
            gs_new.update_thermo_from_pT()
        return gs_new

    def __deepcopy__(self, memo):
        gs_new = GasState(self.gmodel)
        gs_new.p = self.p
        gs_new.T = self.T
        if self.gmodel.n_modes > 0:
            gs_new.T_modes = self.T_modes
        gs_new.massf = self.massf
        if (self.p > 0.0) and (self.T > 0.0):
            gs_new.update_thermo_from_pT()
        return gs_new

    def __str__(self):
        text = 'GasState(rho=%g' % self.rho
        text += ', p=%g' % self.p
        text += ', T=%g' % self.T
        text += ', u=%g' % self.u
        if self.gmodel.n_modes > 0:
            text += ', T_modes=%s' % self.T_modes
            text += ', u_modes=%s' % self.u_modes
        text += ', a=%g' % self.a
        if self.gmodel.n_species > 1:
            text += ', massf=%s' % str(self.massf)
        text += ', id=%d, gmodel.id=%d)' % (self.id, self.gmodel.id)
        return text

    @property
    def rho(self):
        flag = so.gas_state_get_scalar_field(self.id, b"rho", self._valuep)
        if flag < 0: raise Exception("could not get density.")
        return self._valuep[0]
    @rho.setter
    def rho(self, value):
        flag = so.gas_state_set_scalar_field(self.id, b"rho", value)
        if flag < 0: raise Exception("could not set density.")
        return

    @property
    def p(self):
        flag = so.gas_state_get_scalar_field(self.id, b"p", self._valuep)
        if flag < 0: raise Exception("could not get pressure.")
        return self._valuep[0]
    @p.setter
    def p(self, value):
        flag = so.gas_state_set_scalar_field(self.id, b"p", value)
        if flag < 0: raise Exception("could not set pressure.")
        return

    @property
    def T(self):
        flag = so.gas_state_get_scalar_field(self.id, b"T", self._valuep)
        if flag < 0: raise Exception("could not get temperature.")
        return self._valuep[0]
    @T.setter
    def T(self, value):
        flag = so.gas_state_set_scalar_field(self.id, b"T", value)
        if flag < 0: raise Exception("could not set temperature.")
        return

    @property
    def u(self):
        flag = so.gas_state_get_scalar_field(self.id, b"u", self._valuep)
        if flag < 0: raise Exception("could not get internal energy.")
        return self._valuep[0]
    @u.setter
    def u(self, value):
        flag = so.gas_state_set_scalar_field(self.id, b"u", value)
        if flag < 0: raise Exception("could not set internal energy.")
        return

    @property
    def a(self):
        flag = so.gas_state_get_scalar_field(self.id, b"a", self._valuep)
        if flag < 0: raise Exception("could not get sound speed.")
        return self._valuep[0]

    @property
    def k(self):
        flag = so.gas_state_get_scalar_field(self.id, b"k", self._valuep)
        if flag < 0: raise Exception("could not get conductivity.")
        return self._valuep[0]

    @property
    def mu(self):
        flag = so.gas_state_get_scalar_field(self.id, b"mu", self._valuep)
        if flag < 0: raise Exception("could not get viscosity.")
        return self._valuep[0]

    @property
    def n_species(self):
        return self.gmodel.n_species

    @property
    def massf(self):
        """
        Mass fractions are always returned as a list.
        """
        nsp = self.gmodel.n_species
        flag = so.gas_state_get_array_field(self.id, b"massf", self._mf, nsp)
        if flag < 0: raise Exception("could not get mass-fractions.")
        return [self._mf[i] for i in range(nsp)]
    @property
    def massf_as_dict(self):
        nsp = self.gmodel.n_species
        names = self.gmodel.species_names
        mf = self.massf
        result = {}
        for i in range(nsp): result[names[i]] = mf[i]
        return result
    @massf.setter
    def massf(self, mf_given):
        """
        Mass fractions may be provided as a list or dictionary.
        """
        nsp = self.gmodel.n_species
        if type(mf_given) == type([]):
            if len(mf_given) != nsp:
                raise Exception(f"mass fraction list is not correct length. nsp={nsp}; len(massf)={len(mf_given)}")
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
        for i in range(nsp): self._mf[i] = mf_list[i]
        flag = so.gas_state_set_array_field(self.id, b"massf", self._mf, nsp)
        if flag < 0: raise Exception("could not set mass-fractions.")
        return mf_list

    @property
    def molef(self):
        """
        Mole fractions are always returned as a list.
        """
        nsp = self.gmodel.n_species
        flag = so.gas_model_gas_state_get_molef(self.gmodel.id, self.id, self._mf)
        if flag < 0: raise Exception("could not get mole-fractions.")
        return [self._mf[i] for i in range(nsp)]
    @property
    def molef_as_dict(self):
        nsp = self.gmodel.n_species
        names = self.gmodel.species_names
        mf = self.molef
        result = {}
        for i in range(nsp): result[names[i]] = mf[i]
        return result
    @molef.setter
    def molef(self, molef_given):
        """
        Mole fractions may be provided as a list or dictionary.
        """
        nsp = self.gmodel.n_species
        mf_list = self.gmodel.molef2massf(molef_given)
        for i in range(nsp): self._mf[i] = mf_list[i]
        flag = so.gas_state_set_array_field(self.id, b"massf", self._mf, nsp)
        if flag < 0: raise Exception("could not set mass-fractions from mole-fractions.")
        # At this point, we may not have the mole-fractions as a list
        # because it may have been provided as a dictionary.
        # So, don't return anything.
        return None

    @property
    def conc(self):
        nsp = self.gmodel.n_species
        flag = so.gas_model_gas_state_get_conc(self.gmodel.id, self.id, self._mf)
        if flag < 0: raise Exception("could not get concentrations.")
        return [self._mf[i] for i in range(nsp)]
    @property
    def conc_as_dict(self):
        nsp = self.gmodel.n_species
        names = self.gmodel.species_names
        conc_list = self.conc
        result = {}
        for i in range(nsp): result[names[i]] = conc_list[i]
        return result

    @property
    def n_modes(self):
        return self.gmodel.n_modes

    @property
    def u_modes(self):
        n = self.gmodel.n_modes
        if n == 0: return []
        flag = so.gas_state_get_array_field(self.id, b"u_modes", self._modes, n)
        if flag < 0: raise Exception("could not get u_modes.")
        return [self._modes[i] for i in range(n)]
    @u_modes.setter
    def u_modes(self, um_given):
        n = self.gmodel.n_modes
        if n == 0: return []
        if type(um_given) != type([]):
            raise Exception("u_modes needs to be supplied as a list.")
        if len(um_given) != n:
            raise Exception(f"u_modes list is not correct length. nmodes={n}; len(u_modes)={len(um_given)}")
        for i in range(n): self._modes[i] = um_given[i]
        flag = so.gas_state_set_array_field(self.id, b"u_modes", self._modes, n)
        if flag < 0: raise Exception("could not set u_modes.")
        return um_given

    @property
    def T_modes(self):
        n = self.gmodel.n_modes
        if n == 0: return []
        flag = so.gas_state_get_array_field(self.id, b"T_modes", self._modes, n)
        if flag < 0: raise Exception("could not get T_modes.")
        return [self._modes[i] for i in range(n)]
    @T_modes.setter
    def T_modes(self, Tm_given):
        n = self.gmodel.n_modes
        if n == 0: return []
        if type(Tm_given) != type([]):
            raise Exception("T_modes needs to be supplied as a list.")
        if len(Tm_given) != n:
            raise Exception(f"T_modes list is not correct length. nmodes={n}; len(T_modes)={len(Tm_given)}")
        for i in range(n): self._modes[i] = Tm_given[i]
        flag = so.gas_state_set_array_field(self.id, b"T_modes", self._modes, n)
        if flag < 0: raise Exception("could not set T_modes.")
        return Tm_given

    @property
    def k_modes(self):
        n = self.gmodel.n_modes
        if n == 0: return []
        flag = so.gas_state_get_array_field(self.id, b"k_modes", self._modes, n)
        if flag < 0: raise Exception("could not get k_modes.")
        return [self._modes[i] for i in range(n)]

    @property
    def ceaSavedData(self):
        my_data = {}
        scalar_fields = ["p", "rho", "u", "h", "T", "a",
                         "Mmass", "Rgas", "gamma", "Cp",
                         "s", "k", "mu"]
        valuep = ffi.new("double *")
        for name in scalar_fields:
            flag = so.gas_state_get_ceaSavedData_field(self.id, bytes(name, 'utf-8'), valuep)
            if flag < 0: raise Exception("could not get ceaSavedData field %s." % name)
            my_data[name] = valuep[0]
        buf = ffi.new("char[]", b'\000'*1024)
        flag = so.gas_state_get_ceaSavedData_species_names(self.id, buf, 1024)
        if flag < 0: raise Exception("could not get ceaSavedData species_names.")
        cea_species_names = ffi.string(buf).decode('utf-8').split("\t")
        massf_data = {}
        for name in cea_species_names:
            flag = so.gas_state_get_ceaSavedData_massf(self.id, bytes(name, 'utf-8'), valuep)
            if flag < 0: raise Exception("could not get ceaSavedData massf[%s]." % name)
            massf_data[name] = valuep[0]
        my_data["massf"] = massf_data
        return my_data

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

    def update_python_data(self):
        """
        There are places in the gas-dynamic functions where a new GasState
        has been filled in over in the DLang domain and
        all of the Python-domain data needs to be updated.
        """
        # In this wrapper class, all of the state data is in the Dlang domain.
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
    def gamma(self):
        return self.gmodel.gamma(self)
    @property
    def Prandtl(self):
        return self.gmodel.Prandtl(self)
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


# -----------------------------------------------------------------------------------
class PyGasState(object):
    """
    A GasState that has most of its data in the Python domain.

    This may be more efficient for Python programs
    that access the properties a lot in between updates.
    """
    __slots__ = ('gmodel', 'dgs', 'rho', 'p', 'T', 'u',
                 'n_modes', 'T_modes', 'u_modes', 'k_modes',
                 'a', 'n_species', '_massf', 'k', 'mu',
                 '_valuep', '_mf', '_modes', '_thermo_values')
    # Beyond the slots listed above, there are a number of properties defined below.
    # Together, these attributes allow the PyGasState object to look and behave like
    # a corresponding GasState object.

    def __init__(self, gmodel):
        self.gmodel = gmodel # Reference to the underlying Dlang GasModel
        self.dgs = GasState(gmodel) # Underlying Dlang GasState object
        # Python-domain attributes, for quick reference.
        self.rho = 0.0
        self.p = 0.0
        self.T = 0.0
        self.u = 0.0
        self.n_modes = gmodel.n_modes
        self.T_modes = [0.0,]*self.n_modes
        self.u_modes = [0.0,]*self.n_modes
        self.k_modes = [0.0,]*self.n_modes
        self.a = 0.0
        self.n_species = gmodel.n_species
        self.massf = [0.0,]*self.n_species if self.n_species > 1 else [1.0,]
        self.k = 0.0
        self.mu = 0.0
        # Allocate ffi buffers here, so that they can be reused.
        self._valuep = ffi.new("double *")
        self._mf = ffi.new("double[]", [0.0]*gmodel.n_species)
        self._modes = ffi.new("double[]", [0.0]*gmodel.n_modes)
        self._thermo_values = ffi.new("double[]", [0.0,]*5)
        return

    @property
    def id(self):
        return self.dgs.id

    def __copy__(self):
        gs_new = PyGasState(self.gmodel)
        gs_new.p = self.p
        gs_new.T = self.T
        if self.gmodel.n_modes > 0:
            gs_new.T_modes = self.T_modes
        gs_new.massf = self.massf
        if (self.p > 0.0) and (self.T > 0.0):
            gs_new.update_thermo_from_pT()
        return gs_new

    def __deepcopy__(self, memo):
        gs_new = PyGasState(self.gmodel)
        gs_new.p = self.p
        gs_new.T = self.T
        if self.gmodel.n_modes > 0:
            gs_new.T_modes = self.T_modes
        gs_new.massf = self.massf
        if (self.p > 0.0) and (self.T > 0.0):
            gs_new.update_thermo_from_pT()
        return gs_new

    def __str__(self):
        text = 'PyGasState(rho=%g' % self.rho
        text += ', p=%g' % self.p
        text += ', T=%g' % self.T
        text += ', u=%g' % self.u
        if len(self.T_modes) > 0:
            text += ', T_modes=%s' % self.T_modes
            text += ', u_modes=%s' % self.u_modes
            text += ', k_modes=%s' % self.k_modes
        text += ', a=%g' % self.a
        text += ', mu=%g' % self.mu
        text += ', k=%g' % self.k
        if self.n_species > 1:
            text += ', massf=%s' % str(self.massf)
        text += ', id=%d, gmodel.id=%d)' % (self.id, self.gmodel.id)
        return text

    def copy_thermo_properties_from_dgs(self):
        id = self.id
        flag = so.gas_state_get_thermo_scalars(id, self._thermo_values)
        if flag < 0: raise Exception("could not get thermo scalars from Dlang GasState")
        self.rho = self._thermo_values[0]
        self.p = self._thermo_values[1]
        self.T = self._thermo_values[2]
        self.u = self._thermo_values[3]
        self.a = self._thermo_values[4]
        # We assume that the mass fractions have not have changed.
        if self.n_modes > 0:
            flag = so.gas_state_get_array_field(self.id, b"T_modes", self._modes, self.n_modes)
            if flag < 0: raise Exception("could not get T_modes from Dlang GasState")
            self.T_modes = [self._modes[i] for i in range(self.n_modes)]
            flag = so.gas_state_get_array_field(self.id, b"u_modes", self._modes, self.n_modes)
            if flag < 0: raise Exception("could not get u_modes from Dlang GasState")
            self.u_modes = [self._modes[i] for i in range(self.n_modes)]
        return

    def copy_mass_fractions_into_dgs(self):
        # The mass fractions may have been specified as a dictionary
        # with species names as keys and Float values
        # or as a list of Float values.
        mf_list = self.massf.copy()
        if abs(sum(mf_list) - 1.0) > 1.0e-6:
            raise Exception("mass fractions do not sum to 1.")
        for i in range(self.n_species): self._mf[i] = mf_list[i]
        flag = so.gas_state_set_array_field(self.id, b"massf", self._mf, self.n_species)
        if flag < 0: raise Exception("could not set mass-fractions in Dlang GasState")
        return

    def update_thermo_from_pT(self):
        id = self.id
        flag = so.gas_state_set_scalar_field(id, b"p", self.p)
        if flag < 0: raise Exception("could not set pressure in Dlang GasState")
        flag = so.gas_state_set_scalar_field(id, b"T", self.T)
        if flag < 0: raise Exception("could not set temperature in Dlang GasState")
        if self.n_species > 1: self.copy_mass_fractions_into_dgs()
        if self.n_modes > 0:
            if type(self.T_modes) != type([]):
                raise Exception("T_modes needs to be supplied as a list.")
            if len(self.T_modes) != self.n_modes:
                raise Exception(f"T_modes list is not correct length. nmodes={n}; len(T_modes)={len(self.T_modes)}")
            for i in range(self.n_modes): self._modes[i] = self.T_modes[i]
            flag = so.gas_state_set_array_field(self.id, b"T_modes", self._modes, self.n_modes)
            if flag < 0: raise Exception("could not set T_modes in Dlang GasState")
        self.gmodel.update_thermo_from_pT(self.dgs)
        self.copy_thermo_properties_from_dgs()
        return

    def update_thermo_from_rhou(self):
        id = self.id
        flag = so.gas_state_set_scalar_field(id, b"rho", self.rho)
        if flag < 0: raise Exception("could not set density in Dlang GasState")
        flag = so.gas_state_set_scalar_field(id, b"u", self.u)
        if flag < 0: raise Exception("could not set internal-energy in Dlang GasState")
        if self.n_species > 1: self.copy_mass_fractions_into_dgs()
        if self.n_modes > 0:
            if type(self.u_modes) != type([]):
                raise Exception("u_modes needs to be supplied as a list.")
            if len(self.u_modes) != self.n_modes:
                raise Exception(f"u_modes list is not correct length. nmodes={n}; len(u_modes)={len(self.u_modes)}")
            for i in range(self.n_modes): self._modes[i] = self.u_modes[i]
            flag = so.gas_state_set_array_field(self.id, b"u_modes", self._modes, self.n_modes)
            if flag < 0: raise Exception("could not set u_modes in Dlang GasState")
        self.gmodel.update_thermo_from_rhou(self.dgs)
        self.copy_thermo_properties_from_dgs()
        return

    def update_thermo_from_rhoT(self):
        id = self.id
        flag = so.gas_state_set_scalar_field(id, b"rho", self.rho)
        if flag < 0: raise Exception("could not set density in Dlang GasState")
        flag = so.gas_state_set_scalar_field(id, b"T", self.T)
        if flag < 0: raise Exception("could not set temperature in Dlang GasState")
        if self.n_species > 1: self.copy_mass_fractions_into_dgs()
        if self.n_modes > 0:
            if type(self.T_modes) != type([]):
                raise Exception("T_modes needs to be supplied as a list.")
            if len(self.T_modes) != self.n_modes:
                raise Exception(f"T_modes list is not correct length. nmodes={n}; len(T_modes)={len(self.T_modes)}")
            for i in range(self.n_modes): self._modes[i] = self.T_modes[i]
            flag = so.gas_state_set_array_field(self.id, b"T_modes", self._modes, self.n_modes)
            if flag < 0: raise Exception("could not set T_modes in Dlang GasState")
        self.gmodel.update_thermo_from_rhoT(self.dgs)
        self.copy_thermo_properties_from_dgs()
        return

    def update_thermo_from_rhop(self):
        id = self.id
        flag = so.gas_state_set_scalar_field(id, b"rho", self.rho)
        if flag < 0: raise Exception("could not set density in Dlang GasState")
        flag = so.gas_state_set_scalar_field(id, b"p", self.p)
        if flag < 0: raise Exception("could not set pressure in Dlang GasState")
        if self.n_species > 1: self.copy_mass_fractions_into_dgs()
        if self.n_modes > 0: raise NotImplementedError('T_modes')
        self.gmodel.update_thermo_from_rhop(self.dgs)
        self.copy_thermo_properties_from_dgs()
        return

    def update_thermo_from_ps(self, s):
        id = self.id
        flag = so.gas_state_set_scalar_field(id, b"p", self.p)
        if flag < 0: raise Exception("could not set pressure in Dlang GasState")
        if self.n_species > 1: self.copy_mass_fractions_into_dgs()
        if self.n_modes > 0: raise NotImplementedError('T_modes')
        self.gmodel.update_thermo_from_ps(self.dgs, s)
        self.copy_thermo_properties_from_dgs()
        return

    def update_thermo_from_hs(self, h, s):
        if self.n_species > 1: self.copy_mass_fractions_into_dgs()
        if self.n_modes > 0: raise NotImplementedError('T_modes')
        self.gmodel.update_thermo_from_hs(self.dgs, h, s)
        self.copy_thermo_properties_from_dgs()
        return

    def update_sound_speed(self):
        self.dgs.gmodel.update_sound_speed(self.dgs)
        id = self.dgs.id
        flag = so.gas_state_get_scalar_field(id, b"a", self._valuep)
        if flag < 0: raise Exception("could not get sound-speed from Dlang GasState")
        self.a = self._valuep[0]
        return

    def update_trans_coeffs(self):
        self.dgs.gmodel.update_trans_coeffs(self.dgs)
        id = self.dgs.id
        flag = so.gas_state_get_scalar_field(id, b"mu", self._valuep)
        if flag < 0: raise Exception("could not get viscosity from Dlang GasState")
        self.mu = self._valuep[0]
        flag = so.gas_state_get_scalar_field(id, b"k", self._valuep)
        if flag < 0: raise Exception("could not get thermal-conductivity from Dlang GasState")
        self.k = self._valuep[0]
        if self.n_modes > 0:
            flag = so.gas_state_get_array_field(self.id, b"k_modes", self._modes, self.n_modes)
            if flag < 0: raise Exception("could not get k_modes from Dlang")
            self.k_modes = [self._modes[i] for i in range(self.n_modes)]
        return

    def update_python_data(self):
        """
        There are places in the gas-dynamic functions where a new GasState
        has been filled in over in the DLang domain and
        all of the Python-domain data needs to be updated.
        """
        id = self.id
        flag = so.gas_state_get_thermo_scalars(id, self._thermo_values)
        if flag < 0: raise Exception("could not get thermo scalars from Dlang GasState")
        self.rho = self._thermo_values[0]
        self.p = self._thermo_values[1]
        self.T = self._thermo_values[2]
        self.u = self._thermo_values[3]
        self.a = self._thermo_values[4]
        #
        flag = so.gas_state_get_array_field(id, b"massf", self._mf, self.n_species)
        if flag < 0: raise Exception("could not get mass-fractions from Dlang GasState")
        self.massf = [self._mf[i] for i in range(self.n_species)]
        #
        flag = so.gas_state_get_scalar_field(id, b"mu", self._valuep)
        if flag < 0: raise Exception("could not get viscosity from Dlang GasState")
        self.mu = self._valuep[0]
        flag = so.gas_state_get_scalar_field(id, b"k", self._valuep)
        if flag < 0: raise Exception("could not get thermal-conductivity from Dlang GasState")
        self.k = self._valuep[0]
        #
        if self.n_modes > 0:
            flag = so.gas_state_get_array_field(self.id, b"T_modes", self._modes, self.n_modes)
            if flag < 0: raise Exception("could not get T_modes from Dlang GasState")
            self.T_modes = [self._modes[i] for i in range(self.n_modes)]
            flag = so.gas_state_get_array_field(self.id, b"u_modes", self._modes, self.n_modes)
            if flag < 0: raise Exception("could not get u_modes from Dlang GasState")
            self.u_modes = [self._modes[i] for i in range(self.n_modes)]
            #
            flag = so.gas_state_get_array_field(self.id, b"k_modes", self._modes, self.n_modes)
            if flag < 0: raise Exception("could not get k_modes from Dlang GasState")
            self.k_modes = [self._modes[i] for i in range(self.n_modes)]
        return

    @property
    def massf(self):
        return self._massf
    @property
    def massf_as_dict(self):
        nsp = self.n_species
        names = self.gmodel.species_names
        mf = self.massf
        result = {}
        for i in range(nsp): result[names[i]] = mf[i]
        return result
    @massf.setter
    def massf(self, massf_given):
        if isinstance(massf_given, list):
            nsp = self.n_species
            if len(massf_given) != nsp:
                raise Exception(f"mass fraction list is not correct length. nsp={nsp}; len(massf)={len(mf_given)}")
            massf_list = massf_given.copy()
        elif isinstance(massf_given, dict):
            massf_list = [0.0] * self.n_species
            for ind, spcs in enumerate(self.gmodel.species_names):
                if spcs in massf_given.keys():
                    massf_list[ind] = massf_given[spcs]
        self._massf = massf_list

    @property
    def molef(self):
        """
        Mole fractions are always returned as a list.
        """
        nsp = self.n_species
        flag = so.gas_model_gas_state_get_molef(self.gmodel.id, self.id, self._mf)
        if flag < 0: raise Exception("could not get mole-fractions from Dlang GasState.")
        return [self._mf[i] for i in range(nsp)]
    @property
    def molef_as_dict(self):
        nsp = self.n_species
        names = self.gmodel.species_names
        mf = self.molef
        result = {}
        for i in range(nsp): result[names[i]] = mf[i]
        return result
    @molef.setter
    def molef(self, molef_given):
        """
        Mole fractions may be provided as a list or dictionary.

        The object does not actually store the mole fractions but,
        instead, the corresponding mass-fractions
        """
        self.massf = self.gmodel.molef2massf(molef_given)
        return None

    @property
    def conc(self):
        nsp = self.n_species
        flag = so.gas_model_gas_state_get_conc(self.gmodel.id, self.id, self._mf)
        if flag < 0: raise Exception("could not get concentrations from Dlang GasState.")
        return [self._mf[i] for i in range(nsp)]
    @property
    def conc_as_dict(self):
        nsp = self.n_species
        names = self.gmodel.species_names
        conc_list = self.conc
        result = {}
        for i in range(nsp): result[names[i]] = conc_list[i]
        return result

    def copy_values(self, gstate):
        self.gmodel = gstate.gmodel
        flag = so.gas_state_copy_values(self.id, gstate.id)
        if flag < 0: raise Exception("could not copy values in Dlang GasStates")
        self.update_python_data()
        return

    # For the following properties, we just delegate the work
    # to the underlying Dlang GasState object, presuming that
    # it has consistent data.
    @property
    def Cv(self):
        return self.gmodel.Cv(self.dgs)
    @property
    def Cp(self):
        return self.gmodel.Cp(self.dgs)
    @property
    def dpdrho_const_T(self):
        return self.gmodel.dpdrho_const_T(self.dgs)
    @property
    def R(self):
        return self.gmodel.R(self.dgs)
    @property
    def gamma(self):
        return self.gmodel.gamma(self.dgs)
    @property
    def Prandtl(self):
        return self.gmodel.Prandtl(self.dgs)
    @property
    def internal_energy(self):
        return self.gmodel.internal_energy(self.dgs)
    @property
    def enthalpy(self):
        return self.gmodel.enthalpy(self.dgs)
    @property
    def entropy(self):
        return self.gmodel.entropy(self.dgs)
    @property
    def molecular_mass(self):
        return self.gmodel.molecular_mass(self.dgs)

    def enthalpy_isp(self, isp):
        return self.gmodel.enthalpy_isp(self.dgs, isp)
    def entropy_isp(self, isp):
        return self.gmodel.entropy_isp(self.dgs, isp)
    def gibbs_free_energy_isp(self, isp):
        return self.gmodel.gibbs_free_energy_isp(self.dgs, isp)


# -----------------------------------------------------------------------------------
class ThermochemicalReactor(object):
    def __init__(self, gmodel, filename1, filename2=""):
        self.filename1 = filename1
        self.filename2 = filename2
        self.gmodel = gmodel
        self.id = so.thermochemical_reactor_new(self.gmodel.id,
                                                bytes(self.filename1, 'utf-8'),
                                                bytes(self.filename2, 'utf-8'))
        self.dt_suggestp = ffi.new("double *")

    def __copy__(self):
        reactor_new = ThermochemicalReactor(self.gmodel, self.filename1, self.filename2)
        return reactor_new

    def __deepcopy__(self, memo):
        reactor_new = ThermochemicalReactor(self.gmodel, self.filename1, self.filename2)
        return reactor_new

    def __str__(self):
        text = 'ThermochemicalReactor(id=%d, gmodel.id=%d, file1="%s", file2="%s")' % \
            (self.id, self.gmodel.id, self.filename1, self.filename2)
        return text

    def update_state(self, gstate, t_interval, dt_suggest=-1.0):
        self.dt_suggestp[0] = dt_suggest
        flag = so.thermochemical_reactor_gas_state_update(self.id, gstate.id,
                                                          t_interval, self.dt_suggestp)
        if flag < 0: raise Exception("could not update state.")
        gstate.update_python_data()
        return self.dt_suggestp[0]

    def source_terms(self, gstate):
        source = ffi.new("double[]", [0.0]*(self.gmodel.n_species + self.gmodel.n_modes))
        flag = so.thermochemical_reactor_eval_source_terms(self.id, self.gmodel.id, gstate.id,
                                                           self.gmodel.n_species, self.gmodel.n_modes, source)
        if flag < 0: raise Exception("Could not compute source terms.")
        return [source[i] for i in range(self.gmodel.n_species+self.gmodel.n_modes)]


# -----------------------------------------------------------------------------------
class ReactionMechanism(object):
    def __init__(self, gmodel, filename):
        self.filename = filename
        self.gmodel = gmodel
        self.id = so.reaction_mechanism_new(self.gmodel.id,
                                                bytes(self.filename, 'utf-8'))
    def __str__(self):
        text = 'ReactionMechanism(id=%d, gmodel.id=%d, file="%s")' % \
            (self.id, self.gmodel.id, self.filename)
        return text

    @property
    def n_reactions(self):
        return so.reaction_mechanism_n_reactions(self.id)

    def reaction_tickrates(self, gstate):
        forwardrates = ffi.new("double[]", [0.0]*self.n_reactions)
        backwardrates = ffi.new("double[]", [0.0]*self.n_reactions)
        flag = so.reaction_mechanism_tickrates(self.id, self.gmodel.id, gstate.id, forwardrates, backwardrates);

        if flag < 0: raise Exception("Could not compute reaction tickrates.")
        return [forwardrates[i] - backwardrates[i] for i in range(self.n_reactions)]

    def forward_tickrates(self, gstate):
        forwardrates = ffi.new("double[]", [0.0]*self.n_reactions)
        backwardrates = ffi.new("double[]", [0.0]*self.n_reactions)
        flag = so.reaction_mechanism_tickrates(self.id, self.gmodel.id, gstate.id, forwardrates, backwardrates);

        if flag < 0: raise Exception("Could not compute reaction tickrates.")
        return [forwardrates[i] for i in range(self.n_reactions)]

    def backward_tickrates(self, gstate):
        forwardrates = ffi.new("double[]", [0.0]*self.n_reactions)
        backwardrates = ffi.new("double[]", [0.0]*self.n_reactions)
        flag = so.reaction_mechanism_tickrates(self.id, self.gmodel.id, gstate.id, forwardrates, backwardrates);

        if flag < 0: raise Exception("Could not compute reaction tickrates.")
        return [backwardrates[i] for i in range(self.n_reactions)]


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
        state2.update_python_data()
        v2 = my_results[0]
        vg = my_results[1]
        return [v2, vg]

    def normal_shock(self, state1, vs, state2, rho_tol=1.0e-6, T_tol=0.1):
        my_results = ffi.new("double[]", [0.0]*2)
        flag = so.gasflow_normal_shock(state1.id, vs, state2.id, self.gmodel.id, my_results,
                                       rho_tol, T_tol)
        if flag < 0: raise Exception("failed to compute normal shock jump.")
        state2.update_python_data()
        v2 = my_results[0]
        vg = my_results[1]
        return [v2, vg]

    def normal_shock_1(self, state1, vs, state2, p_tol=0.5, T_tol=0.1):
        my_results = ffi.new("double[]", [0.0]*2)
        flag = so.gasflow_normal_shock_1(state1.id, vs, state2.id, self.gmodel.id, my_results,
                                         p_tol, T_tol)
        if flag < 0: raise Exception("failed to compute normal shock jump.")
        state2.update_python_data()
        v2 = my_results[0]
        vg = my_results[1]
        return [v2, vg]

    def normal_shock_p2p1(self, state1, p2p1, state2):
        my_results = ffi.new("double[]", [0.0]*3)
        flag = so.gasflow_normal_shock_p2p1(state1.id, p2p1, state2.id, self.gmodel.id,
                                            my_results)
        if flag < 0: raise Exception("failed to compute normal shock jump from p2p1.")
        state2.update_python_data()
        vs = my_results[0]
        v2 = my_results[1]
        vg = my_results[2]
        return [vs, v2, vg]

    def reflected_shock(self, state2, vg, state5):
        my_results = ffi.new("double[]", [0.0])
        flag = so.gasflow_reflected_shock(state2.id, vg, state5.id, self.gmodel.id,
                                          my_results)
        if flag < 0: raise Exception("failed to compute reflected shock.")
        state5.update_python_data()
        vr = my_results[0]
        return vr

    def expand_from_stagnation(self, state0, p_over_p0, state1):
        my_results = ffi.new("double[]", [0.0])
        flag = so.gasflow_expand_from_stagnation(state0.id, p_over_p0, state1.id,
                                                 self.gmodel.id, my_results)
        if flag < 0: raise Exception("failed to compute expansion from stagnation.")
        state1.update_python_data()
        v = my_results[0]
        return v

    def expand_to_mach(self, state0, mach, state1):
        my_results = ffi.new("double[]", [0.0])
        flag = so.gasflow_expand_to_mach(state0.id, mach, state1.id,
                                         self.gmodel.id, my_results)
        if flag < 0: raise Exception("failed to compute expansion to mach number.")
        state1.update_python_data()
        v = my_results[0]
        return v

    def total_condition(self, state1, v1, state0):
        flag = so.gasflow_total_condition(state1.id, v1, state0.id, self.gmodel.id)
        if flag < 0: raise Exception("failed to compute total condition.")
        state0.update_python_data()
        return

    def pitot_condition(self, state1, v1, state2pitot):
        flag = so.gasflow_pitot_condition(state1.id, v1, state2pitot.id, self.gmodel.id)
        if flag < 0: raise Exception("failed to compute pitot condition.")
        state2pitot.update_python_data()
        return

    def steady_flow_with_area_change(self, state1, v1, area2_over_area1, state2,
                                     tol=1.0e-4, p2p1_min=0.0001):
        my_results = ffi.new("double[]", [0.0])
        flag = so.gasflow_steady_flow_with_area_change(state1.id, v1, area2_over_area1,
                                                       state2.id, self.gmodel.id, tol,
                                                       p2p1_min, my_results)
        if flag < 0: raise Exception("failed to compute steady flow with area change.")
        state2.update_python_data()
        v2 = my_results[0]
        return v2

    def finite_wave_dp(self, state1, v1, characteristic, p2, state2, steps=100):
        char_name = bytes(characteristic, 'utf-8')
        my_results = ffi.new("double[]", [0.0])
        flag = so.gasflow_finite_wave_dp(state1.id, v1, char_name, p2,
                                         state2.id, self.gmodel.id, steps,
                                         my_results)
        if flag < 0: raise Exception("failed to compute (unsteady) finite wave dp.")
        state2.update_python_data()
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
        state2.update_python_data()
        v2 = my_results[0]
        return v2

    def osher_riemann(self, stateL, stateR, velL, velR, stateLstar, stateRstar, stateX0):
        my_results = ffi.new("double[]", [0.0]*5)
        flag = so.gasflow_osher_riemann(stateL.id, stateR.id, velL, velR,
                                        stateLstar.id, stateRstar.id,
                                        stateX0.id, self.gmodel.id, my_results)
        if flag < 0: raise Exception("failed to compute solution to Riemann problem.")
        stateLstar.update_python_data()
        stateRstar.update_python_data()
        stateX0.update_python_data()
        pstar = my_results[0]
        wstar = my_results[1]
        wL = my_results[2]
        wR = my_results[3]
        velX0 = my_results[4]
        return [pstar, wstar, wL, wR, velX0]

    def osher_flux(self, stateL, stateR, velL, velR):
        my_results = ffi.new("double[]", [0.0]*3)
        flag = so.gasflow_osher_flux(stateL.id, stateR.id, velL, velR, self.gmodel.id, my_results)
        if flag < 0: raise Exception("failed to compute Osher flux.")
        F_mass = my_results[0]
        F_x_momentum = my_results[1]
        F_energy = my_results[2]
        return [F_mass, F_x_momentum, F_energy]

    def roe_flux(self, stateL, stateR, velL, velR):
        my_results = ffi.new("double[]", [0.0]*3)
        flag = so.gasflow_roe_flux(stateL.id, stateR.id, velL, velR, self.gmodel.id, my_results)
        if flag < 0: raise Exception("failed to compute Roe flux.")
        F_mass = my_results[0]
        F_x_momentum = my_results[1]
        F_energy = my_results[2]
        return [F_mass, F_x_momentum, F_energy]

    def lrivp(self, stateL, stateR, velL, velR):
        my_pstar = ffi.new("double[]", [0.0])
        my_wstar = ffi.new("double[]", [0.0])
        flag = so.gasflow_lrivp(stateL.id, stateR.id, velL, velR,
                                stateL.gmodel.id, stateR.gmodel.id,
                                my_wstar, my_pstar)
        if flag < 0: raise Exception("failed to compute solution for lrivp.")
        pstar = my_pstar[0]
        wstar = my_wstar[0]
        return [pstar, wstar]

    def piston_at_left(self, stateR, velR, wstar):
        my_pstar = ffi.new("double[]", [0.0])
        flag = so.gasflow_piston_at_left(stateR.id, velR, stateR.gmodel.id,
                                         wstar, my_pstar)
        if flag < 0: raise Exception("failed to compute solution for piston_at_left.")
        return my_pstar[0]

    def piston_at_right(self, stateL, velL, wstar):
        my_pstar = ffi.new("double[]", [0.0])
        flag = so.gasflow_piston_at_right(stateL.id, velL, stateL.gmodel.id,
                                          wstar, my_pstar)
        if flag < 0: raise Exception("failed to compute solution for piston_at_right.")
        return my_pstar[0]

    def theta_oblique(self, state1, v1, beta, state2):
        my_results = ffi.new("double[]", [0.0, 0.0])
        flag = so.gasflow_theta_oblique(state1.id, v1, beta,
                                        state2.id, self.gmodel.id, my_results)
        if flag < 0: raise Exception("failed to compute theta oblique.")
        state2.update_python_data()
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

    def theta_cone(self, state1, v1, beta, state_c, dtheta=-0.01*math.pi/180.0):
        my_results = ffi.new("double[]", [0.0, 0.0])
        flag = so.gasflow_theta_cone(state1.id, v1, beta,
                                     state_c.id, self.gmodel.id,
                                     dtheta, my_results)
        if flag < 0: raise Exception("failed to compute theta cone.")
        state_c.update_python_data()
        theta_c = my_results[0]
        v2_c = my_results[1]
        return theta_c, v2_c

    def beta_cone(self, state1, v1, theta, dtheta=-0.01*math.pi/180.0):
        my_results = ffi.new("double[]", [0.0])
        flag = so.gasflow_beta_cone(state1.id, v1, theta,
                                    self.gmodel.id, dtheta, my_results)
        if flag < 0: raise Exception("failed to compute beta cone.")
        beta = my_results[0]
        return beta
