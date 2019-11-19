# gasmodule.py
# A Python wrapper for GasModel and GasState classes.
# PJ 2019-07-24: start of experiment with FFI.
#    2019-07-25: added Python wrapper
#
from cffi import FFI

ffi = FFI()
ffi.cdef("""
    int cwrap_gas_module_init();

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

    int gas_model_gas_state_update_thermo_from_pT(int gm_i, int gs_i);
    int gas_model_gas_state_update_thermo_from_rhou(int gm_i, int gs_i);
    int gas_model_gas_state_update_thermo_from_rhoT(int gm_i, int gs_i);
    int gas_model_gas_state_update_thermo_from_rhop(int gm_i, int gs_i);
    int gas_model_gas_state_update_thermo_from_ps(int gm_i, int gs_i, double s);
    int gas_model_gas_state_update_thermo_from_hs(int gm_i, int gs_i, double h, double s);
    int gas_model_gas_state_update_sound_speed(int gm_i, int gs_i);
    int gas_model_gas_state_update_trans_coeffs(int gm_i, int gs_i);
""")
so = ffi.dlopen("libgasmodule.so")
so.cwrap_gas_module_init()

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
        text = 'GasModel(file="%s", id=%d)' % (self.file_name, self.id)
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
        flag = so.gas_model_gas_state_update_thermo_from_rhop(self.id, gstate.id, s)
        if flag < 0: raise Exception("could not update thermo from p,s.")
        return

    def update_thermo_from_hs(self, gstate, h, s):
        flag = so.gas_model_gas_state_update_thermo_from_rhop(self.id, gstate.id, h, s)
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


class GasState(object):
    def __init__(self, gmodel):
        self.gmodel = gmodel
        self.id = so.gas_state_new(self.gmodel.id)

    def __str__(self):
        text = 'GasState(rho=%g' % self.rho
        text += ', p=%g' % self.p
        text += ', T=%g' % self.T
        text += ', u=%g' % self.u
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
