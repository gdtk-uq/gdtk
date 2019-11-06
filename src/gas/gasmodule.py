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
    char* gas_model_species_name(int gm_i, int isp, char* name, int n);
    int gas_model_mol_masses(int gm_i, double* mm);

    int gas_state_new(int gm_i);
    int gas_state_set_scalar_field(int gs_i, char* field_name, double value);
    double gas_state_get_scalar_field(int gs_i, char* field_name);
    int gas_state_set_array_field(int gs_i, char* field_name, double* values, int n);
    int gas_state_get_array_field(int gs_i, char* field_name, double* values, int n);

    int gas_model_gas_state_update_thermo_from_pT(int gm_i, int gs_i);
    int gas_model_gas_state_update_thermo_from_rhou(int gm_i, int gs_i);
    int gas_model_gas_state_update_thermo_from_rhoT(int gm_i, int gs_i);
    int gas_model_gas_state_update_thermo_from_rhop(int gm_i, int gs_i);
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
        return so.gas_state_get_scalar_field(self.id, b"rho")
    @rho.setter
    def rho(self, value):
        flag = so.gas_state_set_scalar_field(self.id, b"rho", value)
        if flag < 0: raise Exception("could not set density.")
        return
    
    @property
    def p(self):
        return so.gas_state_get_scalar_field(self.id, b"p")
    @p.setter
    def p(self, value):
        flag = so.gas_state_set_scalar_field(self.id, b"p", value)
        if flag < 0: raise Exception("could not set pressure.")
        return
    
    @property
    def T(self):
        return so.gas_state_get_scalar_field(self.id, b"T")
    @T.setter
    def T(self, value):
        flag = so.gas_state_set_scalar_field(self.id, b"T", value)
        if flag < 0: raise Exception("could not set temperature.")
        return
    
    @property
    def u(self):
        return so.gas_state_get_scalar_field(self.id, b"u")
    @u.setter
    def u(self, value):
        flag = so.gas_state_set_scalar_field(self.id, b"u", value)
        if flag < 0: raise Exception("could not set internal energy.")
        return

    @property
    def massf(self):
        nsp = self.gmodel.n_species
        mf = ffi.new("double[]", [0.0]*nsp)
        so.gas_state_get_array_field(self.id, b"massf", mf, nsp)
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
        so.gas_state_set_array_field(self.id, b"massf", mf, nsp)
        return mf_list
            
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
