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

    int gas_state_new(int gm_i);
    int gas_state_set_scalar_field(int gs_i, char* field_name, double value);
    double gas_state_get_scalar_field(int gs_i, char* field_name);

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

    def update_thermo_from_pT(self, gstate):
        flag = so.gas_model_gas_state_update_thermo_from_pT(self.id, gstate.id)
        if flag < 0: raise Exception("Oops, could not update thermo from p,T.")
        return

    def update_thermo_from_rhou(self, gstate):
        flag = so.gas_model_gas_state_update_thermo_from_rhou(self.id, gstate.id)
        if flag < 0: raise Exception("Oops, could not update thermo from rho,u.")
        return

    def update_thermo_from_rhoT(self, gstate):
        flag = so.gas_model_gas_state_update_thermo_from_rhoT(self.id, gstate.id)
        if flag < 0: raise Exception("Oops, could not update thermo from rho,T.")
        return

    def update_thermo_from_rhop(self, gstate):
        flag = so.gas_model_gas_state_update_thermo_from_rhop(self.id, gstate.id)
        if flag < 0: raise Exception("Oops, could not update thermo from rho,p.")
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
        text += ', id=%d, gmodel.id=%d)' % (self.id, self.gmodel.id)
        return text
    
    @property
    def rho(self):
        return so.gas_state_get_scalar_field(self.id, b"rho")
    @rho.getter
    def rho(self):
        return so.gas_state_get_scalar_field(self.id, b"rho")
    @rho.setter
    def rho(self, value):
        flag = so.gas_state_set_scalar_field(self.id, b"rho", value)
        if flag < 0: raise Exception("Oops, could not set density.")
        return
    
    @property
    def p(self):
        return so.gas_state_get_scalar_field(self.id, b"p")
    @p.getter
    def p(self):
        return so.gas_state_get_scalar_field(self.id, b"p")
    @p.setter
    def p(self, value):
        flag = so.gas_state_set_scalar_field(self.id, b"p", value)
        if flag < 0: raise Exception("Oops, could not set pressure.")
        return
    
    @property
    def T(self):
        return so.gas_state_get_scalar_field(self.id, b"T")
    @T.getter
    def T(self):
        return so.gas_state_get_scalar_field(self.id, b"T")
    @T.setter
    def T(self, value):
        flag = so.gas_state_set_scalar_field(self.id, b"T", value)
        if flag < 0: raise Exception("Oops, could not set temperature.")
        return
    
    @property
    def u(self):
        return so.gas_state_get_scalar_field(self.id, b"u")
    @u.getter
    def u(self):
        return so.gas_state_get_scalar_field(self.id, b"u")
    @u.setter
    def u(self, value):
        flag = so.gas_state_set_scalar_field(self.id, b"u", value)
        if flag < 0: raise Exception("Oops, could not set internal energy.")
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
