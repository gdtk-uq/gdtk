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

    int gas_state_new(int gm_i);
    int gas_state_set_scalar_field(int gs_i, char* field_name, double value);
    double gas_state_get_scalar_field(int gs_i, char* field_name);
""")
so = ffi.dlopen("libgasmodule.so")
so.cwrap_gas_module_init()

# Service classes that wrap the C-API in a nice Pythonic API...

class GasModel(object):
    def __init__(self, file_name):
        self.file_name = file_name
        self.id = so.gas_model_new(bytes(self.file_name, 'utf-8'))

    def __str__(self):
        text = 'GasModel(file="%s", id=%d)' % (self.file_name, self.id)
        return text
    
    @property
    def n_species(self):
        return so.gas_model_n_species(self.id)


class GasState(object):
    def __init__(self, gmodel):
        self.gmodel = gmodel
        self.id = so.gas_state_new(self.gmodel.id)

    def __str__(self):
        text = 'GasState(rho=%g' % self.rho
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
