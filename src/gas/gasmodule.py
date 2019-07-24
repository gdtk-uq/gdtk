# gasmodule.py
# A Python wrapper for GasModel and GasState classes.
# PJ 2019-07-24: start of experiment.
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

# Should now introduce a service class that wraps the C-API
# in a nice Pythonic API...

