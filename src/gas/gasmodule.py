# gasmodule.py
# A Python wrapper for GasModel and GasState classes.
# PJ 2019-07-24: start of experiment.
#
from cffi import FFI
ffi = FFI()
ffi.cdef("""
    int cwrap_gas_module_init();
""")
so = ffi.dlopen("libgasmodule.so")
so.cwrap_gas_module_init()

# Should now introduce a service class that wraps the C-API
# in a nice Pythonic API...

