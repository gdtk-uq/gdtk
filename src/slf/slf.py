"""
Python code for interfacing with slf.

@author: Nick Gibbons
"""

from cffi import FFI
from numpy import array, zeros
from os import environ, path

ffi = FFI()
ffi.cdef("""
    int cwrap_init();
    int init_slf(const char* file_name);
    int set_initial_condition();
    int run();
    int get_size_of_solution();
    int get_solution(double* U);
""")
#so = ffi.dlopen("libslf.so")
#so.cwrap_init()

class Flame(object):
    def __init__(self, file):
        self.file = file
        #check = so.init_slf(bytes(file, 'utf-8'))

        dgdpath = environ.get('DGD')
        if dgdpath==None: raise Exception("DGD install not found!")
        libpath = path.join(dgdpath, 'lib', 'libslf.so')
        print("got libpath: ", libpath)
        self.lib = ffi.dlopen(libpath)
        self.lib.cwrap_init()
        check = self.lib.init_slf(bytes(file, 'utf-8'))
        return

    def set_initial_condition(self):
        result = self.lib.set_initial_condition()
        return result

    def run(self):
        result = self.lib.run()
        return result

    def get_solution(self):
        n = self.lib.get_size_of_solution()
        U = zeros(n)
        Up = ffi.from_buffer("double[]", U)
        self.lib.get_solution(Up)
        return U

if __name__=='__main__':
    print("Put a test here, I guess...")
