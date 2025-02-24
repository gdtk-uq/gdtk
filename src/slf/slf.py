"""
Python code for interfacing with slf.

@author: Nick Gibbons
"""

from cffi import FFI
from numpy import array, zeros
from os import environ, path

header = """
    int cwrap_init();
    int init_slf(const char* file_name);
    int set_initial_condition();
    int run();

    int get_nsp();
    int get_neq();
    int get_N();
    int get_n();

    double get_D();
    double get_p();
    double get_dZ();
    double get_T0();
    double get_T1();

    int get_Z(double*);
    int get_Y0(double*);
    int get_Y1(double*);
    int get_U0(double*);
    int get_U1(double*);
    int get_U(double*);
    int get_R(double*);
    int get_T(double*);
    int get_Y(double*);
"""

ffi = FFI()
ffi.cdef(header)

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

    @property
    def nsp(self):
        return self.lib.get_nsp()

    @property
    def neq(self):
        return self.lib.get_neq()

    @property
    def N(self):
        return self.lib.get_N()

    @property
    def n(self):
        return self.lib.get_n()

    # Put some setters for these, I think
    @property
    def D(self):
        return self.lib.get_D()

    @property
    def p(self):
        return self.lib.get_p()

    @property
    def dZ(self):
        return self.lib.get_dZ()

    @property
    def T0(self):
        return self.lib.get_T0()

    @property
    def T1(self):
        return self.lib.get_T1()

    @property
    def Z(self):
        Z = zeros(self.N)
        code = self.lib.get_Z(Z)
        if code!=0: raise Exception("Failed to get array Z from slf")
        return Z

    @property
    def Y0(self):
        Y0 = zeros(self.nsp)
        Y0p = ffi.cast("double *", ffi.from_buffer(Y0))
        code = self.lib.get_Y0(Y0p)
        if code!=0: raise Exception("Failed to get array Y0 from slf")
        return Y0

    @property
    def Y1(self):
        Y1 = zeros(self.nsp)
        Y1p = ffi.cast("double *", ffi.from_buffer(Y1))
        code = self.lib.get_Y1(Y1p)
        if code!=0: raise Exception("Failed to get array Y1 from slf")
        return Y1

    @property
    def U0(self):
        U0 = zeros(self.nsp+1)
        U0p = ffi.cast("double *", ffi.from_buffer(U0))
        code = self.lib.get_U0(U0p)
        if code!=0: raise Exception("Failed to get array U0 from slf")
        return U0

    @property
    def U1(self):
        U1 = zeros(self.nsp+1)
        U1p = ffi.cast("double *", ffi.from_buffer(U1))
        code = self.lib.get_U1(U1p)
        if code!=0: raise Exception("Failed to get array U1 from slf")
        return U1

    @property
    def U(self):
        U = zeros(self.n)
        Up = ffi.cast("double *", ffi.from_buffer(U))
        code = self.lib.get_U(Up)
        if code!=0: raise Exception("Failed to get array U from slf")
        return U

    @property
    def R(self):
        R = zeros(self.n)
        Rp = ffi.cast("double *", ffi.from_buffer(R))
        code = self.lib.get_R(Rp)
        if code!=0: raise Exception("Failed to get array R from slf")
        return R

    @property
    def T(self):
        T = zeros(self.N)
        Tp = ffi.cast("double *", ffi.from_buffer(T))
        code = self.lib.get_T(Tp)
        if code!=0: raise Exception("Failed to get array T from slf")
        return T

    @property
    def Y(self):
        Y = zeros(self.N*self.nsp)
        Yp = ffi.cast("double *", ffi.from_buffer(Y))
        code = self.lib.get_Y(Yp)
        if code!=0: raise Exception("Failed to get array Y from slf")
        return Y

if __name__=='__main__':
    print("Put a test here, I guess...")
