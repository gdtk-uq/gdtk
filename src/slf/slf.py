"""
Python code for interfacing with slf.

@author: Nick Gibbons
"""

from cffi import FFI
from numpy import array, zeros, interp, frombuffer
from os import environ, path
from numpy.polynomial import Polynomial
import struct
from io import BytesIO, BufferedWriter

header = """
    int cwrap_init();
    int init_slf(const char* file_name);
    int set_initial_condition();
    int run();
    int save_solution();
    int save_log();
    int load_solution();
    int back_out_scalar_dissipation();
    int extract_residual(double*);
    int get_species_name(int, char*, int);

    int get_nsp();
    int get_neq();
    int get_N();
    int get_n();

    double get_D();
    double get_p();
    double get_dZ();
    double get_T0();
    double get_T1();
    double get_targetGRR();
    double get_lewis_number();

    void set_D(double D);
    void set_p(double p);
    void set_targetGRR(double targetGRR);
    void set_lewis_number(double lewis_number);

    int get_chi(double*);
    int get_Z(double*);
    int get_Y0(double*);
    int get_Y1(double*);
    int get_U0(double*);
    int get_U1(double*);
    int get_U(double*);
    int get_R(double*);
    int get_T(double*);
    int get_Y(double*);

    int set_U(double*);
    int set_chi(double*);
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
        if check!=0: raise Exception("Failed to initialize slf flame object")

        self.species_names = []
        buf = ffi.new("char[]", b'\000'*32)
        for isp in range(self.nsp):
            self.lib.get_species_name(isp, buf, 32)
            self.species_names.append(ffi.string(buf).decode('utf-8'))
        return

    def set_initial_condition(self):
        result = self.lib.set_initial_condition()
        return result

    def run(self):
        result = self.lib.run()
        if result!=0: raise Exception("Run failed for some reason.")
        return result

    def save_solution(self):
        result = self.lib.save_solution()
        if result!=0: raise Exception("Failed to save slf solution to disk")
        return

    def save_log(self):
        result = self.lib.save_log()
        if result!=0: raise Exception("Failed to save slf log file to disk")
        return

    def load_solution(self):
        result = self.lib.load_solution()
        if result!=0: raise Exception("Failed to load slf solution from disk")

    def back_out_scalar_dissipation(self):
        result = self.lib.back_out_scalar_dissipation() 
        if result!=0: raise Exception("Failed to reverse engineer a scalar dissipation")
        return

    def extract_residual(self):
        R = zeros(self.n)
        Rp = ffi.cast("double *", ffi.from_buffer(R))
        result = self.lib.extract_residual(Rp)
        if result!=0: raise Exception("Failed to extract a residual vector")
        return R.reshape((self.N, self.neq))

    def create_polynomial_fit(self):
        """
        To get smooth gradients and double gradients, we need to fit some curves
        to the flame solution.
        """
        self.polyfits = {}

        Y = self.Y
        for isp in range(self.nsp):
            sp = self.species_names[isp]
            self.polyfits['massf-{}'.format(sp)] = Polynomial.fit(self.Z, Y[:,isp], deg=4)

        self.polyfits['T'] = Polynomial.fit(self.Z, self.T, deg=4)
        return

    def return_linear_interps(self, Zs):
        interps = {}
        Y = self.Y
        Z = self.Z
        Y0 = self.Y0
        Y1 = self.Y1

        for isp in range(self.nsp):
            sp = self.species_names[isp]
            key = 'massf-{}'.format(sp)
            interps[key] = interp(Zs, Z, Y[:,isp], left=Y0[isp], right=Y1[isp])
        interps['T'] = interp(Zs, Z, self.T, left=self.T0, right=self.T1)

        return interps


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

    @property
    def D(self):
        return self.lib.get_D()

    @D.setter
    def D(self, value):
        self.lib.set_D(value)

    @property
    def p(self):
        return self.lib.get_p()

    @p.setter
    def p(self, value):
        self.lib.set_p(value)

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
    def targetGRR(self):
        return self.lib.get_targetGRR()

    @targetGRR.setter
    def targetGRR(self, value):
        self.lib.set_targetGRR(value)

    @property
    def lewis_number(self):
        return self.lib.get_lewis_number()

    @lewis_number.setter
    def lewis_number(self, value):
        self.lib.set_lewis_number(value)

    @property
    def Z(self):
        Z = zeros(self.N)
        Zp = ffi.cast("double *", ffi.from_buffer(Z))
        code = self.lib.get_Z(Zp)
        if code!=0: raise Exception("Failed to get array Z from slf")
        return Z

    @property
    def chi(self):
        chi = zeros(self.N)
        chip = ffi.cast("double *", ffi.from_buffer(chi))
        code = self.lib.get_chi(chip)
        if code!=0: raise Exception("Failed to get array chi from slf")
        return chi

    @chi.setter
    def chi(self, arr):
        chip = ffi.cast("double *", ffi.from_buffer(arr))
        code = self.lib.set_chi(chip)
        if code!=0: raise Exception("Failed to set array chi from slf")
        return

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
        return U.reshape((self.N, self.neq))

    @U.setter
    def U(self, arr):
        Up = ffi.cast("double *", ffi.from_buffer(arr))
        code = self.lib.set_U(Up)
        if code!=0: raise Exception("Failed to set array U in slf")
        return

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
        return Y.reshape((self.N, self.nsp))

def number(thing):
    try:
        outnumber = int(thing)
    except(ValueError):
        try:
            outnumber = float(thing)
        except(ValueError):
            outnumber = complex(thing.replace('i','j')).real
    return outnumber

def read_log_file(filename):
    with open(filename) as fp:
       lines = [line.strip().split() for line in fp if line.startswith('iter')]
    
    tokens = lines[0][::2]
    data = {token:[] for token in tokens}

    for line in lines:
        tokens = line[::2]
        values  = list(map(number, line[1::2]))
        for t,d in zip(tokens, values):
            data[t].append(d)
    return data

def read_solution_file(filename):
    with open(filename, 'rb') as fp:
        bytes = fp.read()

    stream = BytesIO(bytes)

    buff = stream.read(8*4)
    nsp, neq, N, n = struct.unpack("Q"*4, buff)

    buff = stream.read(8*5)
    D, p, dZ, T0, T1 = struct.unpack("d"*5, buff)

    buff = stream.read(8*N)
    Z = frombuffer(buff)

    buff = stream.read(8*nsp)
    Y0 = frombuffer(buff)

    buff = stream.read(8*nsp)
    Y1 = frombuffer(buff)

    buff = stream.read(8*n)
    U = frombuffer(buff)
    U = U.reshape((N, neq))
    Y = U[:,:-1]
    T = U[:,-1]

    data = {}
    data['nsp'] = nsp; data['neq'] = neq; data['N'] = N; data['n'] = n
    data['D'] = p; data['p'] = p; data['dZ'] = dZ; data['T0'] = T0; data['T1'] = T1
    data['Z'] = Z; data['Y0'] = Y0; data['Y1'] = Y1
    data['Y'] = Y; data['T'] = T
    return data

if __name__=='__main__':
    print("Put a test here, I guess...")
