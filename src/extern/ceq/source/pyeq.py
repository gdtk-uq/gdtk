"""
Python interface for ceq chemical equilibrium calculator

References:
    "Computer Program for Calculation of Complex Equilibrium Compositions and Applications"
    Nasa Reference Publication 1311, October 1995
    Sanford Gordon and Bonnie J. McBride

    "NASA Glenn Coefficients for Calculating Thermodynamic Properties of Individual Species"
    NASA/TP - 2002-211556, September 2002
    Bonnie J. McBride, Michael J. Zehe, and Sanford Gordon


@author: Nick Gibbons (n.gibbons(at)uq.edu.au)
"""

from string import ascii_letters
from numpy import array, zeros, log
from ctypes import cdll,c_double,POINTER,c_int,byref

letters = set(ascii_letters)
DBPATH='../thermo.inp'
LIBPATH='./libceq.so'

class EqCalculator(object):
    """ Python interface to low level ceq routines """
    def __init__(self, spnames):
        self.spnames = spnames
        self.nsp = len(spnames)
        self.lib = self.load_ceq_library()

        atoms = []
        M = []
        lewis = []
        for sp in spnames:
            asp, Msp, lsp = self.get_thermo_data(sp)
            atoms.append(asp)
            M.append(Msp)
            lewis.append(lsp)

        self.atoms = atoms
        self.M = array(M)
        self.lewis = array(lewis)

        elements = set()
        for a in atoms:
            for k in a.keys(): elements.add(k)
        elements = list(elements)
        elements.sort()

        self.elements = elements
        self.nel = len(elements)

        a = zeros((len(elements),len(spnames)))
        for i,s in enumerate(atoms):
            for k,v in s.items():
                j = elements.index(k)
                a[j,i] = v

        self.a = a

        return

    def get_thermo_data(self, name):
        """ Stripped down version of get_species from lewis_thermo.py """
        data = self.readdb(name)
        header = data[0].strip().split()
        if header[0]!=name:
            raise IOError("Database read failed! {}!={}".format(name, header[0]))

        info = data[1].strip().split()
        intervals = int(info[0])
        M = float(info[-2])/1e3 # convert to kg/mol from kg/kmol

        atomstrings = [data[1][10+8*i:10+8*(i+1)] for i in range(5)]
        atoms = {}
        for s in atomstrings:
            if s[0]==' ': continue
            elementname,amount = s.split()
            atoms[elementname] = float(amount)
        
        lewis = []
        for i in range(intervals):
            li = []
            line1 = data[3+3*i].replace('D','E') # Fix Fortran Double nonsense
            li.extend([float(line1[16*j:16*(j+1)]) for j in range(5)])

            line2 = data[4+3*i].replace('D','E') # Fix Fortran Double nonsense
            li.extend([float(line2[16*j:16*(j+1)]) for j in range(2)])
            li.extend([float(line2[16*j:16*(j+1)]) for j in range(3,5)])
            lewis.append(li)

        if len(lewis)!=3:
            lsp3 = self.fix_missing_thermo_segment(lewis, M)
            lewis.append(lsp3)
        assert len(lewis)==3 and len(lewis[0])==9 and len(lewis[1])==9 and len(lewis[2])==9
        return atoms, M, lewis

    def readdb(self,name):
        """ Retrieve species 'name' from the lewis_thermo.db file """
        with open(DBPATH) as fp:

            iterline=iter(fp) # Make an iterator for inner loop behaviour

            # Start looking for the beginning of a species entry
            for line in iterline:
                if not line[0] in letters: continue # Skip if not a letter in first position

                if line.startswith(name):        # This could be the one, check the rest
                    header = line.strip().split()
                    if header[0]!=name: continue # Nope false alarm

                    lines = [line]               # We've found it!
                    for nextline in iterline:
                        if nextline[0] in letters: break
                        lines.append(nextline)
                    break                        # Break the outer for loop and jump to return lines
            else:
                raise Exception("Name: {} not found!".format(name))
        return lines

    def fix_missing_thermo_segment(self, lsp, Mi):
        """ Create a NASA9 polynominal assuming constant Cp above 6000 K"""
        assert len(lsp)==2
        lspa = zeros((3,9))
        lspa[0] = lsp[0] # This should copy data automatically
        lspa[1] = lsp[1] # This should copy data automatically

        T = 5999.99999
        Ru = 8.3144621
        X = array([1.0])
        M = array([Mi])
        Td = c_double(T)

        c_double_p = POINTER(c_double)
        Xp    = X.ctypes.data_as(c_double_p)
        Mp    = M.ctypes.data_as(c_double_p)
        lewisp= lspa.ctypes.data_as(c_double_p)
        cp = self.lib.get_cp(Td, Xp, 1, lewisp, Mp)
        h  = self.lib.get_h(Td, Xp, 1, lewisp, Mp)
        s0 = self.lib.get_s0(Td, Xp, 1, lewisp, Mp)

        a2 = cp*Mi/Ru
        b1 = h*Mi/Ru - a2*T
        b2 = s0*Mi/Ru - a2*log(T)
        lsp2 = [0.0, 0.0, a2, 0.0, 0.0, 0.0, 0.0, b1, b2]
        return lsp2


    @staticmethod
    def load_ceq_library(LIBPATH=LIBPATH):
        """ Load the c library and set return types """
        lib = cdll.LoadLibrary(LIBPATH)
        lib.compute_Cp0_R.restype = c_double
        lib.pt.restype = c_int
        lib.rhou.restype = c_int
        lib.get_u.restype = c_double
        lib.get_h.restype = c_double
        lib.get_cp.restype = c_double
        lib.get_s0.restype = c_double
        lib.get_s.restype = c_double
        lib.batch_pt.restype = c_int
        lib.batch_rhou.restype = c_int
        lib.batch_u.restype = c_int
        return lib

    def pt(self, p, T, Xs0, verbose=0):
        """ Call c library to compute equilibrium concentrations at fixed p, T """
        if Xs0.size!=self.nsp: raise Exception('Mismatched array size {}!={}'.format(Xs0.size, self.nsp))
        Xs1 = zeros(Xs0.shape)
        pp = c_double(p)
        Tp = c_double(T)

        c_double_p = POINTER(c_double)
        Xs0p  = Xs0.ctypes.data_as(c_double_p)
        Mp    = self.M.ctypes.data_as(c_double_p)
        lewisp= self.lewis.ctypes.data_as(c_double_p)
        ap    = self.a.ctypes.data_as(c_double_p)
        Xs1p  = Xs1.ctypes.data_as(c_double_p)

        recode = self.lib.pt(pp, Tp, Xs0p, self.nsp, self.nel, lewisp, Mp, ap, Xs1p, verbose)
        if recode!=0: raise Exception("Equilibrium Calc Failed.")
        return Xs1

    def rhou(self, rho, u, Xs0, verbose=0):
        """ Call c library to compute equilibrium concentrations at fixed rho, u """
        if Xs0.size!=self.nsp: raise Exception('Mismatched array size {}!={}'.format(Xs0.size, self.nsp))
        Xs1 = zeros(Xs0.shape)
        rhop = c_double(rho)
        up = c_double(u)
        Tp = c_double()
        Tref = byref(Tp)

        c_double_p = POINTER(c_double)
        Xs0p  = Xs0.ctypes.data_as(c_double_p)
        Mp    = self.M.ctypes.data_as(c_double_p)
        lewisp= self.lewis.ctypes.data_as(c_double_p)
        ap    = self.a.ctypes.data_as(c_double_p)
        Xs1p  = Xs1.ctypes.data_as(c_double_p)

        recode = self.lib.rhou(rhop, up, Xs0p, self.nsp, self.nel, lewisp, Mp, ap, Xs1p, Tref, verbose)
        if recode!=0: raise Exception("Equilibrium Calc Failed.")
        T = Tp.value
        return Xs1, T

    def ps(self, pt, st, Xs0, verbose=0):
        """ Call c library to compute equilibrium concentrations at fixed p, s """
        if Xs0.size!=self.nsp: raise Exception('Mismatched array size {}!={}'.format(Xs0.size, self.nsp))
        Xs1 = zeros(Xs0.shape)
        pp = c_double(pt)
        sp = c_double(st)
        Tp = c_double()
        Tref = byref(Tp)

        c_double_p = POINTER(c_double)
        Xs0p  = Xs0.ctypes.data_as(c_double_p)
        Mp    = self.M.ctypes.data_as(c_double_p)
        lewisp= self.lewis.ctypes.data_as(c_double_p)
        ap    = self.a.ctypes.data_as(c_double_p)
        Xs1p  = Xs1.ctypes.data_as(c_double_p)

        recode = self.lib.ps(pp, sp, Xs0p, self.nsp, self.nel, lewisp, Mp, ap, Xs1p, Tref, verbose)
        if recode!=0: raise Exception("Equilibrium Calc Failed.")
        T = Tp.value
        return Xs1, T

    def get_u(self, X, T):
        """ Call c library to compute internal energy at fixed composition and temperature """
        Td = c_double(T)

        c_double_p = POINTER(c_double)
        Xp    = X.ctypes.data_as(c_double_p)
        Mp    = self.M.ctypes.data_as(c_double_p)
        lewisp= self.lewis.ctypes.data_as(c_double_p)

        u = self.lib.get_u(Td, Xp, self.nsp, lewisp, Mp)
        return u

    def get_h(self, X, T):
        """ Call c library to compute internal energy at fixed composition and temperature """
        Td = c_double(T)

        c_double_p = POINTER(c_double)
        Xp    = X.ctypes.data_as(c_double_p)
        Mp    = self.M.ctypes.data_as(c_double_p)
        lewisp= self.lewis.ctypes.data_as(c_double_p)

        h = self.lib.get_h(Td, Xp, self.nsp, lewisp, Mp)
        return h

    def get_cp(self, X, T):
        """ Call c library to compute internal energy at fixed composition and temperature """
        Td = c_double(T)

        c_double_p = POINTER(c_double)
        Xp    = X.ctypes.data_as(c_double_p)
        Mp    = self.M.ctypes.data_as(c_double_p)
        lewisp= self.lewis.ctypes.data_as(c_double_p)

        cp = self.lib.get_cp(Td, Xp, self.nsp, lewisp, Mp)
        return cp

    def get_s0(self, X, T):
        """ Call c library to compute specific entropy at standard state and arbitrary temperature """
        Td = c_double(T)

        c_double_p = POINTER(c_double)
        Xp    = X.ctypes.data_as(c_double_p)
        Mp    = self.M.ctypes.data_as(c_double_p)
        lewisp= self.lewis.ctypes.data_as(c_double_p)

        s0 = self.lib.get_s0(Td, Xp, self.nsp, lewisp, Mp)
        return s0

    def get_s(self, X, T, p):
        """ Call c library to compute internal entropy at an arbitrary pressure and temperature """
        Td = c_double(T)
        pd = c_double(p)

        c_double_p = POINTER(c_double)
        Xp    = X.ctypes.data_as(c_double_p)
        Mp    = self.M.ctypes.data_as(c_double_p)
        lewisp= self.lewis.ctypes.data_as(c_double_p)

        s = self.lib.get_s(Td, pd, Xp, self.nsp, lewisp, Mp)
        return s


    def batch_pt(self, p, T, Xs0, verbose=0):
        """ Call c library to compute equilibrium concentrations at fixed p, T """
        N, nspcheck = Xs0.shape
        if not Xs0.flags['OWNDATA']: raise Exception("Xs0 Memory Error: Array must own its data")
        if nspcheck!=self.nsp: raise Exception("nsp ({}) != Xs0.shape[1] ({})".format(self.nsp, nspcheck))
        if N!=p.size: raise Exception("p.size ({}) != Xs0.shape[0] ({})".format(p.size, N))
        if N!=T.size: raise Exception("T.size ({}) != Xs0.shape[0] ({})".format(T.size, N))

        Xs1 = zeros(Xs0.shape)

        c_double_p = POINTER(c_double)
        pp    = p.ctypes.data_as(c_double_p)
        Tp    = T.ctypes.data_as(c_double_p)
        Xs0p  = Xs0.ctypes.data_as(c_double_p)
        Mp    = self.M.ctypes.data_as(c_double_p)
        lewisp= self.lewis.ctypes.data_as(c_double_p)
        ap    = self.a.ctypes.data_as(c_double_p)
        Xs1p  = Xs1.ctypes.data_as(c_double_p)

        recode = self.lib.batch_pt(N, pp, Tp, Xs0p, self.nsp, self.nel, lewisp, Mp, ap, Xs1p, verbose)
        if recode!=0: raise Exception("Equilibrium Calc Failed.")
        return Xs1

    def batch_rhou(self, rho, u, Xs0, verbose=0):
        """ Call c library to compute equilibrium concentrations at fixed rho, u """
        N, nspcheck = Xs0.shape
        if not Xs0.flags['OWNDATA']: raise Exception("Xs0 Memory Error: Array must own its data")
        if nspcheck!=self.nsp: raise Exception("nsp ({}) != Xs0.shape[1] ({})".format(self.nsp, nspcheck))
        if N!=rho.size: raise Exception("rho.size ({}) != Xs0.shape[0] ({})".format(rho.size, N))
        if N!=u.size: raise Exception("u.size ({}) != Xs0.shape[0] ({})".format(u.size, N))

        Xs1 = zeros(Xs0.shape)
        T = zeros(rho.shape)

        c_double_p = POINTER(c_double)
        rhop  = rho.ctypes.data_as(c_double_p)
        up    = u.ctypes.data_as(c_double_p)
        Tp    = T.ctypes.data_as(c_double_p)
        Xs0p  = Xs0.ctypes.data_as(c_double_p)
        Mp    = self.M.ctypes.data_as(c_double_p)
        lewisp= self.lewis.ctypes.data_as(c_double_p)
        ap    = self.a.ctypes.data_as(c_double_p)
        Xs1p  = Xs1.ctypes.data_as(c_double_p)

        recode = self.lib.batch_rhou(N,rhop,up,Xs0p,self.nsp,self.nel,lewisp,Mp,ap,Xs1p,Tp,verbose)
        if recode!=0: raise Exception("Equilibrium Calc Failed.")
        return Xs1, T

    def batch_u(self, X, T):
        """ Call c library to compute internal energy at fixed composition and temperature """
        N, nspcheck = X.shape
        if not X.flags['OWNDATA']: raise Exception("X Memory Error: Array must own its data")
        if nspcheck!=self.nsp: raise Exception("nsp ({}) != X.shape[1] ({})".format(self.nsp, nspcheck))
        if N!=T.size: raise Exception("T.size ({}) != X.shape[0] ({})".format(T.size, N))

        u = zeros(T.shape)

        c_double_p = POINTER(c_double)
        Xp    = X.ctypes.data_as(c_double_p)
        Tp    = T.ctypes.data_as(c_double_p)
        up    = u.ctypes.data_as(c_double_p)
        Mp    = self.M.ctypes.data_as(c_double_p)
        lewisp= self.lewis.ctypes.data_as(c_double_p)

        recode = self.lib.batch_u(N, Tp, Xp, self.nsp, lewisp, Mp, up)
        if recode!=0: raise Exception("u calc failed.")
        return u

    def YtoX(self, Y):
        Mmix = 1.0/((Y/self.M).sum())
        X = Y*Mmix/self.M
        return X
    
    def XtoY(self, X):
        Mmix = (X*self.M).sum()
        Y = X*self.M/Mmix
        return Y

if __name__=='__main__':
    print("Called pyeq main!")
