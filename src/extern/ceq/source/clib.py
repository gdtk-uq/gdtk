"""
Automatic boilerplate for c library integration

"""
from ctypes import c_double,POINTER,c_int,byref,CDLL
from re import findall 
from collections import namedtuple
import numpy as np
from glob import glob

cfunction = namedtuple('cfunction', ['name','restype','argtypes','argnames'])
ndarray_d = np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS')
ndarray_i = np.ctypeslib.ndpointer(dtype=np.int32, flags='C_CONTIGUOUS')
ctypesmap = {'int':c_int, 'double':c_double, 'double*':ndarray_d, 'int*':ndarray_i}


class CLib(CDLL):
    def __init__(self, libname, headerfiles):
        super().__init__(libname)

        cfunctions = []
        for headerfile in headerfiles:
            with open(headerfile) as fp:
                text = fp.read() 

            matches = findall(r'extern (\w+) (\w+)\((.*)\);', text);
            for m in matches:
                argtypes,argnames=zip(*[i.strip().split() for i in m[2].split(',')])
                cfunctions.append(cfunction(m[1], m[0], argtypes, argnames))

        for function in cfunctions:
            getattr(self,function.name).restype = ctypesmap[function.restype]
            getattr(self,function.name).argtypes = [ctypesmap[i] for i in function.argtypes]

if __name__=='__main__':
    clib = CLib('libceq.so', glob('*.h')) 
