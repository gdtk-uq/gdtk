"""
Python wrapper for nenzf1d.

@author: Nick Gibbons
"""

from cffi import FFI
from numpy import linspace 


ffi = FFI()
ffi.cdef("""
    int cwrap_run();
    typedef struct Foo Foo;
    struct Foo
    {
          int a;
          double b;
          int n_c;
          double* c;
          int n_name;
          char* name;
    };
    void pass_a_struct(Foo bar);
""")
C = ffi.dlopen('libnenzf1d.so')
C.cwrap_run()        # call run 

c = linspace(1.5,5.5,4)
print("c: ", c)

data = c.__array_interface__['data'][0]
cptr = ffi.cast ( "double*" , data )

fooStruct = ffi.new("struct Foo *")
fooStruct.a = 12;
fooStruct.b = 3.5;
fooStruct.n_c = c.size;
fooStruct.c = cptr;
cname = ffi.new("char[]", bytes("asdf", encoding="utf-8"))
fooStruct.n_name = 4
fooStruct.name = cname
C.pass_a_struct(fooStruct[0])
