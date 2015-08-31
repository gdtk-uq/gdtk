# Author: Rowan J. Gollan
# Place: The University of Queensland, Brisbane, Australia
# Date: 06-Jun-2014
#
# This script is used to generate the analytical source
# terms required to run the Method of Manufactured Solutions
# test case. The generated code is in Fortran95 and it can
# be converted to Lua with a separate script.
#
# This is an exercise in using sympy to generate the source
# terms. It is a transliteration of PJ's original work
# done using Maxima.

from sympy import *
from analytic_solution import *

Rgas, g, Prandtl, Cv, Cp = symbols('Rgas g Prandtl Cv Cp')
Rgas = 287.0
g = 1.4
Prandtl = 1.0
Cv = Rgas/(g-1)
Cp = g*Cv

mu, k = symbols('mu k')
k = 10000.0
mu = k*Prandtl/Cp

# Thermodynamic behaviour, equation of state and energy equation
e, p, et = symbols('e p et')
p = rho*Rgas*T
e = Cv*T
et = e + u*u/2 + v*v/2

# Heat flux
qx, qy = symbols('qx qy')
qx = -k*diff(T, x)
qy = -k*diff(T, y)

# Shear stress
tauxx, tauyy, tauxy = symbols('tauxx tauyy tauxy')
tauxx = 2./3*mu*(2*diff(u, x) - diff(v, y))
tauyy = 2./3*mu*(2*diff(v, y) - diff(u, x))
tauxy = mu*(diff(u, y) + diff(v, x))

# Navier-Stokes equations in conservative form
t, fmass, fxmom, fymom, fe = symbols('t fmass fxmom fymom fe')
fmass = diff(rho, t) + diff(rho*u, x) + diff(rho*v, y)
fxmom = diff(rho*u, t) + diff(rho*u*u+p-tauxx, x) + diff(rho*u*v-tauxy, y)
fymom = diff(rho*v, t) + diff(rho*v*u-tauxy, x) + diff(rho*v*v+p-tauyy, y)
fe = diff(rho*et, t) + diff(rho*u*et+p*u-u*tauxx-v*tauxy+qx, x) + diff(rho*v*et+p*v-u*tauxy-v*tauyy+qy, y)

if __name__ == '__main__':
    import re
    from sympy.utilities.codegen import codegen
    print 'Generating manufactured source terms.'
    [(f_name, f_code), (h_name, f_header)] = codegen(
        [("fmass", fmass), ("fxmom", fxmom), ("fymom", fymom), ("fe", fe)],
        "F95", "test", header=False)
    # Convert F95 to Lua code
    # This is heavily borrowed PJ's script: f90_to_lua.py
    # First we'll do some replacements
    f_code = f_code.replace('**', '^')
    f_code = f_code.replace('d0', '')
    # Now we'll break into lines so that we can completely remove
    # some lines and tidy others
    lines = f_code.split('\n')
    lines[:] = [l.lstrip() for l in lines if ( not l.startswith('REAL*8') and
                                               not l.startswith('implicit') and
                                               not l.startswith('end')) ]
    # Now reassemble but collect the split lines into a large line
    buf = ""
    f_code = ""
    for i,l in enumerate(lines):
        if l.endswith('&'):
            buf = buf + l[:-2]
        else:
            if buf == "":
                f_code = f_code + l + '\n'
            else:
                f_code = f_code + buf + l + '\n'
                buf = ""

    fin = open('udf-source-template.lua', 'r')
    template_text = fin.read()
    fin.close()
    lua_text = template_text.replace('<insert-source-terms-here>',
                                     f_code)

    fout = open('udf-source-terms.lua', 'w')
    fout.write(lua_text)
    fout.close()
    print 'Done converting to Lua.'

    



