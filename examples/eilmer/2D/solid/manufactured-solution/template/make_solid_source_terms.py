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

rho_s = 10000
Cp_s = 100.0

e_s, fe_s, t = symbols('e_s fe_s t')
# What should go here is the LHS of the energy equation
# for the solid phase
e_s = rho_s*Cp_s*T_s
fe_s = diff(e_s, t) - diff(k_s*diff(T_s, x), x) - diff(k_s*diff(T_s, y), y)

if __name__ == '__main__':
    import re
    from sympy.utilities.codegen import codegen
    print 'Generating manufactured source terms.'
    [(f_name, f_code), (h_name, f_header)] = codegen(
        [("fe_s", fe_s)], "F95", "test", header=False)
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

    fin = open('udf-solid-source-template.lua', 'r')
    template_text = fin.read()
    fin.close()
    lua_text = template_text.replace('<insert-source-terms-here>',
                                     f_code)

    fout = open('udf-solid-source-terms.lua', 'w')
    fout.write(lua_text)
    fout.close()
    print 'Done converting to Lua.'

    



