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
#
# UPDATE: 25-Apr-2016
#         Script now generates source terms file
#                              BCs file
#                              reference solution file.
# UPDATE: 2018 (Kyle A. Damm)
#         Script now handles three-dimensional source terms

from sympy import *
from analytic_solution import *
import re
from sympy.utilities.codegen import codegen

# Thermodynamic behvaiour, equation of state and energy equation
e, T, et = symbols('e T et')
e = p/rho/(gamma-1)
T = e/Cv
et = e + u*u/2 + v*v/2 + w*w/2

# Heat flux
qx, qy, qz= symbols('qx qy qz')
qx = -k*diff(T, x)
qy = -k*diff(T, y)
qz = -k*diff(T, z)

# Shear stress
tauxx, tauyy, tauzz, tauxy, tauyz, tauzx= symbols('tauxx tauyy tauzz tauxy tauyz tauzx')
tauxx = 2./3*mu*(2*diff(u, x) - diff(v, y) - diff(w, z))
tauyy = 2./3*mu*(2*diff(v, y) - diff(u, x) - diff(w, z))
tauzz = 2./3*mu*(2*diff(w, z) - diff(u, x) - diff(v, y))
tauxy = mu*(diff(u, y) + diff(v, x))
tauyz = mu*(diff(w, y) + diff(v, z))
tauzx = mu*(diff(u, z) + diff(w, x))

# Navier-Stokes equations in conservative form
t, fmass, fxmom, fymom, fzmom, fe = symbols('t fmass fxmom fymom fzmom fe')
fmass = diff(rho, t) + diff(rho*u, x) + diff(rho*v, y) + diff(rho*w, z)
fxmom = diff(rho*u, t) + diff(rho*u*u+p-tauxx, x) + diff(rho*u*v-tauxy, y) + diff(rho*u*w-tauzx, z)
fymom = diff(rho*v, t) + diff(rho*v*u-tauxy, x) + diff(rho*v*v+p-tauyy, y) + diff(rho*v*w-tauyz, z)
fzmom = diff(rho*w, t) + diff(rho*w*u-tauzx, x) + diff(rho*w*v-tauyz, y) + diff(rho*w*w+p-tauzz, z)
fe = diff(rho*et, t) + diff(rho*u*et+p*u-u*tauxx-v*tauxy-w*tauzx+qx, x) + diff(rho*v*et+p*v-u*tauxy-v*tauyy-w*tauyz+qy, y) \
     + diff(rho*w*et+p*w-u*tauzx-v*tauyz-w*tauzz+qz, z)

def generateSource(varList):
    """This function generates Lua source code based on a selection
    of the expressions defined above. The selection is defined in the
    varList. The varList is a list of tuples of the form:

    varList = [("fmass", fmass), ("fxmom", fxmom)]

    This functions returns the generated source as a string.
    """
    [(f_name, f_code), (h_name, f_header)] = codegen(
        varList, "F95", "test", header=False)
    # Convert F95 to Lua code
    # This is heavily borrowed PJ's script: f90_to_lua.py
    # First we'll do some replacements
    f_code = f_code.replace('**', '^')
    f_code = f_code.replace('d0', '')
    f_code = f_code.replace('d-','e-')
    # Now we'll break into lines so that we can completely remove
    # some lines and tidy others
    lines = f_code.split('\n')
    lines[:] = [l.lstrip() for l in lines if ( not l.startswith('REAL*8') and
                                               not l.startswith('INTEGER*4') and
                                               not l.startswith('implicit') and
                                               not l.startswith('end')) ]
    # Now reassemble but collect the split lines into a large line
    buf = ""
    sourceCode = ""
    for i,l in enumerate(lines):
        if l.endswith('&'):
            buf = buf + l[:-2]
        else:
            if buf == "":
                sourceCode = sourceCode + l + '\n'
            else:
                sourceCode = sourceCode + buf + l + '\n'
                buf = ""
    return sourceCode

def createFileFromTemplate(sourceCode, templateName, fileName):
    """Given some source code and template file, do the text substitution
    and create the real file."""
    fin = open(templateName, 'r')
    templateText = fin.read()
    fin.close()
    luaText = templateText.replace('<insert-expressions-here>',
                                   sourceCode)

    fout = open(fileName, 'w')
    fout.write(luaText)
    fout.close()
    return

if __name__ == '__main__':
    taskList = [ {'fName': "udf-source-terms.lua", 'tName': "udf-source-template.lua",
                  'varList': [("fmass", fmass), ("fxmom", fxmom), ("fymom", fymom), ("fzmom", fzmom), ("fe", fe)]},
                 {'fName': "udf-bc.lua", 'tName': "udf-bc-template.lua",
                  'varList': [("tab.p", p), ("tab.T", T), ("tab.velx", u), ("tab.vely", v), ("tab.velz", w)]},
                 {'fName': "fill-fn.lua", 'tName': "fill-fn-template.lua",
                  'varList': [("p", p), ("T", T), ("velx", u), ("vely", v), ("velz", w)]},
                 {'fName': "ref-soln.lua", 'tName': "ref-soln-template.lua",
                  'varList': [("tab.rho", rho), ("tab.p", p), ("tab['T']", T),
                              ("tab['vel.x']", u), ("tab['vel.y']", v), ("tab['vel.z']", w)]},
        ]

    for task in taskList:
        sourceCode = generateSource(task['varList'])
        createFileFromTemplate(sourceCode, task['tName'], task['fName'])
    



