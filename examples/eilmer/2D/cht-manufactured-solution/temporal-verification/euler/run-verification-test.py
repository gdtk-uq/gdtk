#!/usr/bin/env python
#
# This Python script is used to coordinate the
# execution of several grid refinements of
# a Method of Manufactured Solutions case such
# that an observed order of accuracy can be
# extracted. The user can configure the details
# of the case in a supplied Python file.
# For example, if the configuration file is
# called config.py, then this script may launched
# from the command line as:
#
# > python run-verification-test.py config.py
#
# Author: Rowan J. Gollan
# Date: Jun-2015
# Place: Tucson, Arizona
#
# Update: 2015-08-12
#         This version is for the coupled MMS case.
#         2017-04-15
#         All lua files based on analytic solution
#         are now auto-generated.

import os
import shutil
from math import log, sqrt

simFiles = ['coupled-mms.lua',
            'very-viscous-air.lua',
            'constants.txt',
            'analytic_solution.py',
            'make_lua_files.py',
            'udf-source-template.lua',
            'udf-solid-source-template.lua',
            'udf-bc-template.lua',
            'udf-solid-bc-template.lua',
            'fill-fn-template.lua',
            'fill-solid-fn-template.lua',
            'ref-soln-template.lua']

tmpltDir = 'template'

def buildCaseStr(dt):
    str = "%.12e\n" % dt
    return str

def buildRunStr():
    str = "python make_lua_files.py\n"
    str += "e4shared --job=coupled-mms --prep\n"
    str += "e4shared --job=coupled-mms --run --max-cpus=1\n"
    str += 'e4shared --job=coupled-mms --post --tindx-plot=last --ref-soln=ref-soln.lua  --norms="Ttr,T" --verbosity=0 | sed -n -e 6p -e 15p > T-norms.txt\n'
    return str

def prepareCases(dtList):
    cwd = os.getcwd()
    for i, dt in enumerate(dtList):
        subDir = "dt-%d" % i
        try:
            os.mkdir(subDir)
        except OSError:
            pass
        os.chdir(subDir)
        f = open('case.txt', 'w')
        f.write(buildCaseStr(dt))
        f.close()
        f = open('run.sh', 'w')
        f.write(buildRunStr())
        f.close()
        for f in simFiles:
            fname = "../template/%s" % f
            shutil.copy(fname, ".")
        os.chdir(cwd)
    return

def runCases(dtList):
    cwd = os.getcwd()
    for i, dt in enumerate(dtList):
        print "========= Working on case %d: dt= %.3e ===========" % (i, dt)
        subDir = "dt-%d" % i
        os.chdir(subDir)
        cmd = "sh run.sh"
        os.system(cmd)
        os.chdir(cwd)
    return

def gatherResults(ncellsList):
    f0 = open('T-L2-norm.dat', 'w')
    f1 = open('T-L2-observed-order.dat', 'w')
    g0 = open('T-Linf-norm.dat', 'w')
    g1 = open('T-Linf-observed-order.dat', 'w')

    L2 = []; Linf = []; dts = []

    for i, dt in enumerate(dtList):
        subDir = "dt-%d" % i
        dName = subDir + "/T-norms.txt"
        f = open(dName, 'r')
        tks = f.readline().split()
        L2G = float(tks[3])
        LinfG = float(tks[5])
        tks = f.readline().split()
        L2S = float(tks[3])
        LinfS = float(tks[5])
        f.close()
        dts.append(dt)
        L2.append(sqrt(L2G*L2G + L2S*L2S))
        Linf.append(max(LinfG, LinfS))

    for i in range(len(dtList)):
        f0.write("%20.12e %20.12e\n" % (dts[i], L2[i]))
        g0.write("%20.12e %20.12e\n" % (dts[i], Linf[i]))
        if i == 0:
            continue
        logr = log(dts[i-1]/dts[i])
        p_L2 = log(L2[i-1]/L2[i])/logr
        p_Linf = log(Linf[i-1]/Linf[i])/logr
        f1.write("%20.12e %20.12e\n" % (dts[i-1], p_L2))
        g1.write("%20.12e %20.12e\n" % (dts[i-1], p_Linf))

    f0.close()
    f1.close()
    g0.close()
    g1.close()
    return


if __name__ == "__main__":
    import sys
    caseOptFile = sys.argv[1]
    execfile(caseOptFile)
    prepareCases(dtList)
    runCases(dtList)
    gatherResults(dtList)
    


    
                
