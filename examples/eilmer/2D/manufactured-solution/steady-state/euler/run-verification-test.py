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

import os
import shutil
from math import log

simFiles = ['mms.lua',
            'constants.txt',
            'very-viscous-air.lua',
            'analytic_solution.py',
            'make_lua_files.py',
            'udf-source-template.lua',
            'udf-bc-template.lua',
            'ref-soln-template.lua',
            'fill-fn-template.lua',
	    'live-residuals.gplot']

tmpltDir = 'template'

def buildCaseStr(scale, ncells, fluxCalc, derivCalc, blocking):
    str = "%d\n" % scale
    str += "%s\n" % fluxCalc
    str += "%s\n" % derivCalc
    str += "%s\n" % xOrder
    str += "%s\n" % blocking
    str += "%d\n" % ncells
    return str

def buildRunStr(threading):
    str = "python make_lua_files.py\n"
    str += "e4shared --job=mms --prep\n"
    str += "e4zsss --job=mms"
    if threading == 'single':
        str += " --max-cpus=1"
    str += "\n"
    str += 'e4shared --job=mms --post --tindx-plot=last --ref-soln=ref-soln.lua  --norms="rho" | tail -3 > rho-norms.txt\n'
    return str

def prepareCases(scale, ncellsList, fluxCalc, derivCalc, blocking, threading):
    cwd = os.getcwd()
    for ncells in ncellsList:
        subDir = "%dx%d" % (ncells, ncells)
        try:
            os.mkdir(subDir)
        except OSError:
            pass
        os.chdir(subDir)
        f = open('case.txt', 'w')
        f.write(buildCaseStr(scale, ncells, fluxCalc, derivCalc, blocking))
        f.close()
        f = open('run.sh', 'w')
        f.write(buildRunStr(threading))
        f.close()
        for f in simFiles:
            fname = "../template/%s" % f
            shutil.copy(fname, ".")
        os.chdir(cwd)
    return

def runCases(ncellsList):
    cwd = os.getcwd()
    for ncells in ncellsList:
        print "========= Working on grid: %dx%d ===========" % (ncells, ncells)
        subDir = "%dx%d" % (ncells, ncells)
        os.chdir(subDir)
        cmd = "sh run.sh"
        os.system(cmd)
        os.chdir(cwd)
    return

def gatherResults(ncellsList):
    f0 = open('rho-L2-norm.dat', 'w')
    f1 = open('rho-L2-observed-order.dat', 'w')
    g0 = open('rho-Linf-norm.dat', 'w')
    g1 = open('rho-Linf-observed-order.dat', 'w')

    L2 = []; Linf = []; dx = []

    for i, ncells in enumerate(ncellsList):
        subDir = "%dx%d" % (ncells, ncells)
        dName = subDir + "/rho-norms.txt"
        f = open(dName, 'r')
        tks = f.readline().split()
        L2.append(float(tks[3]))
        Linf.append(float(tks[5]))
        f.close()
        dx.append(1.0/ncells)

    for i in range(len(ncellsList)):
        f0.write("%20.12e %20.12e\n" % (dx[i], L2[i]))
        g0.write("%20.12e %20.12e\n" % (dx[i], Linf[i]))
        if i == 0:
            continue
        logr = log(dx[i-1]/dx[i])
        p_L2 = log(L2[i-1]/L2[i])/logr
        p_Linf = log(Linf[i-1]/Linf[i])/logr
        f1.write("%20.12e %20.12e\n" % (dx[i-1], p_L2))
        g1.write("%20.12e %20.12e\n" % (dx[i-1], p_Linf))

    f0.close()
    f1.close()
    g0.close()
    g1.close()
    return


if __name__ == "__main__":
    import sys
    caseOptFile = sys.argv[1]
    execfile(caseOptFile)
    prepareCases(scale, ncellsList, fluxCalc, derivCalc, blocking, threading)
    runCases(ncellsList)
    gatherResults(ncellsList)
    
                
