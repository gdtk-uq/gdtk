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
            'very-viscous-air.lua',
            'constants.txt',
            'make_lua_files.py',
            'analytic_solution.py',
            'udf-source-template.lua',
            'udf-bc-template.lua',
            'ref-soln-template.lua',
	    'fill-fn-template.lua',
	    'live-residuals.gplot']

tmpltDir = 'template'

def buildCaseStr(ncells, fluxCalc, derivCalc, derivLcn, blocking):
    str = "%s\n" % fluxCalc
    str += "%s\n" % derivCalc
    str += "%s\n" % derivLcn
    str += "2\n" # xOrder = 2
    str += "%s\n" % blocking
    str += "%d\n" % ncells
    return str

def buildRunStr(threading):
    str = "python make_lua_files.py\n"
    str += "e4shared --job=mms --prep\n"
    str += "e4-nk-shared --job=mms"
    if threading == 'single':
        str += " --max-cpus=1"
    str += "\n"
    str += 'e4shared --job=mms --post --tindx-plot=last --ref-soln=ref-soln.lua  --norms="rho,tke,omega" | sed -n -e 15p -e 19p -e 23p > norms.txt\n'
    return str

def prepareCases(ncellsList, fluxCalc, derivCalc, derivLcn, blocking, threading):
    cwd = os.getcwd()
    for ncells in ncellsList:
        subDir = "%dx%d" % (ncells, ncells)
        try:
            os.mkdir(subDir)
        except OSError:
            pass
        os.chdir(subDir)
        f = open('case.txt', 'w')
        f.write(buildCaseStr(ncells, fluxCalc, derivCalc, derivLcn, blocking))
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
    a0 = open('rho-L2-norm.dat', 'w')
    a1 = open('rho-L2-observed-order.dat', 'w')
    s0 = open('rho-Linf-norm.dat', 'w')
    s1 = open('rho-Linf-observed-order.dat', 'w')
    f0 = open('tke-L2-norm.dat', 'w')
    f1 = open('tke-L2-observed-order.dat', 'w')
    g0 = open('tke-Linf-norm.dat', 'w')
    g1 = open('tke-Linf-observed-order.dat', 'w')
    h0 = open('omega-L2-norm.dat', 'w')
    h1 = open('omega-L2-observed-order.dat', 'w')
    i0 = open('omega-Linf-norm.dat', 'w')
    i1 = open('omega-Linf-observed-order.dat', 'w')

    L2_rho = []; Linf_rho = [];
    L2_tke = []; Linf_tke = []; dx = []
    L2_omega = []; Linf_omega = []

    for i, ncells in enumerate(ncellsList):
        subDir = "%dx%d" % (ncells, ncells)
        dName = subDir + "/norms.txt"
        f = open(dName, 'r')
        tks = f.readline().split()
        L2_rho.append(float(tks[3]))
        Linf_rho.append(float(tks[5]))
        tks = f.readline().split()
        L2_tke.append(float(tks[3]))
        Linf_tke.append(float(tks[5]))
        tks = f.readline().split()
        L2_omega.append(float(tks[3]))
        Linf_omega.append(float(tks[5]))
        f.close()
        dx.append(1.0/ncells)

    for i in range(len(ncellsList)):
        a0.write("%20.12e %20.12e\n" % (dx[i], L2_rho[i]))
        s0.write("%20.12e %20.12e\n" % (dx[i], Linf_rho[i]))
        f0.write("%20.12e %20.12e\n" % (dx[i], L2_tke[i]))
        g0.write("%20.12e %20.12e\n" % (dx[i], Linf_tke[i]))
        h0.write("%20.12e %20.12e\n" % (dx[i], L2_omega[i]))
        i0.write("%20.12e %20.12e\n" % (dx[i], Linf_omega[i]))
        if i == 0:
            continue
        logr = log(dx[i-1]/dx[i])
        p_L2 = log(L2_rho[i-1]/L2_rho[i])/logr
        p_Linf = log(Linf_rho[i-1]/Linf_rho[i])/logr
        a1.write("%20.12e %20.12e\n" % (dx[i-1], p_L2))
        s1.write("%20.12e %20.12e\n" % (dx[i-1], p_Linf))
        p_L2 = log(L2_tke[i-1]/L2_tke[i])/logr
        p_Linf = log(Linf_tke[i-1]/Linf_tke[i])/logr
        f1.write("%20.12e %20.12e\n" % (dx[i-1], p_L2))
        g1.write("%20.12e %20.12e\n" % (dx[i-1], p_Linf))
        p_L2 = log(L2_omega[i-1]/L2_omega[i])/logr
        p_Linf = log(Linf_omega[i-1]/Linf_omega[i])/logr
        h1.write("%20.12e %20.12e\n" % (dx[i-1], p_L2))
        i1.write("%20.12e %20.12e\n" % (dx[i-1], p_Linf))

    a0.close()
    a1.close()
    s0.close()
    s1.close()
    f0.close()
    f1.close()
    g0.close()
    g1.close()
    h0.close()
    h1.close()
    i0.close()
    i1.close()

    return

if __name__ == "__main__":
    import sys
    caseOptFile = sys.argv[1]
    execfile(caseOptFile)
    prepareCases(ncellsList, fluxCalc, derivCalc, derivLcn, blocking, threading)
    runCases(ncellsList)
    gatherResults(ncellsList)
    


    
                
