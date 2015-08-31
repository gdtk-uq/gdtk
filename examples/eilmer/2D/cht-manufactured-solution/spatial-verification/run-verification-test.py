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

import os
import shutil
from math import log

simFiles = ['coupled-mms.lua',
            'very-viscous-air.lua',
            'analytic_solution.py',
            'make_source_terms.py',
            'make_solid_source_terms.py',
            'udf-source-template.lua',
            'udf-solid-source-template.lua',
            'udf-bc.lua']

tmpltDir = 'template'

def buildCaseStr(ncells):
    str = "%d\n" % ncells
    return str

def buildRunStr():
    str = "python make_source_terms.py\n"
    str += "python make_solid_source_terms.py\n"
    str += "e4shared --job=coupled-mms --prep\n"
    str += "e4shared --job=coupled-mms --run --max-cpus=1\n"
    str += 'e4shared --job=coupled-mms --post --tindx-plot=40 --ref-soln=udf-bc.lua  --norms="T[0],T" | tail -10 > T-norms.txt\n'
    return str

def prepareCases(ncellsList):
    cwd = os.getcwd()
    for ncells in ncellsList:
        subDir = "%dx%d" % (ncells, 1.5*ncells)
        try:
            os.mkdir(subDir)
        except OSError:
            pass
        os.chdir(subDir)
        f = open('case.txt', 'w')
        f.write(buildCaseStr(ncells))
        f.close()
        f = open('run.sh', 'w')
        f.write(buildRunStr())
        f.close()
        for f in simFiles:
            fname = "../template/%s" % f
            shutil.copy(fname, ".")
        os.chdir(cwd)
    return

def runCases(ncellsList):
    cwd = os.getcwd()
    for ncells in ncellsList:
        print "========= Working on grid: %dx%d ===========" % (ncells, 1.5*ncells)
        subDir = "%dx%d" % (ncells, 1.5*ncells)
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

    L2 = []; Linf = []; dx = []

    for i, ncells in enumerate(ncellsList):
        subDir = "%dx%d" % (ncells, 1.5*ncells)
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
    prepareCases(ncellsList)
    runCases(ncellsList)
#    gatherResults(ncellsList)
    


    
                
