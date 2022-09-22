#! /usr/bin/env python3
"""
Python program to post-process the simulation data from the Chicken 3D Flow Solver.

Usage:
  $ chkn-post --job=<jobName>

Author:
  P.A. Jacobs

Versions:
  2022-09-23  First Python code adpated from chkn_prep.py
"""

# ----------------------------------------------------------------------
#
import sys
import os
from getopt import getopt
import math
from copy import copy
import numpy as np
import json
from zipfile import ZipFile

from gdtk.geom.vector3 import Vector3, hexahedron_properties
from gdtk.geom.sgrid import StructuredGrid


shortOptions = "hf:"
longOptions = ["help", "job="]

def printUsage():
    print("Post-process a chicken run to produce VTK format files.")
    print("Usage: chkn-post" + \
          " [--help | -h]" + \
          " [--job=<jobName> | -f <jobName>]")
    print("")
    return

# --------------------------------------------------------------------

config = {}
times = {}
grids = {}
flows = {}

def read_config(jobDir):
    global config, times
    text = open(jobDir + '/config.json', 'r').read()
    config = json.loads(text)
    text = open(jobDir + '/times.data', 'r').readlines()
    for line in text:
        if len(line.strip()) == 0: continue
        if line[0] == '#': continue
        items = line.strip().split()
        times[int(items[0])] = float(items[1])
    return

def read_grids(jobDir):
    global config
    gridDir = jobDir+'/grid'
    if not os.path.exists(gridDir):
        raise RuntimeError('Cannot find grid directory: ' + gridDir)
    for k in range(config['nkb']):
        for j in range(config['njb']):
            for i in range(config['nib']):
                fileName = gridDir + ('/grid-%04d-%04d-%04d.gz' % (i, j, k))
                print("grid file:", fileName)
                if os.path.exists(fileName):
                    grids['%d,%d,%d'%(i,j,k)] = StructuredGrid(gzfile=fileName)
    return

# --------------------------------------------------------------------

if __name__ == '__main__':
    print("Begin chkn-post...")

    userOptions = getopt(sys.argv[1:], shortOptions, longOptions)
    uoDict = dict(userOptions[0])
    if len(userOptions[0]) == 0 or \
           "--help" in uoDict or \
           "-h" in uoDict:
        printUsage()
    else:
        if "--job" in uoDict:
            jobName = uoDict.get("--job", "")
        elif "-f" in uoDict:
            jobName = uoDict.get("-f", "")
        else:
            raise Exception("Job name is not specified.")
        jobDir, ext = os.path.splitext(jobName)
        #
        read_config(jobDir)
        print("times=", times)
        print("config=", config)
        read_grids(jobDir)
        print('grids=', grids)
        print("Done.")
    #
    sys.exit(0)
