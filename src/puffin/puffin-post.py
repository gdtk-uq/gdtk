#! /usr/bin/env python3
"""
Python program to reformat the flow data produced by the steady-state 2D flow solver.

Usage:
  $ puffin-post --job=<jobName> --output=<plotFormat>

Presently, only the legacy VTK format is implemented.
Eventually, we will want profiles and summaries.

PA Jacobs

2022-01-26: First code, now that we have some solutions to look at.
"""

# ----------------------------------------------------------------------
#
import sys
import os
from getopt import getopt
import math
import numpy as np
import json


shortOptions = "hf:o"
longOptions = ["help", "job=", "output="]

def printUsage():
    print("")
    print("Usage: puffin-prep" + \
          " [--help | -h]" + \
          " [--job=<jobName> | -f <jobName>]" + \
          " [--output=<plotFormat>] | -o <plotFormat>]")
    print("")
    return

# ----------------------------------------------------------------------
#
def writeVTKfiles(jobName):
    if not (os.path.exists(jobName) and os.path.isdir(jobName)):
        print("Could not find directory for job.")
        return
    startDir = os.getcwd()
    os.chdir(jobName)
    config = json.loads(open("config.json", "r").read())
    n_streams = config["n_streams"]
    ncells_all = [config["ncells_%d" % i] for i in range(n_streams)]
    print("n_streams=", n_streams)
    for i in range(n_streams):
        dataFile = open("flow-%d.data" % i, "r")
        txt = dataFile.readline()
        if txt.startswith('#'):
            variableNames = txt.strip('# ').split()
            print('variableNames=', variableNames)
            ncells = ncells_all[i]
            streamData = []
            while True:
                sliceData = []
                for j in range(ncells):
                    txt = dataFile.readline()
                    if txt:
                        items = [float(item) for item in txt.strip().split()]
                        if len(items) > 0: sliceData.append(items)
                    else:
                        break
                # At this point, we have the data for a full slice of cells
                # and/or we have arrived at the end of the file.
                if len(sliceData) == 0: break
                streamData.append(sliceData)
            # At this point, we have the full stream data.
            nxs = len(streamData)
            print("nxs=", nxs, "ncells=", ncells)
            plotFile = open("%s-flow-%d.vtk" % (jobName, i), "w")
            plotFile.write("# vtk DataFile Version 2.0\n")
            plotFile.write("# job: %s stream: %d\n" % (jobName, i))
            plotFile.write("ASCII\n")
            plotFile.write("DATASET STRUCTURED_GRID\n")
            plotFile.write("DIMENSIONS %d %d 1\n" % (ncells, nxs))
            plotFile.write("POINTS %d double\n" % (ncells*nxs))
            for ix in range(nxs):
                for j in range(ncells):
                    plotFile.write("%g %g 0.0\n" % (streamData[ix][j][0], streamData[ix][j][1]))
            plotFile.write("POINT_DATA %d\n" % (ncells*nxs))
            for iv in range(2, len(variableNames)):
                plotFile.write("SCALARS %s double 1\n" % (variableNames[iv],))
                plotFile.write("LOOKUP_TABLE default\n")
                for ix in range(nxs):
                    for j in range(ncells):
                        plotFile.write("%g\n" % (streamData[ix][j][iv],))
            plotFile.close()
        else:
            print("First line of stream flow file did not start with #")
        dataFile.close()
    os.chdir(startDir)
    return

# ----------------------------------------------------------------------
#
if __name__ == '__main__':
    print("Begin puffin-prep...")

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
        jobName, ext = os.path.splitext(jobName)
        print("jobName=", jobName)
        #
        if "--output" in uoDict:
            plotFormat = uoDict.get("--output", "")
        elif "-o" in uoDict:
            plotFormat = uoDict.get("-o", "")
        else:
            raise Exception("Plot format is not specified.")
        print("plotFormat=", plotFormat)
        #
        if (plotFormat == "vtk"): writeVTKfiles(jobName)
    print("Done.")
    sys.exit(0)
