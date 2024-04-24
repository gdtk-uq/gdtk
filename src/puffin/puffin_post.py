#! /usr/bin/env python3
"""
Python program to reformat the flow data produced by the steady-state 2D flow solver.

PA Jacobs

2022-01-26: First code, now that we have some solutions to look at.
2022-02-04: Refactor to introduce stream line writer.
2024-02-18: Filled in cross-stream writer code.
"""

# ----------------------------------------------------------------------
#
import sys
import os
import time
from getopt import getopt
import math
import numpy as np
import json


shortOptions = "hf:o:s:c:x:"
longOptions = ["help", "job=", "output=", "stream-index=", "cell-index=", "x-location="]

def printUsage():
    print("")
    print("Usage: puffin-post [--help | -h] [--job=<jobName> | -f <jobName>]" + \
          " [--output=<plotFormat>] | -o <plotFormat>]")
    print("   where plotFormat is one of vtk|stream|cross")
    print("   For plotFormat==stream, also specify --stream-index=<i> | -s <i> and --cell-index=<j> | -c <j>")
    print("   For plotFormat==cross, also specify --x-location=<x> | -x <x>")
    print("")
    return

# ----------------------------------------------------------------------
#
def readStreamFlowFile(fileName, ncells):
    """
    Read the flow data for a single streamtube.
    """
    dataFile = open(fileName, "r")
    txt = dataFile.readline()
    if txt.startswith('#'):
        variableNames = txt.strip('# ').split()
        print('variableNames=', variableNames)
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
    else:
        print("First line of stream flow file did not start with #")
        streamData = None
        variableNames = None
    dataFile.close()
    return streamData, variableNames


def writeVTKfiles(jobName, n_streams, ncells_all):
    """
    Write a VTK-format file for the flow data in each stream tube.
    """
    for i in range(n_streams):
        ncells = ncells_all[i]
        streamData, variableNames = readStreamFlowFile("flow-%d.data" % i, ncells)
        if streamData:
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
    return


def writeStreamLineFile(jobName, streamIndex, ncells, cellIndex):
    """
    Write the data for a particular line of points along a stream tube.
    """
    ncells = ncells_all[streamIndex]
    streamData, variableNames = readStreamFlowFile("flow-%d.data" % streamIndex, ncells)
    if streamData:
        nxs = len(streamData)
        print("nxs=", nxs, "ncells=", ncells)
        dataFile = open("%s-streamline-%d-cell-%d.data" % (jobName, streamIndex, cellIndex), "w")
        dataFile.write("#")
        for name in variableNames: dataFile.write(" %s" % name)
        dataFile.write('\n')
        for ix in range(nxs):
            for iv in range(len(variableNames)):
                fmtStr = "%g" if iv == 0 else " %g"
                dataFile.write(fmtStr % streamData[ix][cellIndex][iv])
            dataFile.write('\n')
        dataFile.close()
    return


def writeCrossFile(jobName, n_streams, ncells_all, xLocation):
    """
    Write the GNUPlot data for a line of points across all stream tubes.
    """
    for i in range(n_streams):
        ncells = ncells_all[i]
        streamData, variableNames = readStreamFlowFile("flow-%d.data" % i, ncells)
        if streamData:
            nxs = len(streamData)
            print("nxs=", nxs, "ncells=", ncells)
            plotFile = open("%s-cross-%d.txt" % (jobName, i), "w")
            plotFile.write("# %s\n" % (' '.join(variableNames)))
            xFound = False
            for ix in range(nxs):
                x = streamData[ix][0][0]
                if x < xLocation: continue
                # Once we encounter the xLocation, write slice and stop.
                for j in range(ncells):
                    for iv in range(len(variableNames)):
                        fmt = '%g' if iv == 0 else ' %g'
                        plotFile.write(fmt % (streamData[ix][j][iv],))
                    plotFile.write("\n")
                xFound = True
                break
            if not xFound: print("Warning: did not find data for x=", xLocation)
            plotFile.close()
    return

# ----------------------------------------------------------------------
#
if __name__ == '__main__':
    print("Begin puffin postprocessing...")

    userOptions = getopt(sys.argv[1:], shortOptions, longOptions)
    uoDict = dict(userOptions[0])
    if len(userOptions[0]) == 0 or \
           "--help" in uoDict or \
           "-h" in uoDict:
        printUsage()
        sys.exit(0)
    #
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
    # Let's try to do some work in the job directory.
    #
    if not (os.path.exists(jobName) and os.path.isdir(jobName)):
        raise Exception("Could not find directory for job.")
    startDir = os.getcwd()
    os.chdir(jobName)
    config = json.loads(open("config.json", "r").read())
    n_streams = config["n_streams"]
    ncells_all = [config["ncells_%d" % i] for i in range(n_streams)]
    print("n_streams=", n_streams)
    if plotFormat == "vtk":
        writeVTKfiles(jobName, n_streams, ncells_all)
    elif plotFormat == "stream":
        if "--stream-index" in uoDict:
            streamIndex = int(uoDict.get("--stream-index", ""))
        elif "-s" in uoDict:
            streamIndex = int(uoDict.get("-s", ""))
        else:
            raise Exception("Stream index is not specified.")
        print("streamIndex=", streamIndex)
        if "--cell-index" in uoDict:
            cellIndexStr = uoDict.get("--cell-index", "")
        elif "-c" in uoDict:
            cellIndexStr = uoDict.get("-c", "")
        else:
            raise Exception("Cell index is not specified.")
        if (cellIndexStr == "$"):
            cellIndex = ncells_all[streamIndex] - 1
        else:
            cellIndex = int(cellIndexStr)
        print("cellIndex=", cellIndex)
        writeStreamLineFile(jobName, streamIndex, ncells_all[streamIndex], cellIndex)
    elif plotFormat == "cross":
        if "--x-location" in uoDict:
            xLocation = float(uoDict.get("--x-location", ""))
        elif "-x" in uoDict:
            xLocation = float(uoDict.get("-x", ""))
        else:
            raise Exception("x-location is not specified.")
        writeCrossFile(jobName, n_streams, ncells_all, xLocation)
    else:
        raise Exception("Not a valid plot format (vtk, stream, cross).")
    os.chdir(startDir)
    print("Done in {:.3f} seconds.".format(time.process_time()))
    sys.exit(0)
