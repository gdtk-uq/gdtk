#! /usr/bin/env python3
"""
Python program to post-process the simulation data from the Chicken 3D Flow Solver.

Usage:
  $ chkn-post --job=<jobName>

Author:
  P.A. Jacobs

Versions:
  2022-09-23  First Python code adpated from chkn_prep.py
  2022-11-09  Rebuilt flow file reader.
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
import gzip
import time
from collections import namedtuple

from gdtk.geom.vector3 import Vector3, hexahedron_properties
from gdtk.geom.sgrid import StructuredGrid


shortOptions = "hf:t:bvp:s:o:"
longOptions = ["help", "job=", "tindx=", "binary", "vtk-xml", "probe=", "slice=", "output-file="]

def printUsage():
    print("Post-process a chicken run to produce VTK format files.")
    print("Usage: chkn-post" +
          " [--help | -h]" +
          " [--job=<jobName> | -f <jobName>]" +
          " [--binary | -b]" +
          " [--tindx=<tindxSpec> | -t <tindxSpec>]" +
          " [--vtk-xml | -v]" +
          " [--slice=<sliceSpec> | -s <sliceSpec>]" +
          " [--output-file=<name> | -o <name>]" +
          " [--probe=<x,y,z> | -p <x,y,z>]"
    )
    print("  Default is to write VTK files. --vtk-xml")
    print("  You may elect to select slices of data, instead, to be written to a GNUPlot file.")
    print("  Or you may probe the data close to a specific point in space.")
    print("")
    print("  tindxSpec may specify a single index or a range of indices.")
    print("  Some examples: 0  1  $  -1  all  0:$  :$  :  0:-1  :-1")
    print("  If tindxSpec is not given, just the final-time snapshot is written.")
    print("")
    print("  When writing GNUPlot files of the flow data,")
    print("  sliceSpec may specify one or more slices of cells in particular blocks.")
    print("  See the code for more details.")
    print("  If the output file is not specified, it defaults to slice-tnnnn.data")
    print("")
    print("  For --probe, <x,y,z> represents 3 float numbers separated by commas.")
    print("  Output is just to the console.")
    print("")
    return

# --------------------------------------------------------------------

config = {}
times = {}
grids = {}
flows = {}

def read_config(jobDir):
    """
    Read the config and times files.
    """
    global config, times
    with open(jobDir + '/config.json', 'r') as fp:
        text = fp.read()
    config = json.loads(text)
    with open(jobDir + '/times.data', 'r') as fp:
        text = fp.readlines()
    for line in text:
        if len(line.strip()) == 0: continue
        if line[0] == '#': continue
        items = line.strip().split()
        times[int(items[0])] = float(items[1])
    return

def read_grids(jobDir, binaryData):
    """
    Read the full set of grids.
    """
    global config
    gridDir = jobDir+'/grid'
    if not os.path.exists(gridDir):
        raise RuntimeError('Cannot find grid directory: ' + gridDir)
    for k in range(config['nkb']):
        for j in range(config['njb']):
            for i in range(config['nib']):
                fileName = gridDir + ('/grid-%04d-%04d-%04d' % (i, j, k))
                if binaryData:
                    fileName += '.bin'
                    if os.path.exists(fileName):
                        grids['%d,%d,%d'%(i,j,k)] = StructuredGrid(binaryfile=fileName)
                else:
                    fileName += '.gz'
                    if os.path.exists(fileName):
                        grids['%d,%d,%d'%(i,j,k)] = StructuredGrid(gzfile=fileName)
    return

def read_flow_blocks(jobDir, tindx, binaryData):
    """
    Read the flow blocks for an individual tindx.
    """
    global config
    flowDir = jobDir + ('/flow/t%04d' % tindx)
    if not os.path.exists(flowDir):
        raise RuntimeError('Cannot find flow directory: ' + flowDir)
    for k in range(config['nkb']):
        for j in range(config['njb']):
            for i in range(config['nib']):
                fileName = flowDir + ('/flow-%04d-%04d-%04d' % (i, j, k))
                fileName += '.bin' if binaryData else '.gz'
                if os.path.exists(fileName):
                    flows['%d,%d,%d'%(i,j,k)] = read_block_of_flow_data(fileName, binaryData)
    return

def read_block_of_flow_data(fileName, binaryData):
    """
    The flow field data comes as a 1D array of float numbers.
    The data for each flow variable is appended end-to-end.

    Return the flow data in flattened arrays, one column for each flow-variable,
    because that is the arrangement that suits the VTK format files.
    """
    global config
    f = open(fileName, 'rb') if binaryData else gzip.open(fileName, 'rt')
    combinedData = np.frombuffer(f.read(), dtype=float) if binaryData else np.loadtxt(f)
    f.close()
    ntotal = combinedData.shape[0]
    nvars = len(config["iovar_names"])
    ncells = ntotal // nvars
    combinedData = combinedData.reshape((nvars,ncells))
    flowData = {}
    for j,var in enumerate(config["iovar_names"]):
        flowData[var] = combinedData[j,:]
    return flowData

def reshape_flow_data_arrays():
    """
    Reshape the individual arrays of data into the original 3D arrangement for each block.
    """
    global config, flows
    nics = config['nics']; njcs = config['njcs']; nkcs = config['nkcs']
    blk_ids = config['blk_ids']
    for ib in range(config['nib']):
        for jb in range(config['njb']):
            for kb in range(config['nkb']):
                if blk_ids[ib][jb][kb] >= 0:
                    # Block data should exist.
                    flowData = flows['%d,%d,%d'%(ib,jb,kb)]
                    for var in flowData.keys():
                        flowData[var] = flowData[var].reshape((nkcs[kb],njcs[jb],nics[ib])).transpose()
                    flows['%d,%d,%d'%(ib,jb,kb)] = flowData
    return

def find_nearest_cell(x, y, z):
    """
    Locate the nearest cell to the location (x,y,z), in any block.

    Returns the block and cell indices of the closest cell centre in a dictionary.
    """
    global config, flows
    closest = {'ib':-1, 'jb':-1, 'kb':-1, 'ic':-1, 'jc':-1, 'kc':-1, 'distance':sys.float_info.max}
    nics = config['nics']; njcs = config['njcs']; nkcs = config['nkcs']
    blk_ids = config['blk_ids']
    for ib in range(config['nib']):
        for jb in range(config['njb']):
            for kb in range(config['nkb']):
                if blk_ids[ib][jb][kb] >= 0:
                    # Block data should exist.
                    flowData = flows['%d,%d,%d'%(ib,jb,kb)]
                    for ic in range(nics[ib]):
                        for jc in range(njcs[jb]):
                            for kc in range(nkcs[kb]):
                                xc = flowData['posx'][ic][jc][kc]
                                yc = flowData['posy'][ic][jc][kc]
                                zc = flowData['posz'][ic][jc][kc]
                                distance = math.sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2)
                                if (distance < closest['distance']):
                                    closest['ib'] = ib; closest['jb'] = jb; closest['kb'] = kb
                                    closest['ic'] = ic; closest['jc'] = jc; closest['kc'] = kc
                                    closest['distance'] = distance
    return closest

# --------------------------------------------------------------------

def write_slice_list_to_gnuplot_file(sliceList, fileName):
    nics = config['nics']; njcs = config['njcs']; nkcs = config['nkcs']
    f = open(fileName, 'wt')
    f.write('# ')
    for var in config['iovar_names']: f.write(var+' ')
    f.write('\n')
    for slice in sliceList:
        ib = slice['ib']; jb = slice['jb']; kb = slice['kb']
        if config['blk_ids'][ib][jb][kb] < 0: continue
        # At this point, the block was active in the simulation and has valid flow data.
        flowData = flows['%d,%d,%d'%(ib, jb, kb)]
        if slice['dir'] == 'i':
            jc = slice['jc']; kc = slice['kc']
            for ic in range(nics[ib]):
                for var in config['iovar_names']:
                    f.write('%e ' % flowData[var][ic][jc][kc])
                f.write('\n')
        if slice['dir'] == 'j':
            ic = slice['ic']; kc = slice['kc']
            for jc in range(njcs[jb]):
                for var in config['iovar_names']:
                    f.write('%e ' % flowData[var][ic][jc][kc])
                f.write('\n')
        if slice['dir'] == 'k':
            ic = slice['ic']; jc = slice['jc']
            for kc in range(nkcs[kb]):
                for var in config['iovar_names']:
                    f.write('%e ' % flowData[var][ic][jc][kc])
                f.write('\n')
    f.close()
    return

# --------------------------------------------------------------------

def write_pvd_file(jobDir, tindxList, timesList):
    """
    Write a single .pvd file that ties all of the snapshot files together
    into a single time-stamped dataset that Paraview can read.
    """
    plotDir = jobDir + '/plot'
    if not os.path.exists(plotDir): os.mkdir(plotDir)
    fileName = plotDir + '/flow.pvd'
    fp = open(fileName, mode='w')
    fp.write('<?xml version="1.0"?>\n')
    fp.write('<VTKFile type="Collection" version="0.1" byte_order="BigEndian">\n')
    fp.write('<Collection>\n')
    for tindx,timeStamp in zip(tindxList, timesList):
        pvtsFileName = 'flow-t%04d.pvts' % (tindx,)
        fp.write('<DataSet timestep="%.18e" group="" part="0" file="%s"/>\n' % (timeStamp, pvtsFileName))
    fp.write('</Collection>\n')
    fp.write('</VTKFile>\n')
    fp.close()
    return

def write_vtk_files(jobDir, tindx):
    """
    Write snapshot of the flow data as a collection of VTK files for a single tindx.
    This collection of files will describe the pieces of an overall StructuredGrid.
    """
    global config, grids, flows
    plotDir = jobDir + '/plot'
    if not os.path.exists(plotDir): os.mkdir(plotDir)
    whole_niv = sum(config['nics']) + 1
    whole_njv = sum(config['njcs']) + 1
    whole_nkv = sum(config['nkcs']) + 1
    # The coordinating .pvts file is written as we write the individual .vts files.
    fileName = plotDir + ('/flow-t%04d.pvts' % (tindx,))
    fp = open(fileName, mode='w')
    fp.write('<VTKFile type="PStructuredGrid" version="0.1" byte_order="BigEndian">\n')
    fp.write('<PStructuredGrid WholeExtent="%d %d %d %d %d %d" GhostLevel="0">\n' %
             (0, whole_niv-1, 0, whole_njv-1, 0, whole_nkv-1))
    fp.write('<PCellData>\n')
    for var in config["iovar_names"]:
        fp.write('<PDataArray Name="%s" type="Float64" NumberOfComponents="1" format="ascii" />\n' % (var,))
    if "velx" in config["iovar_names"] and "a" in config["iovar_names"]:
        fp.write('<PDataArray Name="Mach" type="Float64" NumberOfComponents="1" format="ascii" />\n');
    if "velx" in config["iovar_names"]:
        fp.write('<PDataArray Name="vel.vector" type="Float64" NumberOfComponents="3" format="ascii" />\n');
    fp.write('</PCellData>\n')
    fp.write('<PPoints>\n')
    fp.write('<PDataArray type="Float64" NumberOfComponents="3" format="ascii" />\n')
    fp.write('</PPoints>\n')
    start_nkv = 0
    for k in range(config['nkb']):
        start_njv = 0
        for j in range(config['njb']):
            start_niv = 0
            for i in range(config['nib']):
                fileName = 'flow-t%04d-%04d-%04d-%04d.vts' % (tindx, i, j, k)
                key = '%d,%d,%d'%(i,j,k)
                if key in grids.keys():
                    # Write this piece only if the block exists.
                    grid = grids[key]
                    fp.write('<Piece Extent="%d %d %d %d %d %d" Source="%s" />\n' %
                             (start_niv, start_niv+grid.niv-1,
                              start_njv, start_njv+grid.njv-1,
                              start_nkv, start_nkv+grid.nkv-1, fileName))
                    write_vtk_structured_grid_file(plotDir+'/'+fileName, grid, flows[key],
                                                   whole_niv, whole_njv, whole_nkv,
                                                   start_niv, start_njv, start_nkv)
                start_niv += config['nics'][i]
            start_njv += config['njcs'][j]
        start_nkv += config['nkcs'][k]
    fp.write('</PStructuredGrid>\n')
    fp.write('</VTKFile>\n')
    fp.close()
    return

def write_vtk_structured_grid_file(fileName, grid, flowData,
                                   whole_niv, whole_njv, whole_nkv,
                                   start_niv, start_njv, start_nkv):
    """
    Combine the grid and flow data for one block into a VTK StructuredGrid file

    for one piece of the overall grid..
    """
    with open(fileName, mode='w') as fp:
        fp.write('<VTKFile type="StructuredGrid" version="0.1" byte_order="BigEndian">\n')
        fp.write('<StructuredGrid WholeExtent="%d %d %d %d %d %d">\n' %
                 (0, whole_niv-1, 0, whole_njv-1, 0, whole_nkv-1))
        fp.write('<Piece Extent="%d %d %d %d %d %d">\n' %
                 (start_niv, start_niv+grid.niv-1,
                  start_njv, start_njv+grid.njv-1,
                  start_nkv, start_nkv+grid.nkv-1))
        fp.write('<CellData>\n')
        for var in config["iovar_names"]:
            fp.write('<DataArray Name="%s" type="Float64" NumberOfComponents="1" format="ascii">\n' % (var,))
            for item in flowData[var]: fp.write('%g\n' % (item,))
            fp.write('</DataArray>\n')
        if "velx" in config["iovar_names"] and  "a" in config["iovar_names"]:
            fp.write('<DataArray Name="Mach" type="Float64" NumberOfComponents="1" format="ascii">\n');
            for a,vx,vy,vz in zip(flowData["a"],flowData["velx"],flowData["vely"],flowData["velz"]):
                fp.write('%g\n' % (math.sqrt(vx*vx+vy*vy+vz*vz)/a))
            fp.write('</DataArray>\n')
        if "velx" in config["iovar_names"]:
            fp.write('<DataArray Name="vel.vector" type="Float64" NumberOfComponents="3" format="ascii">\n');
            for velxyz in zip(flowData["velx"],flowData["vely"],flowData["velz"]): fp.write('%g %g %g\n' % velxyz)
            fp.write('</DataArray>\n')
        fp.write('</CellData>\n')
        fp.write('<Points>\n')
        fp.write('<DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        for k in range(grid.nkv):
            for j in range(grid.njv):
                for i in range(grid.niv):
                    x = grid.vertices.x[i][j][k]
                    y = grid.vertices.y[i][j][k]
                    z = grid.vertices.z[i][j][k]
                    fp.write('%g %g %g\n'%(x, y, z))
        fp.write('</DataArray>\n')
        fp.write('</Points>\n')
        fp.write('</Piece>\n')
        fp.write('</StructuredGrid>\n')
        fp.write('</VTKFile>\n')
    return

# --------------------------------------------------------------------

if __name__ == '__main__':
    print("Begin chkn-post...")

    userOptions = getopt(sys.argv[1:], shortOptions, longOptions)
    uoDict = dict(userOptions[0])
    if len(userOptions[0]) == 0 or ("--help" in uoDict) or ("-h" in uoDict):
        printUsage()
        sys.exit(0)
    #
    # Continue on, to do some real work.
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
    #
    binaryData = ("--binary" in uoDict) or ("-f" in uoDict)
    read_grids(jobDir, binaryData)
    #
    tindxSpec = "$" # default is the final-time index
    tindxList = []
    timesKeys = list(times)
    timesKeys.sort()
    if "--tindx" in uoDict:
        tindxSpec = uoDict.get("--tindx", "$")
    elif "-t" in uoDict:
        tindxSpec = uoDict.get("-t", "$")
    #
    # Decide which saved snapshots to write out, saving them to a list.
    if tindxSpec == '$' or tindxSpec == '-1':
        # Pick out the final-time index.
        tindxList = [timesKeys[-1]]
    elif tindxSpec == 'all':
        tindxList = timesKeys[:]
    elif tindxSpec.isnumeric():
        # A single integer is assumed.
        tindxList = [int(tindxSpec)]
    elif tindxSpec.find(':') >= 0:
        # We have a range specified
        firstSpec, lastSpec = tindxSpec.split(':')
        first = 0
        last = timesKeys[-1]
        #
        if firstSpec.isnumeric():
            first = int(firstSpec)
        elif firstSpec == '':
            first = 0
        elif firstSpec == '$' or firstSpec == '-1':
            first = timesKeys[-1]
        else:
            raise Exception("Cannot convert firstSpec={}".format(firstSpec))
        #
        if lastSpec.isnumeric():
            last = int(lastSpec)
        elif lastSpec == '':
            last = timesKeys[-1]
        elif lastSpec == '$' or lastSpec == '-1':
            last = timesKeys[-1]
        else:
            raise Exception("Cannot convert lastSpec={}".format(lastSpec))
        #
        tindxList = []
        for tindx in timesKeys:
            if tindx >= first and tindx <= last: tindxList.append(tindx)
    else:
        raise Exception("Did not know what to do with tindxSpec="+tindxSpec)
    #
    action = "vtk-xml" # Default action
    if ("--vtk-xml" in uoDict) or ("-v" in uoDict):
        action = "vtk-xml"
    elif ("--probe" in uoDict) or ("-p" in uoDict):
        action = "probe"
    elif ("--slice" in uoDict) or ("-s" in uoDict):
        action = "slice"
    #
    if action == "vtk-xml":
        print("Write out flow data snapshots as VTK files.")
        for tindx in tindxList:
            print("Writing tindx={}".format(tindx))
            read_flow_blocks(jobDir, tindx, binaryData)
            write_vtk_files(jobDir, tindx)
        #
        timesList = [times[tindx] for tindx in tindxList]
        write_pvd_file(jobDir, tindxList, timesList)
    #
    if action == "probe":
        if "--probe" in uoDict:
            xyzSpec = uoDict.get("--probe", "")
        elif "-p" in uoDict:
            xyzSpec = uoDict.get("-p", "")
        x, y, z = xyzSpec.strip().split(',')
        x = float(x); y = float(y); z = float(z)
        print("Probe the flow data close to ({}, {}, {}).".format(x,y,z))
        for tindx in tindxList:
            print("Probing flow data at tindx={}".format(tindx))
            read_flow_blocks(jobDir, tindx, binaryData)
            reshape_flow_data_arrays()
            closest = find_nearest_cell(x, y, z)
            print("closest cell at", closest)
            flowData = flows['%d,%d,%d'%(closest['ib'],closest['jb'],closest['kb'])]
            for var in config['iovar_names']:
                print('  {}: {}'.format(var, flowData[var][closest['ic'],closest['jc'],closest['kc']]))
    #
    if action == "slice":
        print("Select just a 1D of flow data.")
        if "--slice" in uoDict:
            sliceSpec = uoDict.get("--slice", "")
        elif "-s" in uoDict:
            sliceSpec = uoDict.get("-s", "")
        # Accept multiple slices, separated by semicolons.
        # Individual slices of the form ib,jb,kb,:,jc,kc or ib,jb,kb,ic,:,kc or ib,jb,kb,ic,jc,:
        # where the : is a literal colon and the other items are integers.
        # These select a slice of cells in one index direction (identified by the colon)
        # for the specific block.
        sliceSpecs = sliceSpec.strip().split(';')
        print(sliceSpecs)
        sliceList = []
        for spText in sliceSpecs:
            ib,jb,kb,ic,jc,kc = spText.split(',')
            slice = {'ib':int(ib), 'jb':int(jb), 'kb':int(kb)}
            if ic == ':':
                slice['dir'] = 'i';
                slice['ic'] = 0
                slice['jc'] = int(jc)
                slice['kc'] = int(kc)
            elif jc == ':':
                slice['dir'] = 'j';
                slice['ic'] = int(ic)
                slice['jc'] = 0
                slice['kc'] = int(kc)
            elif kc == ':':
                slice['dir'] = 'k';
                slice['ic'] = int(ic)
                slice['jc'] = int(jc)
                slice['kc'] = 0
            else:
                print("Did not specify a direction in the slice specification: "+spText)
                continue
            sliceList.append(slice)
        #
        if len(sliceList) == 0:
            print("No valid slices specified.")
        else:
            print("sliceList=", sliceList)
            filePrefix = "slice"
            if "--output-file" in uoDict:
                filePrefix = uoDict.get("--output-file", "")
            elif "-o" in uoDict:
                filePrefix = uoDict.get("-o", "")
            for tindx in tindxList:
                print("Slicing flow data at tindx={}".format(tindx))
                read_flow_blocks(jobDir, tindx, binaryData)
                reshape_flow_data_arrays()
                fileName = "%s-t%04d.data" % (filePrefix, tindx)
                write_slice_list_to_gnuplot_file(sliceList, fileName)
    #
    print("Done in {:.3f} seconds.".format(time.process_time()))
    sys.exit(0)
