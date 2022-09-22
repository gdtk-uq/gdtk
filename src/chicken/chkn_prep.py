#! /usr/bin/env python3
"""
Python program to set up a simulation for the Chicken 3D Flow Solver.

It is intended for the user to define their particular facility and
flow in terms of the data objects defined in this module.
As part of its initialization, this program will read and execute
a user-specified input script that contains, in Python,
the code that defines the facility geometry and gas-path details.

Usage:
  $ chkn_prep --job=<jobName>

The simulation control data is then organised via the classes
GlobalConfig, FlowState, BoundaryCondition and its subclasses.
These classes provide places to store the configuration information
and their function/method names appear as commands in the user's
job description file.

Note that config is a reference to the global configuration information
describing the simulation.  There is one such variable set up by the
main program and the user's script should directly set the attributes
of this variable to adjust settings for the simulation.

Author:
  P.A. Jacobs

Versions:
  2022-09-19  First Python code adpated from l1d4-prep.py
"""

# ----------------------------------------------------------------------
#
import sys
import os
from getopt import getopt
from abc import ABC, abstractmethod
import collections
import math
from copy import copy
import numpy as np
import json
from zipfile import ZipFile

from gdtk.geom.vector3 import Vector3, hexahedron_properties
from gdtk.geom.volume import TFIVolume
from gdtk.geom.sgrid import StructuredGrid
from gdtk.numeric import roberts

sys.path.append("") # so that we can find user's scripts in current directory

shortOptions = "hf:"
longOptions = ["help", "job="]

def printUsage():
    print("Prepare a chicken run by reading a job script and writing grid and flow files.")
    print("Usage: chkn-prep" + \
          " [--help | -h]" + \
          " [--job=<jobName> | -f <jobName>]")
    print("")
    return

#----------------------------------------------------------------------
# We will work with a simple ideal gas.
# Presently, this just for inviscid flow. We have yet to add viscous effects.

class IdealGas():
    """
    A simple gas model.
    """
    __slots__ = 'gamma', 'R', 'Cv'

    def __init__(self, gamma=1.4, R=287.1):
        self.gamma = gamma
        self.R = R
        self.Cv = R/(gamma-1.0)

    def to_json(self):
        result = '{"gamma": %g, "R": %g, "Cv": %g}' % (self.gamma, self.R, self.Cv)
        return result

#----------------------------------------------------------------------
# This is where we store the core data for the simulation.

class GlobalConfig():
    """
    Contains the simulation control parameters.

    The user's script should not create one of these
    but should specify the simulation parameters by
    altering the attributes of the global object "config"
    that already exists by the time the user's script executes.
    """
    count = 0

    # We want to prevent the user's script from introducing new attributes
    # via typographical errors.
    __slots__ = 'job_name', 'title', 'gas_model', \
        'nib', 'njb', 'nkb', 'blk_ids', 'nics', 'njcs', 'nkcs', \
        'dt_init', 'cfl_list', 'cfl_count', 'print_count', \
        'dt_plot_list', 'max_time', 'max_step', \
        'x_order', 't_order', 'flow_var_names'

    def __init__(self):
        """Accepts user-specified data and sets defaults. Make one only."""
        if GlobalConfig.count >= 1:
            raise Exception("Already have a GlobalConfig object defined.")

        self.job_name = ""
        self.title = "Chicken Simulation."
        self.gas_model = IdealGas()
        self.nib = None
        self.njb = None
        self.nkb = None
        self.dt_init = 1.0e-6
        self.cfl_list = []
        self.cfl_count = 10;
        self.print_count = 50;
        # If dt_plot_list is still an empty list when we write
        # the parameter file, just use the max_time value for
        # both dt_plot and dt_his.  Fill in later, when necessary.
        self.dt_plot_list = []
        self.max_time = 1.0e-3
        self.max_step = 10
        self.x_order = 2
        self.t_order = 2
        # The following is not meant to be edited for individual simulations.
        self.flow_var_names = ['pos.x', 'pos.y', 'pos.z', 'vol',
                               'p', 'T', 'rho', 'e', 'a',
                               'vel.x', 'vel.y', 'vel.z']
        #
        GlobalConfig.count += 1
        return

    def write(self, fp):
        """
        Writes the configuration data to the specified file in JSON format.
        """
        global config, fluidBlocksList, flowStateList
        #
        fp.write('  "title": "%s",\n' % self.title)
        fp.write('  "gas_model": %s,\n' % self.gas_model.to_json())
        fp.write('  "nib": %d,\n' % self.nib)
        fp.write('  "njb": %d,\n' % self.njb)
        fp.write('  "nkb": %d,\n' % self.nkb)
        fp.write('  "max_time": %e,\n' % self.max_time)
        fp.write('  "max_step": %d,\n' % self.max_step)
        fp.write('  "dt_init": %e,\n' % self.dt_init)
        if len(self.cfl_list) == 0: add_cfl_value(0.0, 0.5) # Add a default value.
        n_cfl = len(self.cfl_list)
        assert n_cfl > 0, "Require at least one cfl_value to be specified."
        fp.write('  "n_cfl": %d,\n' % n_cfl)
        tlist = [self.cfl_list[i][0] for i in range(n_cfl)]
        fp.write('  "cfl_times": %s,\n' % json.dumps(tlist))
        tlist = [self.cfl_list[i][1] for i in range(n_cfl)]
        assert min(tlist) > 0.0, "Require positive nonzero cfl values."
        fp.write('  "cfl_values": %s,\n' % json.dumps(tlist))
        fp.write('  "cfl_count": %d,\n' % self.cfl_count)
        fp.write('  "print_count": %d,\n' % self.print_count)
        fp.write('  "x_order": %d,\n' % self.x_order)
        fp.write('  "t_order": %d,\n' % self.t_order)
        #
        if len(self.dt_plot_list) == 0:
            # Since the user did not specify any, default to the end.
            add_dt_plot(0.0, config.max_time, config.max_time)
        n_dt_plot = len(self.dt_plot_list)
        assert n_dt_plot > 0, "Require at least one dt_plot to be specified."
        fp.write('  "n_dt_plot": %d,\n' % n_dt_plot)
        tlist = [self.dt_plot_list[i][0] for i in range(n_dt_plot)]
        fp.write('  "t_change": %s,\n' % json.dumps(tlist))
        tlist = [self.dt_plot_list[i][1] for i in range(n_dt_plot)]
        assert min(tlist) > 0.0, "Require positive nonzero dt_plot values."
        fp.write('  "dt_plot": %s,\n' % json.dumps(tlist))
        tlist = [self.dt_plot_list[i][2] for i in range(n_dt_plot)]
        assert min(tlist) > 0.0, "Require positive nonzero dt_hist values."
        fp.write('  "dt_hist": %s,\n' % json.dumps(tlist))
        fp.write('  "flow_var_names": %s,\n' % json.dumps(self.flow_var_names))
        return

    def check_array_of_fluid_blocks(self):
        """
        Allocates the individual blocks to locations in the global array of blocks.

        The array of blocks need not have a block at every location but the number
        of cells for all blocks in a line must be consistent, so that simple full-face
        copying can be done to exchange the flow data.
        """
        global fluidBlocksList
        #
        # First, work out the shape of the overall array of blocks.
        imax = 0; jmax = 0; kmax = 0
        for b in fluidBlocksList:
            imax = max(b.i, imax)
            jmax = max(b.j, jmax)
            kmax = max(b.k, kmax)
        self.nib = imax+1
        self.njb = jmax+1
        self.nkb = kmax+1
        self.blk_ids = np.zeros((self.nib, self.njb, self.nkb), dtype=int)
        self.blk_ids -= 1 # Start with a negative id so that we can see missing blocks.
        #
        # Identify the array positions for the fluid blocks that are present and
        # check that their numbers of cells are consistent across the array.
        self.nics = np.zeros(self.nib, dtype=int)
        self.njcs = np.zeros(self.njb, dtype=int)
        self.nkcs = np.zeros(self.nkb, dtype=int)
        for b in fluidBlocksList:
            i = b.i; j = b.j; k = b.k
            if self.blk_ids[i,j,k] >= 0:
                raise RuntimeError('Already a block at i=%d, j=%d, k=%d' % (i, j, k))
            self.blk_ids[i,j,k] = b.indx
            nic = b.grid.niv-1
            njc = b.grid.njv-1
            nkc = b.grid.nkv-1
            if self.nics[i] == 0: self.nics[i] = nic
            if self.njcs[j] == 0: self.njcs[j] = njc
            if self.nkcs[k] == 0: self.nkcs[k] = nkc
            if nic != self.nics[i]:
                raise RuntimeError('Block at i=%d j=%d k=%d has inconsistent nic=%d' % (i, j, k, nic))
            if njc != self.njcs[j]:
                raise RuntimeError('Block at i=%d j=%d k=%d has inconsistent njc=%d' % (i, j, k, njc))
            if nkc != self.nkcs[k]:
                raise RuntimeError('Block at i=%d j=%d k=%d has inconsistent nkc=%d' % (i, j, k, nkc))
        #
        return

# We will create just one GlobalConfig object that the user can alter.
config = GlobalConfig()

# --------------------------------------------------------------------
# The following functions are to provide convenient ways of setting
# some of the GlobalConfig elements.

def add_cfl_value(t_change, cfl_value):
    """
    Add a (t,cfl) tuple to the cfl tuple list in GlobalConfig.

    t_change: (float) The time, in seconds, at which this cfl should take effect.
    cfl_value: (float) Ratio of actual time step divided by allowable time step.
    Returns: the new length of the list
    """
    global config
    if len(config.cfl_list) > 0:
        # Check that we are adding points monotonically in x.
        if t_change > config.cfl_list[-1][0]:
            config.cfl_list.append((t_change, cfl_value))
        else:
            print("Warning: did not add t_change,cfl_value tuple (", t_change, cfl_value, ").")
    else:
        config.cfl_list.append((t_change, cfl_value))
    return len(config.cfl_list)


def add_dt_plot(t_change, dt_plot, dt_his):
    """
    Add a dt tuple to the dt_plot tuple list in GlobalConfig.

    t_change: (float) The time, in seconds,
        at which this dt_plot and dt_his should take effect.
    dt_plot: (float) Time interval between writing whole solutions
        for later plotting.
    dt_his: (float) Time interval between writing data to history file.
    Returns: the new length of the list
    """
    global config
    if len(config.dt_plot_list) > 0:
        # Check that we are adding points monotonically in x.
        if t_change > config.dt_plot_list[-1][0]:
            config.dt_plot_list.append((t_change, dt_plot, dt_his))
        else:
            print("Warning: did not add dt_plot tuple (", \
                  t_change, dt_plot, dt_his, ").")
    else:
        config.dt_plot_list.append((t_change, dt_plot, dt_his))
    return len(config.dt_plot_list)


#----------------------------------------------------------------------
# The following classes define the objects that will be assembled into
# the simulation.
# ---------------------------------------------------------------------
#
# We will accumulate references to defined objects.
gasStatesList = []
flowStatesList = []
fluidBlocksList = []

class GasState():
    """
    The data that defines the state of the gas.
    """
    __slots__ = 'indx', 'p', 'T', 'rho', 'e', 'a'

    def __init__(self, p, T):
        """
        Input pressure (in Pa) and temperature (in degrees K).
        """
        global config, gasStatesList
        self.p = p
        self.T = T
        self.update_from_pT(config.gas_model)
        self.indx = len(gasStatesList)
        gasStatesList.append(self)

    def update_from_pT(self, gmodel):
        """
        Given p and T, update the other thermodynamic properties
        for a particular gas model.
        """
        self.e = gmodel.Cv * self.T
        self.rho = self.p / (gmodel.R * self.T)
        self.a = math.sqrt(gmodel.gamma * gmodel.R * self.T)
        return

    def __repr__(self):
        return "GasState(p={}, T={})".format(self.p, self.T)

    def to_json(self):
        result = '{"p": %g, "T": %g, "rho": %g, "e": %g, "a": %g}' \
            % (self.p, self.T, self.rho, self.e, self.a)
        return result


class FlowState():
    """
    Contains the gas properties and velocities defining a flow state.

    The user may create more than one FlowState to describe the initial
    gas properties and boundary conditions throughout the flow domain.
    """

    __slots__ = 'indx', 'label', 'gas', 'vel'

    def __init__(self, gas=None, p=100.0e3, T=300.0,
                 vel=None, velx=0.0, vely=0.0, velz=0.0,
                 label=""):
        global flowStatesList
        if isinstance(gas, GasState):
            self.gas = copy(gas)
        else:
            self.gas = GasState(p=p, T=T)
        if isinstance(vel, Vector3):
            self.vel = copy(vel)
        else:
            self.vel = Vector3(velx, vely, velz)
        self.label = label
        self.indx = len(flowStatesList)
        flowStatesList.append(self)

    def __repr__(self):
        return "FlowState(gas={}, vel={})".format(self.gas, self.vel)

    def to_json(self):
        result = '{"gas": %s, "vel": [%g, %g, %g]}' \
            % (self.gas.to_json(), self.vel.x, self.vel.y, self.vel.z)
        return result

# There is a bit of boiler place code here but it gives us
# a nice notation for use in the user's input script.

FaceNames = 'iminus iplus jminus jplus kminus kplus'
Face = collections.namedtuple('_', FaceNames)(*range(6))
FaceList = FaceNames.split()

BCNames = 'exchange wall_slip wall_no_slip inflow outflow'
BCCode = collections.namedtuple('_', BCNames)(*range(5))

class BoundaryCondition(ABC):
    """
    Base class for the family of boundary conditions.
    """
    __slots__ = ['tag']

    @abstractmethod
    def __repr__(self):
        pass

class ExchangeBC(BoundaryCondition):
    """
    Block-face-to-block-face-exchange boundary condition.
    """
    __slots__ = ['tag']

    def __init__(self):
        self.tag = 'exchange'

    def __repr__(self):
        return "ExchangeBC()"

    def to_json(self):
        return '{"tag": "%s"}' % (self.tag)

class WallWithSlipBC(BoundaryCondition):
    """
    Wall with slip boundary condition.
    """
    __slots__ = ['tag']

    def __init__(self):
        self.tag = 'wall_slip'

    def __repr__(self):
        return "WallWithSlipBC()"

    def to_json(self):
        return '{"tag": "%s"}' % (self.tag)

class WallNoSlipBC(BoundaryCondition):
    """
    Wall with slip boundary condition.
    """
    __slots__ = ['tag']

    def __init__(self):
        self.tag = 'wall_no_slip'

    def __repr__(self):
        return "WallNoSlipBC()"

    def to_json(self):
        return '{"tag": "%s"}' % (self.tag)

class InflowBC(BoundaryCondition):
    """
    Inflow boundary condition.
    """
    __slots__ = ['tag', 'fs']

    def __init__(self, fs):
        self.tag = 'inflow'
        self.fs = copy(fs)

    def __repr__(self):
        return "InflowBC(fs={})".format(self.fs)

    def to_json(self):
        return '{"tag": "%s", "flow_state": "%s"}' % (self.tag, self.fs)

class OutflowBC(BoundaryCondition):
    """
    Outflow boundary condition.
    """
    __slots__ = ['tag']

    def __init__(self):
        self.tag = 'outflow'

    def __repr__(self):
        return "OutflowBC()"

    def to_json(self):
        return '{"tag": "%s"}' % (self.tag)


class FluidBlock():
    """
    This class is used to build a description of the flow domain.
    """
    __slots__ = 'indx', 'i', 'j', 'k', 'grid', 'initialState', \
        'nic', 'njc', 'nkc', 'bcs', 'label', \
        'cellc', 'cellv'

    def __init__(self, i=0, j=0, k=0, grid=None, initialState=None,
                 bcs={}, label=""):
        global config, fluidBlocksList, flowStatesList
        if isinstance(grid, StructuredGrid):
            self.grid = grid
        else:
            raise RuntimeError('Need to supply a StructuredGrid object.')
        self.nic = grid.niv - 1
        self.njc = grid.njv - 1
        self.nkc = grid.nkv - 1
        if isinstance(initialState, FlowState):
            self.initialState = initialState
        else:
            raise RuntimeError('Need to supply a FlowState object.')
        #
        # Boundary conditions
        # Set default values and then overwrite, if the user has supplied them.
        self.bcs = {}
        for f in FaceList: self.bcs[f] = WallWithSlipBC()
        for key in bcs.keys():
            if key in [Face.iminus, 'iminus']: self.bcs['iminus'] = bcs[key]
            if key in [Face.iplus, 'iplus']: self.bcs['iplus'] = bcs[key]
            if key in [Face.jminus, 'jminus']: self.bcs['jminus'] = bcs[key]
            if key in [Face.iplus, 'jplus']: self.bcs['jplus'] = bcs[key]
            if key in [Face.kminus, 'kminus']: self.bcs['kminus'] = bcs[key]
            if key in [Face.kplus, 'kplus']: self.bcs['kplus'] = bcs[key]
        #
        self.i = i
        self.j = j
        self.k = k
        self.label = label
        self.indx = len(fluidBlocksList)
        fluidBlocksList.append(self)

    def to_json(self):
        """
        A representation of the block's configuration.
        """
        result = '{"i": %d, "j": %d, "k": %d, "initial_flow_state": %d,' % \
            (self.i, self.j, self.k, self.initialState.indx)
        result += ' "bcs": {'
        for i,f in enumerate(FaceList):
            result += '"%s": %s' % (f, self.bcs[f].to_json())
            result += ',' if i < len(FaceList)-1 else ''
        result += '},'
        result += ' "label": "%s"}' % self.label
        return result

    def construct_cell_properties(self):
        """
        We want the cell volumes and centroids available for writing
        into the initial flow file.
        """
        self.cellc = np.zeros((self.nic, self.njc, self.nkc), dtype=Vector3)
        self.cellv = np.zeros((self.nic, self.njc, self.nkc), dtype=float)
        for i in range(0, self.nic):
            for j in range(0, self.njc):
                for k in range(0, self.nkc):
                    # Construct the cell properties.
                    p000 = self.grid.vertices[i][j][k]
                    p100 = self.grid.vertices[i+1][j][k]
                    p110 = self.grid.vertices[i+1][j+1][k]
                    p010 = self.grid.vertices[i][j+1][k]
                    p001 = self.grid.vertices[i][j][k+1]
                    p101 = self.grid.vertices[i+1][j][k+1]
                    p111 = self.grid.vertices[i+1][j+1][k+1]
                    p011 = self.grid.vertices[i][j+1][k+1]
                    centroid, volume = hexahedron_properties(p000, p100, p110, p010,
                                                             p001, p101, p111, p011)
                    self.cellc[i][j][k] = centroid
                    self.cellv[i][j][k] = volume
        return

    def write_flow_data_to_zip_file(self, fileName, varNamesList):
        """
        The format is approximately Eilmer's new IO format.
        """
        fs = self.initialState
        gas = fs.gas
        vel = fs.vel
        with ZipFile(fileName, mode='w') as zf:
            for varName in varNamesList:
                with zf.open(varName, mode='w') as fp:
                    for k in range(0, self.nkc):
                        for j in range(0, self.njc):
                            for i in range(0, self.nic):
                                value = 0.0
                                if varName == 'pos.x': value = self.cellc[i][j][k].x
                                elif varName == 'pos.y': value = self.cellc[i][j][k].y
                                elif varName == 'pos.z': value = self.cellc[i][j][k].z
                                elif varName == 'vol': value = self.cellv[i][j][k]
                                elif varName == 'p': value = gas.p
                                elif varName == 'T': value = gas.T
                                elif varName == 'rho': value = gas.rho
                                elif varName == 'e': value = gas.e
                                elif varName == 'a': value = gas.a
                                elif varName == 'vel.x': value = vel.x
                                elif varName == 'vel.y': value = vel.y
                                elif varName == 'vel.z': value = vel.z
                                else: raise RuntimeError('unhandled variable name: ' + varName)
                                fp.write(f'{value}\n'.encode('utf-8'))
            # end varName loop
        return

# --------------------------------------------------------------------

def write_initial_files():
    """
    Writes the files needed for the main simulation code.

    These files are found in the directory config.job_name.
    """
    global config, fluidBlocksList, flowStatesList
    print("Begin write initial files.")
    #
    if not os.path.exists(config.job_name):
        os.mkdir(config.job_name)
    #
    with open(config.job_name+'/config.json', 'w') as fp:
        fp.write("{\n")
        config.write(fp)
        #
        n_flow_states = len(flowStatesList)
        fp.write('  "n_flow_states": %d,\n' % n_flow_states)
        fp.write('  "flow_states": [\n')
        for i,fs in enumerate(flowStatesList):
            fp.write('    %s' % fs.to_json())
            fp.write(',\n' if i<(n_flow_states-1) else '\n')
        fp.write('  ],\n')
        #
        n_fluid_blocks = len(fluidBlocksList)
        fp.write('  "n_fluid_blocks": %d,\n' % n_fluid_blocks)
        fp.write('  "blk_ids": %s,\n' % json.dumps(config.blk_ids.tolist()))
        fp.write('  "nics": %s,\n' % json.dumps(config.nics.tolist()))
        fp.write('  "njcs": %s,\n' % json.dumps(config.njcs.tolist()))
        fp.write('  "nkcs": %s,\n' % json.dumps(config.nkcs.tolist()))
        fp.write('  "fluid_blocks": [\n')
        for i,fb in enumerate(fluidBlocksList):
            fp.write('    %s' % fb.to_json())
            fp.write(',\n' if i<(n_fluid_blocks-1) else '\n')
        fp.write('  ]\n') # no comma for last item in JSON dict
        #
        fp.write("}\n")
    #
    print('Write the grid files for individual blocks.')
    #
    gridDir = config.job_name+'/grid'
    if not os.path.exists(gridDir):
        os.mkdir(gridDir)
    for fb in fluidBlocksList:
        fileName = gridDir + ('/grid-%04d-%04d-%04d.gz' % (fb.i, fb.j, fb.k))
        fb.grid.write_to_gzip_file(fileName)
    #
    print('Write the initial flow-field files.')
    #
    flowDir = config.job_name+'/flow/t0000'
    if not os.path.exists(flowDir):
        os.makedirs(flowDir)
    for fb in fluidBlocksList:
        fb.construct_cell_properties()
        fileName = flowDir + ('/flow-%04d-%04d-%04d.zip' % (fb.i, fb.j, fb.k))
        fb.write_flow_data_to_zip_file(fileName, config.flow_var_names)
    #
    with open(config.job_name+'/times.data', 'w') as fp:
        fp.write("# tindx t\n")
        fp.write("0 0.0\n")
    #
    print("End write initial files.")
    return

# --------------------------------------------------------------------

if __name__ == '__main__':
    print("Begin chkn-prep...")

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
        config.job_name, ext = os.path.splitext(jobName)
        inputScriptName = config.job_name + ".py"
        print("Input script file name: %s" % inputScriptName)

        # The user-specified input comes in the form of Python code.
        # It is up to the user to be careful; there is no security.
        exec(compile(open(inputScriptName, "rb").read(), inputScriptName, 'exec'))
        print("Summary:")
        print("  flow states       :", len(flowStatesList))
        print("  fluid blocks      :", len(fluidBlocksList))
        if len(fluidBlocksList) < 1:
            print("Warning: no fluid blocks defined; this is unusual.")
        config.check_array_of_fluid_blocks()
        write_initial_files()
    print("Done.")
    sys.exit(0)
