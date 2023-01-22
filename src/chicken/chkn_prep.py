#! /usr/bin/env python3
"""
Python program to set up a simulation for the Chicken 3D Flow Solver.

It is intended for the user to define their particular facility and
flow in terms of the data objects defined in this module.
As part of its initialization, this program will read and execute
a user-specified input script that contains, in Python,
the code that defines the facility geometry and gas-path details.

Usage:
  $ chkn-prep --job=<jobName> --binary

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
import gzip
import time

from gdtk.geom.vector3 import Vector3, hexahedron_properties
from gdtk.geom.volume import TFIVolume
from gdtk.geom.sgrid import StructuredGrid
from gdtk.numeric import roberts

sys.path.append("") # so that we can find user's scripts in current directory

shortOptions = "hf:b"
longOptions = ["help", "job=", "binary"]

def printUsage():
    print("Prepare a chicken run by reading a job script and writing grid and flow files.")
    print("Usage: chkn-prep" + \
          " [--help | -h]" + \
          " [--job=<jobName> | -f <jobName>]" + \
          " [--binary | -b]")
    print("")
    return

#----------------------------------------------------------------------

class GasModel():
    """
    A simple ideal gas model with a possibility of heat release.

    The selection of the gas model properties needs to be consistent with selection made in gas.cu.
    That selection is made at compile time, so remember to rebuild the simulation code if you
    change parameter values.
    """
    __slots__ = 'gamma', 'R', 'Cv', 'q', 'alpha', 'Ti'

    def __init__(self, gamma=1.4, R=287.1, q=0.0, alpha=0.0, Ti=0.0):
        self.gamma = gamma
        self.R = R
        self.Cv = R/(gamma-1.0)
        self.q = q
        self.alpha = alpha
        self.Ti = Ti

    def to_json(self):
        result = '{"gamma": %g, "R": %g, "Cv": %g, "q": %g, "alpha": %g, "Ti": %g}' % \
            (self.gamma, self.R, self.Cv, self.q, self.alpha, self.Ti)
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
        'x_order', 't_order', 'flux_calc', 'viscous', 'reacting', \
        'source_terms', 'iovar_names', 'threads_per_gpu_block'

    def __init__(self):
        """Accepts user-specified data and sets defaults. Make one only."""
        if GlobalConfig.count >= 1:
            raise Exception("Already have a GlobalConfig object defined.")

        self.job_name = ""
        self.title = "Chicken Simulation."
        self.gas_model = GasModel()
        self.nib = None
        self.njb = None
        self.nkb = None
        self.dt_init = 1.0e-6
        self.cfl_list = []
        self.cfl_count = 10;
        self.print_count = 20;
        # If dt_plot_list is still an empty list when we write
        # the parameter file, just use the max_time value for
        # both dt_plot and dt_his.  Fill in later, when necessary.
        self.dt_plot_list = []
        self.max_time = 1.0e-3
        self.max_step = 100
        self.x_order = 2
        self.t_order = 2
        self.flux_calc = "ausmdv"
        self.reacting = False
        self.viscous = False
        self.source_terms = "none"
        self.threads_per_gpu_block = 128
        # The following is not meant to be edited for individual simulations but
        # should be kept consistent with the symbols in IOvar namespace
        # that is defined in cell.cu.
        self.iovar_names = ['posx', 'posy', 'posz', 'vol',
                            'p', 'T', 'rho', 'e', 'YB', 'a',
                            'velx', 'vely', 'velz']
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
        fp.write('  "flux_calc": "%s",\n' % self.flux_calc)
        fp.write('  "reacting": %s,\n' % json.dumps(self.reacting))
        fp.write('  "viscous": %s,\n' % json.dumps(self.viscous))
        fp.write('  "source_terms": "%s",\n' % self.source_terms)
        fp.write('  "threads_per_gpu_block": %d,\n' % self.threads_per_gpu_block)
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
        fp.write('  "iovar_names": %s,\n' % json.dumps(self.iovar_names))
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
        # As each block is defined in the input script, it will get a location in the array.
        # We will use the following array of ids to indicate which blocks in the array
        # have been defined.  We will allow some blocks to be left undefined.
        self.blk_ids = np.zeros((self.nib, self.njb, self.nkb), dtype=int)
        self.blk_ids -= 1 # Start with a negative id so that we can see undefined blocks.
        #
        # Identify the array positions for the fluid blocks that are defined and
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
                raise RuntimeError('FluidBlock at i=%d j=%d k=%d has inconsistent nic=%d' % (i, j, k, nic))
            if njc != self.njcs[j]:
                raise RuntimeError('FluidBlock at i=%d j=%d k=%d has inconsistent njc=%d' % (i, j, k, njc))
            if nkc != self.nkcs[k]:
                raise RuntimeError('FluidBlock at i=%d j=%d k=%d has inconsistent nkc=%d' % (i, j, k, nkc))
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


def add_dt_plot(t_change, dt_plot, dt_his=1.0e-6):
    """
    Add a dt tuple to the dt_plot tuple list in GlobalConfig.

    t_change: (float) The time, in seconds,
        at which this dt_plot and dt_his should take effect.
    dt_plot: (float) Time interval between writing whole solutions
        for later plotting.
    dt_his: (float) Optional time interval between writing data to history file.
        2022-10-08 We are not presently writing history files. This may change.

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
flowStatesList = []
fluidBlocksList = []

class GasState():
    """
    The data that defines the state of the gas.
    """
    __slots__ = 'indx', 'p', 'T', 'rho', 'e', 'YB', 'a'

    def __init__(self, p, T, YB=0.0):
        """
        Input pressure (in Pa) and temperature (in degrees K).
        """
        self.p = p
        self.T = T
        self.YB = YB
        self.update_from_pT(config.gas_model)

    def update_from_pT(self, gmodel):
        """
        Given p and T, update the other thermodynamic properties
        for a particular gas model.
        """
        self.e = gmodel.Cv * self.T + self.YB * gmodel.q
        self.rho = self.p / (gmodel.R * self.T)
        self.a = np.sqrt(gmodel.gamma * gmodel.R * self.T)
        return

    def __repr__(self):
        return "GasState(p={}, T={}, YB={})".format(self.p, self.T, self.YB)

    def to_json(self):
        result = '{"p": %g, "T": %g, "rho": %g, "e": %g, "YB": %g, "a": %g}' \
            % (self.p, self.T, self.rho, self.e, self.YB, self.a)
        return result


class FlowState():
    """
    Contains the gas properties and velocities defining a flow state.

    The user may create more than one FlowState to describe the initial
    gas properties and boundary conditions throughout the flow domain.
    """

    __slots__ = 'indx', 'label', 'gas', 'vel'

    def __init__(self, gas=None, p=100.0e3, T=300.0, YB=0.0,
                 vel=None, velx=0.0, vely=0.0, velz=0.0,
                 label=""):
        if isinstance(gas, GasState):
            self.gas = copy(gas)
        else:
            self.gas = GasState(p=p, T=T, YB=YB)
        if isinstance(vel, Vector3):
            self.vel = copy(vel)
        else:
            self.vel = Vector3(velx, vely, velz)
        self.label = label

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
        self.tag = 'wall_with_slip'

    def __repr__(self):
        return "WallWithSlipBC()"

    def to_json(self):
        return '{"tag": "%s"}' % (self.tag)

class WallNoSlipAdiabaticBC(BoundaryCondition):
    """
    Wall, no-slip boundary condition, no heat flux.
    """
    __slots__ = ['tag']

    def __init__(self):
        self.tag = 'wall_no_slip_adiabatic'

    def __repr__(self):
        return "WallNoSlipAdiabaticBC()"

    def to_json(self):
        return '{"tag": "%s"}' % (self.tag)

class WallNoSlipFixedTBC(BoundaryCondition):
    """
    Wall, no-slip boundary condition, fixed temperature.

    TWall can be either single value or a user-supplied function TWall(x, y, z).
    If the form is a function, specific values for each boundary face will be
    written to a file.
    """
    __slots__ = ['tag', 'form', 'value', 'fun', 'values']

    def __init__(self, TWall):
        self.tag = 'wall_no_slip_fixed_T'
        self.fun = None
        self.value = 0.0
        if callable(TWall):
            self.form = 'fun'
            self.fun = TWall
            # Note that will fill in the specific face values later,
            # in the context of a given grid.
        elif isinstance(TWall, (int, float)):
            self.form = 'value'
            self.value = float(TWall)
        else:
            raise RuntimeError("TWall should be a single number or a function of (x,y,z)")

    def __repr__(self):
        s = "WallNoSlipFixedTBC(TWall="
        if self.form == 'fun':
            s += 'user-supplied-function)'
        else:
            s += '%g)' % (value)
        return s

    def to_json(self):
        return '{"tag": "%s", "form": "%s", "value": %g}' % (self.tag, self.form, self.value)

    def write_values_to_file(self, fileName):
        with gzip.open(fileName, 'wt') as f:
            for v in self.values: f.write('%g\n' % v)
        return

    def generate_values(self, grid, face):
        """
        Generate the temperature value for each face.
        Note that the loop nesting order must match the reader in the C++ world.
        """
        if not isinstance(grid, StructuredGrid):
            raise RuntimeError('Did not get a StructuredGrid object.')
        nic = grid.niv - 1
        njc = grid.njv - 1
        nkc = grid.nkv - 1
        self.values = []
        if face in ['iminus','iplus']:
            i = 0 if face == 'iminus' else nic
            p0 = Vector3(x=grid.vertices.x[i, :-1, :-1],
                         y=grid.vertices.y[i, :-1, :-1],
                         z=grid.vertices.z[i, :-1, :-1])
            p1 = Vector3(x=grid.vertices.x[i, 1:, :-1],
                         y=grid.vertices.y[i, 1:, :-1],
                         z=grid.vertices.z[i, 1:, :-1])
            p2 = Vector3(x=grid.vertices.x[i, 1:, 1:],
                         y=grid.vertices.y[i, 1:, 1:],
                         z=grid.vertices.z[i: 1:, 1:])
            p3 = Vector3(x=grid.vertices.x[i, :-1, 1:],
                         y=grid.vertices.y[i, :-1, 1:],
                         z=grid.vertices.z[i, :-1, 1:])
            pmid = 0.25*(p0+p1+p2+p3)
            for k in range(nkc):
                for j in range(njc):
                    self.values.append(self.fun(pmid.x[j,k], pmid.y[j,k], pmid.z[j,k]))
        elif face in ['jminus','jplus']:
            j = 0 if face == 'jminus' else njc
            p1 = Vector3(x=grid.vertices.x[:-1, j, :-1],
                         y=grid.vertices.y[:-1, j, :-1],
                         z=grid.vertices.z[:-1, j, :-1])
            p0 = Vector3(x=grid.vertices.x[1:, j, :-1],
                         y=grid.vertices.y[1:, j, :-1],
                         z=grid.vertices.z[1:, j, :-1])
            p3 = Vector3(x=grid.vertices.x[1:, j, 1:],
                         y=grid.vertices.y[1:, j, 1:],
                         z=grid.vertices.z[1:, j, 1:])
            p2 = Vector3(x=grid.vertices.x[:-1, j, 1:],
                         y=grid.vertices.y[:-1, j, 1:],
                         z=grid.vertices.z[:-1, j, 1:])
            pmid = 0.25*(p0+p1+p2+p3)
            for k in range(nkc):
                for i in range(nic):
                    self.values.append(self.fun(pmid.x[i,k], pmid.y[i,k], pmid.z[i,k]))
        elif face in ['kminus','kplus']:
            k = 0 if face == 'kminus' else nkc
            p0 = Vector3(x=grid.vertices.x[:-1, :-1, k],
                         y=grid.vertices.y[:-1, :-1, k],
                         z=grid.vertices.z[:-1, :-1, k])
            p1 = Vector3(x=grid.vertices.x[1:, :-1, k],
                         y=grid.vertices.y[1:, :-1, k],
                         z=grid.vertices.z[1:, :-1, k])
            p2 = Vector3(x=grid.vertices.x[:-1, 1:, k],
                         y=grid.vertices.y[:-1, 1:, k],
                         z=grid.vertices.z[:-1, 1:, k])
            p3 = Vector3(x=grid.vertices.x[1:, 1:, k],
                         y=grid.vertices.y[1:, 1:, k],
                         z=grid.vertices.z[1:, 1:, k])
            pmid = 0.25*(p0+p1+p2+p3)
            for j in range(njc):
                for i in range(nic):
                    self.values.append(self.fun(pmid.x[i,j], pmid.y[i,j], pmid.z[i,j]))
        else:
            raise RuntimeError("Unknown face: "+str(face))
        return


class InflowBC(BoundaryCondition):
    """
    Inflow boundary condition.
    """
    __slots__ = ['tag', 'fs', 'fsi']

    def __init__(self, fs):
        self.tag = 'inflow'
        if type(fs) is FlowState:
            self.fs = fs
        else:
            raise RuntimeError("Inflow boundary condition expects a FlowState object.")

    def __repr__(self):
        return "InflowBC(fs={})".format(self.fs)

    def to_json(self):
        # Note that the index to the FlowState will be assigned when
        # a Block is constructed that has InflowBC objects.
        return '{"tag": "%s", "flow_state_index": %d}' % (self.tag, self.fsi)

class InflowFunctionBC(BoundaryCondition):
    """
    Inflow boundary condition based on a function coded into the main program.
    """
    __slots__ = ['tag', 'fun_name']

    def __init__(self, fun_name):
        self.tag = 'inflow_function'
        self.fun_name = fun_name

    def __repr__(self):
        return "InflowFunctionBC(fun_name={})".format(self.fun_name)

    def to_json(self):
        return '{"tag": "%s", "fun_name": "%s"}' % (self.tag, self.fun_name)

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
    __slots__ = 'indx', 'active', 'i', 'j', 'k', 'grid', 'initialState', \
        'nic', 'njc', 'nkc', 'bcs', 'label', \
        'cellc', 'cellv', 'flow'

    def __init__(self, i=0, j=0, k=0, grid=None, initialState=None,
                 bcs={}, active=True, label=""):
        # When constructing a FluidBlock, we fill out the full grid and
        # the flow data in each cell within the block.
        # This allow us to use numpy arrays, for speed.
        global config, fluidBlocksList, flowStatesList
        if isinstance(grid, StructuredGrid):
            self.grid = grid
        else:
            raise RuntimeError('Need to supply a StructuredGrid object.')
        self.construct_cell_properties()
        #
        self.nic = grid.niv - 1
        self.njc = grid.njv - 1
        self.nkc = grid.nkv - 1
        assert(self.nic==self.cellc.x.shape[0])
        assert(self.njc==self.cellc.y.shape[1])
        assert(self.nkc==self.cellc.z.shape[2])
        ncells = self.nic*self.njc*self.nkc
        #
        if isinstance(initialState, FlowState):
            self.flow = {'p': np.full(ncells, initialState.gas.p), 'T':np.full(ncells, initialState.gas.T),
                         'rho':np.full(ncells, initialState.gas.rho), 'e':np.full(ncells, initialState.gas.e),
                         'YB':np.full(ncells, initialState.gas.YB), 'a':np.full(ncells, initialState.gas.a),
                         'velx': np.full(ncells, initialState.vel.x),
                         'vely': np.full(ncells, initialState.vel.y),
                         'velz': np.full(ncells, initialState.vel.z)}
        elif callable(initialState):
            # The user has supplied a function,
            # which we use one cell-centre at a time.
            self.flow = {'p': np.zeros(ncells, dtype=float), 'T': np.zeros(ncells, dtype=float),
                         'rho': np.zeros(ncells, dtype=float), 'e': np.zeros(ncells, dtype=float),
                         'YB': np.zeros(ncells, dtype=float), 'a': np.zeros(ncells, dtype=float),
                         'velx': np.zeros(ncells, dtype=float),
                         'vely': np.zeros(ncells, dtype=float),
                         'velz': np.zeros(ncells, dtype=float)}
            idx = 0
            # Note VTK ordering of loops
            for kk in range(self.nkc):
                for jj in range(self.njc):
                    for ii in range(self.nic):
                        x = self.cellc.x[ii,jj,kk]; y = self.cellc.y[ii,jj,kk]; z = self.cellc.z[ii,jj,kk]
                        state = initialState(x,y,z)
                        self.flow['p'][idx] = state.gas.p
                        self.flow['T'][idx] = state.gas.T
                        self.flow['rho'][idx] = state.gas.rho
                        self.flow['e'][idx] = state.gas.e
                        self.flow['YB'][idx] = state.gas.YB
                        self.flow['a'][idx] = state.gas.a;
                        self.flow['velx'][idx] = state.vel.x
                        self.flow['vely'][idx] = state.vel.y
                        self.flow['velz'][idx] = state.vel.z
                        idx += 1
        else:
            raise RuntimeError('Need to supply a FlowState object or a function that produces'+
                               'a FlowState object as a function of (x,y,z) location.')
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
        # Check in only those FlowStates associated with boundary conditions.
        for bc in bcs.values():
            if type(bc) is InflowBC:
                # [TODO] PJ 2023-01-21: check if a flowState is already present in the list
                # and add it only if not present.
                bc.fsi = len(flowStatesList)
                flowStatesList.append(bc.fs)
        #
        # For boundary conditions that are fixed-temperature and specified via
        # a user-supplied function, generate the specific values for each boundary face.
        for face in FaceList:
            bc = self.bcs[face]
            if type(bc) is WallNoSlipFixedTBC and bc.form == 'fun':
                bc.generate_values(self.grid, face)
        #
        self.i = i
        self.j = j
        self.k = k
        self.active = active;
        self.label = label
        self.indx = len(fluidBlocksList)
        fluidBlocksList.append(self)

    def to_json(self):
        """
        A representation of the block's configuration.
        """
        result = '{"i": %d, "j": %d, "k": %d,' % (self.i, self.j, self.k)
        result += ' "bcs": {'
        for i,f in enumerate(FaceList):
            result += '"%s": %s' % (f, self.bcs[f].to_json())
            result += ',' if i < len(FaceList)-1 else ''
        result += '},'
        result += ' "active": %s,' % json.dumps(self.active)
        result += ' "label": "%s"}' % self.label
        return result

    def construct_cell_properties(self):
        """
        We want the cell volumes and centroids available for writing
        into the initial flow file.
        """
        p000 = Vector3(x=self.grid.vertices.x[:-1, :-1, :-1],
                       y=self.grid.vertices.y[:-1, :-1, :-1],
                       z=self.grid.vertices.z[:-1, :-1, :-1])
        p100 = Vector3(x=self.grid.vertices.x[1:, :-1, :-1],
                       y=self.grid.vertices.y[1:, :-1, :-1],
                       z=self.grid.vertices.z[1:, :-1, :-1])
        p010 = Vector3(x=self.grid.vertices.x[:-1, 1:, :-1],
                       y=self.grid.vertices.y[:-1, 1:, :-1],
                       z=self.grid.vertices.z[:-1, 1:, :-1])
        p110 = Vector3(x=self.grid.vertices.x[1:, 1:, :-1],
                       y=self.grid.vertices.y[1:, 1:, :-1],
                       z=self.grid.vertices.z[1:, 1:, :-1])
        p001 = Vector3(x=self.grid.vertices.x[:-1, :-1, 1:],
                       y=self.grid.vertices.y[:-1, :-1, 1:],
                       z=self.grid.vertices.z[:-1, :-1, 1:])
        p101 = Vector3(x=self.grid.vertices.x[1:, :-1, 1:],
                       y=self.grid.vertices.y[1:, :-1, 1:],
                       z=self.grid.vertices.z[1:, :-1, 1:])
        p011 = Vector3(x=self.grid.vertices.x[:-1, 1:, 1:],
                       y=self.grid.vertices.y[:-1, 1:, 1:],
                       z=self.grid.vertices.z[:-1, 1:, 1:])
        p111 = Vector3(x=self.grid.vertices.x[1:, 1:, 1:],
                       y=self.grid.vertices.y[1:, 1:, 1:],
                       z=self.grid.vertices.z[1:, 1:, 1:])
        self.cellc, self.cellv = hexahedron_properties(p000, p100, p110, p010,
                                                       p001, p101, p111, p011)


    def write_flow_data_to_file(self, fileName, varNamesList, binaryData):
        """
        Write the data into the one file, one IO-variable after the other.

        Within each block, the order of the cells corresponds to
        the order expected for a VTK structured-grid.
        """
        f = open(fileName, 'wb') if binaryData else gzip.open(fileName, 'wt')
        for varName in varNamesList:
            if varName == 'posx': value = self.cellc.x
            elif varName == 'posy': value = self.cellc.y
            elif varName == 'posz': value = self.cellc.z
            elif varName == 'vol': value = self.cellv
            elif varName == 'p': value = self.flow['p']
            elif varName == 'T': value = self.flow['T']
            elif varName == 'rho': value = self.flow['rho']
            elif varName == 'e': value = self.flow['e']
            elif varName == 'YB': value = self.flow['YB']
            elif varName == 'a': value = self.flow['a']
            elif varName == 'velx': value = self.flow['velx']
            elif varName == 'vely': value = self.flow['vely']
            elif varName == 'velz': value = self.flow['velz']
            else: raise RuntimeError('unhandled variable name: ' + varName)
            if binaryData:
                buf = value.flatten().tobytes()
            else:
                text = np.array2string(value.flatten(), separator='\n', threshold=100_000_000)
                buf = text.strip('[]')+'\n'
            f.write(buf)
            # end varName loop
        f.close()
        return

# --------------------------------------------------------------------

def makeFBArray(i0=0, j0=0, k0=0, ni=1, nj=1, nk=1,
                grid=None, initialState=None, bcs={}, active=True, label=""):
    """
    This class is used to help build a description of the flow domain,
    several FluidBlocks at a time.
    A single grid is split into an array of subgrids and
    a FluidBlock is constructed for each subgrid.
    """
    global config, fluidBlocksList, flowStatesList
    #
    if not isinstance(grid, StructuredGrid):
        raise RuntimeError('Need to supply a StructuredGrid object to subdivide.')
    #
    if not (isinstance(initialState, FlowState) or callable(initialState)):
        raise RuntimeError('Need to supply a FlowState object or a function for initialState.')
    #
    # Boundary conditions
    # Set default values and then overwrite, if the user has supplied them.
    _bcs = {}
    for f in FaceList: _bcs[f] = WallWithSlipBC()
    for key in bcs.keys():
        if key in [Face.iminus, 'iminus']: _bcs['iminus'] = bcs[key]
        if key in [Face.iplus, 'iplus']: _bcs['iplus'] = bcs[key]
        if key in [Face.jminus, 'jminus']: _bcs['jminus'] = bcs[key]
        if key in [Face.iplus, 'jplus']: _bcs['jplus'] = bcs[key]
        if key in [Face.kminus, 'kminus']: _bcs['kminus'] = bcs[key]
        if key in [Face.kplus, 'kplus']: _bcs['kplus'] = bcs[key]
    #
    # Set up an array of FluidBlocks, dividing the overall grid into equal parts
    # or as close as reasonable.
    blks = []
    nic = (grid.niv-1) // ni
    njc = (grid.njv-1) // nj
    nkc = (grid.nkv-1) // nk
    newBCs = {}
    for k in range(nk):
        k0v = k * nkc
        nkv = nkc+1 if k < nk-1 else grid.nkv-k0v
        newBCs['kminus'] = _bcs['kminus'] if k==0 else ExchangeBC()
        newBCs['kplus'] = _bcs['kplus'] if k==nk-1 else ExchangeBC()
        for j in range(nj):
            j0v = j * njc
            njv = njc+1 if j < nj-1 else grid.njv-j0v
            newBCs['jminus'] = _bcs['jminus'] if j==0 else ExchangeBC()
            newBCs['jplus'] = _bcs['jplus'] if j==nj-1 else ExchangeBC()
            for i in range(ni):
                i0v = i * nic
                niv = nic+1 if k < ni-1 else grid.niv-i0v
                newBCs['iminus'] = _bcs['iminus'] if i==0 else ExchangeBC()
                newBCs['iplus'] = _bcs['iplus'] if i==ni-1 else ExchangeBC()
                newGrid = grid.subgrid(i0=i0v, j0=j0v, k0=k0v, niv=niv, njv=njv, nkv=nkv)
                newBlock = FluidBlock(i=i0+i, j=j0+j, k=k0+k, grid=newGrid,
                                      initialState=initialState, bcs=newBCs,
                                      active=active,
                                      label=label+("-%d-%d-%d".format(i, j, k)))
                blks.append(newBlock)
    return blks


# --------------------------------------------------------------------

def write_initial_files(binaryData):
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
        fp.write('  "nib": %d,\n' % config.nib)
        fp.write('  "njb": %d,\n' % config.njb)
        fp.write('  "nkb": %d,\n' % config.nkb)
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
        fileName = gridDir + ('/grid-%04d-%04d-%04d' % (fb.i, fb.j, fb.k))
        if binaryData:
            fb.grid.write_to_binary_file(fileName+'.bin')
        else:
            fb.grid.write_to_gzip_file(fileName+'.gz')
    #
    print('Write the initial flow-field files.')
    #
    flowDir = config.job_name+'/flow/t0000'
    if not os.path.exists(flowDir):
        os.makedirs(flowDir)
    for fb in fluidBlocksList:
        fileName = flowDir + ('/flow-%04d-%04d-%04d' % (fb.i, fb.j, fb.k))
        fileName += '.bin' if binaryData else '.gz'
        fb.write_flow_data_to_file(fileName, config.iovar_names, binaryData)
        #
        # For boundary conditions that are fixed-temperature and specified via
        # a user-supplied function, write out the specific values for each boundary face.
        for face in FaceList:
            bc = fb.bcs[face]
            if type(bc) is WallNoSlipFixedTBC and bc.form == 'fun':
                fileName = 'TWall-%04d-%04d-%04d-%s.gz' % (fb.i, fb.j, fb.k, face)
                bc.write_values_to_file(config.job_name+'/'+fileName)
    #
    with open(config.job_name+'/times.data', 'w') as fp:
        fp.write("# tindx t\n")
        fp.write("0 0.0\n")
    #
    print("End write initial files.")
    return

# --------------------------------------------------------------------

if __name__ == '__main__':
    print("Begin chicken preprocessing...")

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
        print("  flow states   :", len(flowStatesList))
        if len(fluidBlocksList) < 1:
            print("Warning: no fluid blocks defined; this is unusual.")
        else:
            nActive = 0;
            for b in fluidBlocksList:
                if b.active: nActive += 1
            print("  fluid blocks  :", nActive, "active", len(fluidBlocksList)-nActive, "inactive")
        config.check_array_of_fluid_blocks()
        binaryData = ("--binary" in uoDict) or ("-b" in uoDict)
        write_initial_files(binaryData)
    print("Done in {:.3f} seconds.".format(time.process_time()))
    sys.exit(0)
