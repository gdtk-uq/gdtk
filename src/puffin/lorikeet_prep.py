#! /usr/bin/env python3
"""
Python program to set up a simulation for the Lorikeet 2D Flow Solver.

It is intended for the user to define their particular domain and
flow in terms of the data objects defined in this module.
As part of its initialization, this program will read and execute
a user-specified input script that contains, in Python,
the code that defines the domain geometry and gas-path details.

Usage:
  $ lrkt-prep --job=<jobName>

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
  2022-12-11  First Python code adpated from chkn_prep.py and puffin-prep.py.
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

from gdtk.geom.vector3 import Vector3, quad_properties
from gdtk.geom.surface import CoonsPatch
from gdtk.geom.sgrid import StructuredGrid
from gdtk.numeric import roberts
from gdtk.gas import GasModel, GasState, ThermochemicalReactor

sys.path.append("") # so that we can find user's scripts in current directory

shortOptions = "hf:"
longOptions = ["help", "job="]

def printUsage():
    print("Prepare a lorikeet run by reading a job script and writing grid and flow files.")
    print("Usage: lorikeet-prep" + \
          " [--help | -h]" + \
          " [--job=<jobName> | -f <jobName>]")
    print("")
    return

#----------------------------------------------------------------------
# This is where we store the core data for the simulation.

FCNames = 'ausmdv hanel riemann ausmdv_plus_hanel riemann_plus_hanel'
FluxCalc = collections.namedtuple('_', FCNames)(*range(5))

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
    __slots__ = 'job_name', 'title', \
        'gas_model_file', 'gmodel', 'iovar_names', \
        'reaction_file_1', 'reaction_file_2', 'reactor', 'reacting', 'T_frozen', \
        'axisymmetric', \
        'nib', 'njb', 'blk_ids', 'nics', 'njcs', \
        'dt_init', 'cfl', 'cfl_count', 'print_count', \
        'plot_dt', 'max_time', 'max_step', \
        'x_order', 't_order', 'flux_calc', 'compression_tol', 'shear_tol'

    def __init__(self):
        """Accepts user-specified data and sets defaults. Make one only."""
        if GlobalConfig.count >= 1:
            raise Exception("Already have a GlobalConfig object defined.")
        #
        self.job_name = ""
        self.title = "Two-Dimensional Flow Simulation."
        self.gas_model_file = ""
        self.gmodel = None
        self.iovar_names = []
        self.reaction_file_1 = ""
        self.reaction_file_2 = ""
        self.reactor = None
        self.reacting = False
        self.T_frozen = 300.0
        self.axisymmetric = False
        self.nib = None
        self.njb = None
        self.dt_init = 1.0e-6
        self.cfl = 0.5
        self.cfl_count = 10;
        self.print_count = 20;
        self.plot_dt = 1.0e-3
        self.max_time = 1.0e-3
        self.max_step = 100
        self.x_order = 2
        self.t_order = 2
        self.flux_calc = FluxCalc.ausmdv
        self.compression_tol = -0.01
        self.shear_tol = 0.2
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
        fp.write('  "gas_model_file": "%s",\n' % self.gas_model_file)
        fp.write('  "iovar_names": %s,\n' % json.dumps(self.iovar_names))
        fp.write('  "reaction_file_1": "%s",\n' % self.reaction_file_1)
        fp.write('  "reaction_file_2": "%s",\n' % self.reaction_file_2)
        fp.write('  "reacting": %s,\n' % json.dumps(self.reacting))
        fp.write('  "T_frozen": %e,\n' % self.T_frozen)
        fp.write('  "axisymmetric": %s,\n' % json.dumps(self.axisymmetric))
        fp.write('  "max_time": %e,\n' % self.max_time)
        fp.write('  "max_step": %d,\n' % self.max_step)
        fp.write('  "dt_init": %e,\n' % self.dt_init)
        fp.write('  "cfl": %e,\n' % self.cfl)
        fp.write('  "cfl_count": %d,\n' % self.cfl_count)
        fp.write('  "print_count": %d,\n' % self.print_count)
        fp.write('  "x_order": %d,\n' % self.x_order)
        fp.write('  "t_order": %d,\n' % self.t_order)
        fp.write('  "flux_calc": "%s",\n' % self.flux_calc)
        fp.write('  "plot_dt": %e,\n' % self.plot_dt)
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
        self.nib = imax+1
        self.njb = jmax+1
        # As each block is defined in the input script, it will get a location in the array.
        # We will use the following array of ids to indicate which blocks in the array
        # have been defined.  We will allow some blocks to be left undefined.
        self.blk_ids = np.zeros((self.nib, self.njb), dtype=int)
        self.blk_ids -= 1 # Start with a negative id so that we can see undefined blocks.
        #
        # Identify the array positions for the fluid blocks that are defined and
        # check that their numbers of cells are consistent across the array.
        self.nics = np.zeros(self.nib, dtype=int)
        self.njcs = np.zeros(self.njb, dtype=int)
        for b in fluidBlocksList:
            i = b.i; j = b.j
            if self.blk_ids[i,j] >= 0:
                raise RuntimeError('Already a block at i=%d, j=%d' % (i, j))
            self.blk_ids[i,j] = b.indx
            nic = b.grid.niv-1
            njc = b.grid.njv-1
            if self.nics[i] == 0: self.nics[i] = nic
            if self.njcs[j] == 0: self.njcs[j] = njc
            if nic != self.nics[i]:
                raise RuntimeError('FluidBlock at i=%d j=%d k=%d has inconsistent nic=%d' % (i, j, nic))
            if njc != self.njcs[j]:
                raise RuntimeError('FluidBlock at i=%d j=%d k=%d has inconsistent njc=%d' % (i, j, njc))
        #
        return

# We will create just one GlobalConfig object that the user can alter.
config = GlobalConfig()

def init_gas_model(fileName, reaction_file_1="", reaction_file_2=""):
    """
    Initialize a GasModel object for use in the user's input script.

    fileName: (string) Name of the detailed-gas-model file.
    reaction_file_1: (string) Name of the detailed chemistry file for reacting gas.
    reaction_file_2: (string) Name of the second thermochemistry file.
    This second thermochemistry file is needed for only a few of the multi-T models.
    """
    global config
    if config.gmodel:
        raise RuntimeError("Already have a gas model initialized.")
    if not os.path.exists(fileName):
        raise Exception("Gas model file not found: " + fileName)
    config.gmodel = GasModel(fileName)
    config.gas_model_file = fileName
    # Now that we have a gas model, we can construct the list of names for
    # the IO variables.  We use underscores for the species names so that
    # we can easily pick out their index values later.
    config.iovar_names = ['posx', 'posy', 'vol', 'p', 'T', 'rho', 'e', 'a']
    for isp in range(config.gmodel.n_species):
        name = 'massf_%d_%s' % (isp, config.gmodel.species_names[isp])
        config.iovar_names.append(name)
    config.iovar_names += ['velx', 'vely']
    #
    # ThermochemicalReactor is always associated with a particular GasModel.
    if reaction_file_1:
        if not os.path.exists(reaction_file_1):
            raise Exception("Reaction file 1 not found: " + reaction_file_1)
        if reaction_file_2 and (not os.path.exists(reaction_file_2)):
            raise Exception("Reaction file 2 not found: " + reaction_file_2)
        reactor = ThermochemicalReactor(gmodel, reaction_file_1, reaction_file_2)
    else:
        reactor = None
    config.reactor = reactor
    config.reaction_file_1 = reaction_file_1
    config.reaction_file_2 = reaction_file_2
    return


#----------------------------------------------------------------------
# The following classes define the objects that will be assembled into
# the simulation.
# ---------------------------------------------------------------------
#
# We will accumulate references to defined objects.
fluidBlocksList = []

class FlowState():
    """
    Contains the gas properties and velocities defining a flow state.

    The user may create more than one FlowState to describe the initial
    gas properties and boundary conditions throughout the flow domain.
    """

    __slots__ = 'indx', 'label', 'gas', 'vel'

    def __init__(self, gas=None, p=100.0e3, T=300.0, massf=[1.0],
                 vel=None, velx=0.0, vely=0.0, velz=0.0,
                 label=""):
        global config
        self.gas = GasState(config.gmodel)
        if isinstance(gas, GasState):
            self.gas.copy_values(gas)
        else:
            self.gas.p=p; self.gas.T=T; self.gas.massf=massf
            self.gas.update_thermo_from_pT()
        if isinstance(vel, Vector3):
            self.vel = copy(vel)
        else:
            self.vel = Vector3(velx, vely)
        self.label = label

    def __repr__(self):
        return "FlowState(gas={}, vel={})".format(self.gas, self.vel)

    def to_json(self):
        gas_dict = {'p': self.gas.p, 'T': self.gas.T, 'massf': self.gas.massf}
        result = '{"gas": %s, "vel": [%g, %g]}' % (json.dumps(gas_dict), self.vel.x, self.vel.y)
        return result

# There is a bit of boiler place code here but it gives us
# a nice notation for use in the user's input script.

FaceNames = 'iminus iplus jminus jplus'
Face = collections.namedtuple('_', FaceNames)(*range(4))
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


class InflowBC(BoundaryCondition):
    """
    Inflow boundary condition.
    """
    __slots__ = ['tag', 'fs']

    def __init__(self, fs):
        self.tag = 'inflow'
        if type(fs) is FlowState:
            self.fs = fs
        else:
            raise RuntimeError("Inflow boundary condition expects a FlowState object.")

    def __repr__(self):
        return "InflowBC(fs={})".format(self.fs)

    def to_json(self):
        return '{"tag": "%s", "flow_state": %s}' % (self.tag, self.fs.to_json())


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
    __slots__ = 'indx', 'active', 'i', 'j', 'grid', 'initialState', \
        'nic', 'njc', 'bcs', 'label', \
        'cellc', 'cellv', 'flow'

    def __init__(self, i=0, j=0, grid=None, initialState=None, bcs={}, active=True, label=""):
        # When constructing a FluidBlock, we fill out the full grid and
        # the flow data in each cell within the block.
        # This allow us to use numpy arrays, for speed.
        global config, fluidBlocksList
        if isinstance(grid, StructuredGrid):
            self.grid = grid
        else:
            raise RuntimeError('Need to supply a StructuredGrid object.')
        self.construct_cell_properties()
        #
        self.nic = grid.niv - 1
        self.njc = grid.njv - 1
        assert(self.nic==self.cellc.x.shape[0])
        assert(self.njc==self.cellc.y.shape[1])
        ncells = self.nic*self.njc
        #
        # The data constructed for the flow dictionary must cover all of the iovar_names.
        #
        if isinstance(initialState, FlowState):
            self.flow = {'posx':self.cellc.x, 'posy':self.cellc.y, 'vol':self.cellv,
                         'p': np.full(ncells, initialState.gas.p), 'T':np.full(ncells, initialState.gas.T),
                         'rho':np.full(ncells, initialState.gas.rho), 'e':np.full(ncells, initialState.gas.u),
                         'a':np.full(ncells, initialState.gas.a),
                         'velx': np.full(ncells, initialState.vel.x), 'vely': np.full(ncells, initialState.vel.y)}
            for isp in range(config.gmodel.n_species):
                name = 'massf_%d_%s' % (isp, config.gmodel.species_names[isp])
                self.flow[name] = np.full(ncells, initialState.gas.massf[isp])
        elif callable(initialState):
            # The user has supplied a function,
            # which we use one cell-centre at a time.
            self.flow = {'posx':np.zeros(ncells, dtype=float), 'posy':np.zeros(ncells, dtype=float),
                         'vol':np.zeros(ncells, dtype=float),
                         'p':np.zeros(ncells, dtype=float), 'T':np.zeros(ncells, dtype=float),
                         'rho':np.zeros(ncells, dtype=float), 'e':np.zeros(ncells, dtype=float),
                         'a':np.zeros(ncells, dtype=float),
                         'velx':np.zeros(ncells, dtype=float), 'vely':np.zeros(ncells, dtype=float)}
            for isp in range(config.gmodel.n_species):
                name = 'massf_%d_%s' % (isp, config.gmodel.species_names[isp])
                self.flow[name] = np.zeros(ncells, dtype=float)
            idx = 0
            # Note VTK ordering of loops; the IO functions assume this ordering.
            for jj in range(self.njc):
                for ii in range(self.nic):
                    x = self.cellc.x[ii,jj]; y = self.cellc.y[ii,jj]
                    self.flow['posx'][idx] = x
                    self.flow['posy'][idx] = y
                    self.flow['vol'][idx] = self.cellv[ii,jj]
                    state = initialState(x,y)
                    self.flow['p'][idx] = state.gas.p
                    self.flow['T'][idx] = state.gas.T
                    self.flow['rho'][idx] = state.gas.rho
                    self.flow['e'][idx] = state.gas.u
                    for isp in range(config.gmodel.n_species):
                        name = 'massf_%d_%s' % (isp, config.gmodel.species_names[isp])
                        self.flow[name][idx] = state.gas.massf[isp]
                    self.flow['a'][idx] = state.gas.a;
                    self.flow['velx'][idx] = state.vel.x
                    self.flow['vely'][idx] = state.vel.y
                    idx += 1
        else:
            raise RuntimeError('Need to supply a FlowState object or a function that produces'+
                               'a FlowState object as a function of (x,y) location.')
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
        #
        self.i = i
        self.j = j
        self.active = active;
        self.label = label
        self.indx = len(fluidBlocksList)
        fluidBlocksList.append(self)

    def to_json(self):
        """
        A representation of the block's configuration.
        """
        result = '{"i": %d, "j": %d, ' % (self.i, self.j)
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
        global config
        p00 = Vector3(x=self.grid.vertices.x[:-1, :-1], y=self.grid.vertices.y[:-1, :-1])
        p10 = Vector3(x=self.grid.vertices.x[1:, :-1], y=self.grid.vertices.y[1:, :-1])
        p01 = Vector3(x=self.grid.vertices.x[:-1, 1:], y=self.grid.vertices.y[:-1, 1:])
        p11 = Vector3(x=self.grid.vertices.x[1:, 1:], y=self.grid.vertices.y[1:, 1:])
        # In 2D planar geometry, the cells are assumed to be defined in the (x,y)-plane
        # and, with unit depth in the z-direction, the volume per unit depth is the x,y area.
        self.cellc, nrm, t1, t2, self.cellv = quad_properties(p00, p10, p11, p01)
        # In 2D axisymmetric geometry, the volume is per radian swept about the x-axis.
        if config.axisymmetric: self.cellv *= self.cellc.y


    def write_flow_data_to_file(self, fileName, varNamesList):
        """
        Write the data into the one file, one IO-variable after the other.

        Within each block, the order of the cells corresponds to
        the order expected for a VTK structured-grid.
        """
        f = gzip.open(fileName, 'wt')
        for varName in varNamesList:
            text = np.array2string(self.flow[varName].flatten(), separator='\n', threshold=100_000_000)
            buf = text.strip('[]')+'\n'
            f.write(buf)
        f.close()
        return

# --------------------------------------------------------------------

def write_initial_files():
    """
    Writes the files needed for the main simulation code.

    These files are found in the directory config.job_name.
    """
    global config, fluidBlocksList
    print("Begin write initial files.")
    #
    if not os.path.exists(config.job_name):
        os.mkdir(config.job_name)
    #
    with open(config.job_name+'/config.json', 'w') as fp:
        fp.write("{\n")
        config.write(fp)
        #
        n_fluid_blocks = len(fluidBlocksList)
        fp.write('  "n_fluid_blocks": %d,\n' % n_fluid_blocks)
        fp.write('  "nib": %d,\n' % config.nib)
        fp.write('  "njb": %d,\n' % config.njb)
        fp.write('  "blk_ids": %s,\n' % json.dumps(config.blk_ids.tolist()))
        fp.write('  "nics": %s,\n' % json.dumps(config.nics.tolist()))
        fp.write('  "njcs": %s,\n' % json.dumps(config.njcs.tolist()))
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
        fileName = gridDir + ('/grid-%04d-%04d.gz' % (fb.i, fb.j))
        fb.grid.write_to_gzip_file(fileName)
    #
    print('Write the initial flow-field files.')
    #
    flowDir = config.job_name+'/flow/t0000'
    if not os.path.exists(flowDir):
        os.makedirs(flowDir)
    for fb in fluidBlocksList:
        fileName = flowDir + ('/flow-%04d-%04d.gz' % (fb.i, fb.j))
        fb.write_flow_data_to_file(fileName, config.iovar_names)
    #
    with open(config.job_name+'/times.data', 'w') as fp:
        fp.write("# tindx t\n")
        fp.write("0 0.0\n")
    #
    print("End write initial files.")
    return

# --------------------------------------------------------------------

if __name__ == '__main__':
    print("Begin lorikeet-prep...")

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
        #
        # The user-specified input comes in the form of Python code.
        # It is up to the user to be careful; there is no security.
        #
        exec(compile(open(inputScriptName, "rb").read(), inputScriptName, 'exec'))
        #
        print("Summary:")
        print("  iovar_names   :", config.iovar_names)
        print("  fluid blocks  :", len(fluidBlocksList))
        if len(fluidBlocksList) < 1:
            print("Warning: no fluid blocks defined; this is unusual.")
        config.check_array_of_fluid_blocks()
        write_initial_files()
    print("Done in {:.3f} seconds.".format(time.process_time()))
    sys.exit(0)
