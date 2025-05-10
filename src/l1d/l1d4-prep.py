#! /usr/bin/env python3
"""
Python program to set up a simulation for the Lagrangian 1D Flow Solver.

It is intended for the user to define their particular facility and
flow in terms of the data objects defined in this module.
As part of its initialization, this program will read and execute
a user-specified input script that contains, in Python,
the code that defines the facility geometry and gas-path details.

Usage:
  $ l1d4-prep --job=<jobName>

The simulation control data is then organised via the classes
GlobalConfig, GasSlug, Piston, EndCondition and its subclasses.
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
  2005-Jun     First Python code
  2006-Jul-24  Ported to Rowan's new C++ chemistry
  2010-Mar-May More reworking for L1d3 and the latest-greatest thermochemistry
  2012-Sep-Oct Much cleaning up and Sphinx docs.
  2020-Apr-04  L1d4 flavour started
  2021-Apr-16  Matt McGilvray's heat-transfer coefficient factor.
  2021-Apr-23  Improve piston+buffer interaction.
  2025-Mar-19  Allow variable piston mass to model the secondary diaphragm
               in an expansion tube, where the diaphragm is initially punched
               out of its restraint and temporarily acts as a piston.
  2025-May-10  Implement the code for writing of history cells.
"""

# ----------------------------------------------------------------------
#
import sys
import os
from getopt import getopt
import math
import numpy as np
import json
from gdtk.gas import GasModel, GasState, ThermochemicalReactor
from gdtk.numeric import roberts

sys.path.append("") # so that we can find user's scripts in current directory

shortOptions = "hf:"
longOptions = ["help", "job="]

def printUsage():
    print("")
    print("Usage: l1d-prep" + \
          " [--help | -h]" + \
          " [--job=<jobName> | -f <jobName>]")
    print("")
    return

#----------------------------------------------------------------------
# This is where we store the core data for the simulation.

class GlobalConfig(object):
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
                'gas_model_files', 'gmodels', \
                'overall_species_count', 'overall_species_index', \
                'overall_modes_count', 'overall_modes_index', \
                'reaction_files_1', 'reaction_files_2', \
                'reactors', 'reacting', 'T_frozen', \
                'dt_init', 'cfl_list', 'cfl_count', 'print_count', 'dt_plot_list', \
                'max_time', 'max_step', 'x_order', 't_order', \
                'hloc_list'

    def __init__(self):
        """Accepts user-specified data and sets defaults. Make one only."""
        if GlobalConfig.count >= 1:
            raise Exception("Already have a GlobalConfig object defined.")

        self.job_name = ""
        self.title = "Another L1d4 Simulation."
        self.gas_model_files = []
        self.gmodels = []
        # History points will have to deal with species from all gas models.
        self.overall_species_count = 0
        self.overall_species_index = []
        self.overall_modes_count = 0
        self.overall_modes_index = []
        self.reaction_files_1 = []
        self.reaction_files_2 = []
        self.reactors = []
        self.reacting = False
        self.T_frozen = 300.0
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
        # History locations are collected in a list.
        self.hloc_list = []
        #
        GlobalConfig.count += 1
        return

    def write(self, fp):
        """
        Writes the configuration data to the specified file in JSON format.
        """
        global config, slugList, pistonList, diaphragmList
        #
        fp.write('"config": {\n')
        fp.write('  "title": "%s",\n' % self.title)
        fp.write('  "gas_model_files": %s,\n' % json.dumps(self.gas_model_files))
        fp.write('  "reaction_files_1": %s,\n' % json.dumps(self.reaction_files_1))
        fp.write('  "reaction_files_2": %s,\n' % json.dumps(self.reaction_files_2))
        fp.write('  "reacting": %s,\n' % json.dumps(self.reacting))
        fp.write('  "T_frozen": %e,\n' % self.T_frozen)
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
        #
        n_hloc = len(config.hloc_list)
        fp.write('  "hloc_n": %d,\n' % n_hloc)
        xlist = [self.hloc_list[i] for i in range(n_hloc)]
        fp.write('  "hloc_x": %s,\n' % json.dumps(xlist))
        #
        fp.write('  "nslugs": %d,\n' % len(slugList))
        fp.write('  "npistons": %d,\n' % len(pistonList))
        fp.write('  "nvalves": %d,\n' % len(valveList))
        fp.write('  "ndiaphragms": %d,\n' % len(diaphragmList))
        fp.write('  "necs": %d\n' % len(ecList))
        # Note, no comma after last item inside JSON dict
        fp.write('},\n') # comma here because not the last item in file
        return

# We will create just one GlobalConfig object that the user can alter.
config = GlobalConfig()

# --------------------------------------------------------------------
# The following functions are to provide convenient ways of setting
# some of the GlobalConfig elements.

def add_gas_model(fileName, reaction_file_1="", reaction_file_2=""):
    """
    Initialize and add a GasModel in the GlobalConfig list of gas models.

    fileName: (string) Name of the detailed-gas-model file.
    reaction_file_1: (string) Name of the detailed chemistry file for reacting gas.
    reaction_file_2: (string) Name of the second thermochemistry file.
    This second thermochemistry file is needed for only a few of the multi-T models.
    """
    global config
    if not os.path.exists(fileName):
        raise Exception("Gas model file not found: " + fileName)
    gmodel = GasModel(fileName)
    gmodel_id = len(config.gmodels)
    config.gmodels.append(gmodel)
    config.gas_model_files.append(fileName)
    nsp = gmodel.n_species
    config.overall_species_index.append([config.overall_species_count+i for i in range(nsp)])
    config.overall_species_count += nsp
    nmodes = gmodel.n_modes
    config.overall_modes_index.append([config.overall_modes_count+i for i in range(nmodes)])
    config.overall_modes_count += nmodes
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
    config.reactors.append(reactor)
    config.reaction_files_1.append(reaction_file_1)
    config.reaction_files_2.append(reaction_file_2)
    #
    return gmodel_id


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


def add_history_loc(x):
    """
    Add a location to the history-location list in GlobalConfig.

    x: (float) x-coordinate, in metres, of the sample point.
    Returns: Number of sample points defined so far.
    """
    global config
    if isinstance(x, list):
        config.hloc_list.extend(x)
    else:
        config.hloc_list.append(float(x))
    return len(config.hloc_list)

# --------------------------------------------------------------------

class Tube(object):
    """
    Contains the tube specification.

    The user's script should not create one of these
    but should specify the tube details by calling the add_xxxx functions.
    """
    count = 0

    __slots__ = 'n', 'x_list', 'd_list', 'T_nominal', 'T_patch_list', \
        'loss_region_list', 'viscous_factor_patch_list', 'htc_factor_patch_list', \
        'xs', 'ds', 'K_over_Ls', 'Ts', 'vfs', 'htcfs'

    def __init__(self):
        """Accepts user-specified data and sets defaults. Make one only."""
        if Tube.count >= 1:
            raise Exception("Already have a Tube object defined.")
        # Tube definition is a list of (x,diameter) tuples
        # defining the break-points of the tube wall.
        # The transition_flag with a value of 1
        # indicates linear transitions from break-point i to point i+1.
        # The alternative is a cubic transition (I think).
        self.x_list = []
        self.d_list = []
        self.n = 4000
        # Wall temperature is specified as a nominal value with
        # patches of other temperatures.
        self.T_nominal = 300.0  # Nominal tube-wall temperature in Kelvin
        self.T_patch_list = []
        # Head-losses are also spread over finite-length patches.
        self.loss_region_list = []
        # Sections of tube may have reduced viscous effects.
        self.viscous_factor_patch_list = []
        # Sections may also have augmentation of the heat transfer.
        self.htc_factor_patch_list = []
        #
        Tube.count += 1
        return

    def set_up_tables(self):
        """
        Set up the lists that define the "discretized" tube.
        """
        nsegments = len(self.x_list) - 1
        assert len(self.d_list) == nsegments+1, "Mismatch in tube.x_list, tube.d_list lengths."
        if nsegments < 1:
            raise Exception("You did not specify at least two points in your tube.")
        self.xs = np.linspace(tube.x_list[0], tube.x_list[-1], num=tube.n+1)
        self.ds = []
        for x in self.xs:
            iseg = nsegments
            for i in range(nsegments):
                if x <= self.x_list[i]:
                    iseg = i
                    break
            # Have selected the segment containing x.
            frac = (x - self.x_list[iseg-1]) / (self.x_list[iseg] - self.x_list[iseg-1])
            d = (1.0-frac)*self.d_list[iseg-1] + frac*self.d_list[iseg]
            self.ds.append(d)
        #
        self.K_over_Ls = []
        for x in self.xs:
            value = 0.0
            for region in self.loss_region_list:
                xL = region['xL']; xR = region['xR']
                if x >= xL and x <= xR: value = region['K']/(xR-xL)
            self.K_over_Ls.append(value)
        #
        self.Ts = []
        for x in self.xs:
            value = self.T_nominal
            for region in self.T_patch_list:
                xL = region['xL']; xR = region['xR']
                if x >= xL and x <= xR: value = region['T']
            self.Ts.append(value)
        #
        self.vfs = []
        for x in self.xs:
            value = 1.0
            for region in self.viscous_factor_patch_list:
                xL = region['xL']; xR = region['xR']
                if x >= xL and x <= xR: value = region['vf']
            self.vfs.append(value)
        #
        self.htcfs = []
        for x in self.xs:
            value = 1.0
            for region in self.htc_factor_patch_list:
                xL = region['xL']; xR = region['xR']
                if x >= xL and x <= xR: value = region['htcf']
            self.htcfs.append(value)
        #
        return

    def eval(self, x):
        """
        Computes tube cross-section properties at position x.
        """
        if x <= self.xs[0]:
            d = self.ds[0]
            K_over_L = self.K_over_Ls[0]
            Twall = self.Ts[0]
            vf = self.vfs[0]
            htcf = self.htcfs[0]
        elif x >= self.xs[-1]:
            d = self.ds[-1]
            K_over_L = self.K_over_Ls[-1]
            Twall = self.Ts[-1]
            vf = self.vfs[-1]
            htcf = self.htcfs[-1]
        else:
            dx = self.xs[1] - self.xs[0]
            i = int((x - self.xs[0])/dx)
            frac = (x - self.xs[i])/dx
            d = (1.0-frac)*self.ds[i] + frac*self.ds[i+1]
            K_over_L = (1.0-frac)*self.K_over_Ls[i] + frac*self.K_over_Ls[i+1]
            Twall = (1.0-frac)*self.Ts[i] + frac*self.Ts[i+1]
            vf = (1.0-frac)*self.vfs[i] + frac*self.vfs[i+1]
            htcf = (1.0-frac)*self.htcfs[i] + frac*self.htcfs[i+1]
        # We compute the area from diameter, assuming a circular cross-section.
        # It is put into the file for later plotting.
        area = math.pi*(d**2)/4
        return (d, area, K_over_L, Twall, vf, htcf)

    def write(self, fp):
        """
        Writes the tube specification to the specified file, in small steps.
        """
        fp.write('# n= %d\n' % self.n) # n+1 points along tube to follow
        fp.write('# 1:x,m  2:d,m  3:area,m^2  4:K_over_L,1/m  5:Twall,K ' +
                 ' 6:viscous-factor  7:htc-factor\n')
        for i in range(len(self.xs)):
            fp.write('%e %e %e %e %e %e %e\n' %
                     (self.xs[i], self.ds[i], math.pi*(self.ds[i]**2)/4,
                      self.K_over_Ls[i], self.Ts[i], self.vfs[i], self.htcfs[i]))
        return

# We will create just one Tube object that the user can alter.
tube = Tube()

# --------------------------------------------------------------------
# The following functions are to provide convenient ways of setting
# some of the Tube elements.

def add_break_point(x, d):
    """
    Add a break-point tuple to the tube-diameter description.

    The tube is described as a set of (x,d)-coordinate pairs that
    define break points in the profile of the tube wall.
    You need at least 2 break points to define the tube.
    Linear variation of diameter between the break points is assumed.

    x: (float) x-coordinate, in metres, of the break point
    d: (float) diameter, in metres, of the tube wall at the break-point.
    Returns: Number of break points defined so far.
    """
    global tube
    if len(tube.x_list) > 0:
        # Check that we are adding points monotonically in x.
        if x > tube.x_list[-1]:
            tube.x_list.append(x); tube.d_list.append(d)
        else:
            print("Warning: did not add new break-point (", x, d, ").")
    else:
        tube.x_list.append(x); tube.d_list.append(d)
    return len(tube.x_list)


def add_loss_region(xL, xR, K):
    """
    Add a head-loss region to the tube description.

    There is a momentum-sink term much like the so-called minor-loss terms
    in the fluid mechanics text books.
    The effect of the loss is spread over a finite region so that the cells
    are gradually affected as they pass through the region

    xL: (float) Left-end location, in metres, of the loss region.
    xR: (float) Right-end location, in metres, of the loss region.
    K: (float) Head-loss coefficient.  A value of 0.25 seems to be good for a
        reasonably smooth contraction such as the T4 main diaphragm station.
    Returns: Number of loss regions defined so far.
    """
    global tube
    if xR < xL:
        # Keep x-values in increasing order
        xL, xR = xR, xL
    if abs(xR - xL) < 1.0e-3:
        print("Warning: loss region is very short: (", xL, xR, ")")
    tube.loss_region_list.append({'xL':xL, 'xR':xR, 'K':K})
    return len(tube.loss_region_list)


def add_T_patch(xL, xR, T):
    """
    Add a temperature patch for a region where the wall temperature
    is different from the nominal value.

    xL: (float) Left-end location, in metres, of the temperature patch.
    xR: (float) Right-end location, in metres, of the temperature patch.
    T: (float) Wall temperature in degrees K.
    Returns: Number of temperature patches defined so far.
    """
    if xR < xL:
        # Keep x-values in increasing order
        xL, xR = xR, xL
    if abs(xR - xL) < 1.0e-3:
        print("Warning: temperature patch is very short: (", xL, xR, ")")
    tube.T_patch_list.append({'xL':xL, 'xR':xR, 'T':T})
    return len(tube.T_patch_list)


def add_vf_patch(xL, xR, vf):
    """
    Add a viscous-factor patch for a region where the viscous effects
    are scaled from the nominal value.

    xL: (float) Left-end location, in metres, of the loss region.
    xR: (float) Right-end location, in metres, of the loss region.
    vf: (float) Viscous factor for limiting viscous effects at wall.
        Nominal value is 1.0.  A completely inviscid wall has a value of 0.0.
    Returns: Number of viscous-factor patches defined so far.
    """
    if xR < xL:
        # Keep x-values in increasing order
        xL, xR = xR, xL
    tube.viscous_factor_patch_list.append({'xL':xL, 'xR':xR, 'vf':vf})
    return len(tube.viscous_factor_patch_list)


def add_htcf_patch(xL, xR, htcf):
    """
    Add a heat-transfer-coefficient factor patch for a region where
    heat-transfer is scaled from the nominal value.
    This may be used to augment the heat hransfer is twisty little
    passages as found around the valve in the Oxford High-Density Tunnel.

    xL: (float) Left-end location, in metres, of the augmentation region.
    xR: (float) Right-end location, in metres, of the augmentation region.
    htcf: (float) Heat-transfer factor for modulating the computed heat transfer.
        Nominal value is 1.0, to get the heat transfer for a simple circular pipe.
    Returns: Number of heat-transfer-coefficient factor patches defined so far.
    """
    if xR < xL:
        # Keep x-values in increasing order
        xL, xR = xR, xL
    tube.htc_factor_patch_list.append({'xL':xL, 'xR':xR, 'htcf':htcf})
    return len(tube.htc_factor_patch_list)


#----------------------------------------------------------------------
# The following classes define the objects that will be assembled into
# the gas path.
# ---------------------------------------------------------------------

# We will accumulate references to defined objects.
slugList = []
pistonList = []
valveList = []
ecList = []
diaphragmList = []
interfaceList = []
freeEndList = []
velocityEndList = []
pistonFaceList = []

class GasSlug():
    """
    Contains the gas properties and discretisation for each gas slug.

    The user may create more than one gas slug to describe the initial
    gas properties throughout the facility.

    Note that a slug needs to have appropriate end-conditions.
    This is achieved by creating end-condition objects such as
    FreeEnd and VelocityEnd objects and then assembling
    the gas-path via a call to the function assemble_gas_path.
    """

    __slots__ = 'indx', 'label', \
                'gas', 'gmodel', 'gmodel_id', \
                'initial_fs_fun', 'vel', 'xL', 'xR', \
                'ecL', 'ecR', \
                'ncells', 'to_end_L', 'to_end_R', 'cluster_strength', \
                'viscous_effects', 'adiabatic', 'hcells', \
                'ifxs'

    def __init__(self,
                 gmodel_id = None,
                 p = 100.0e3,
                 vel = 0.0,
                 T = 300.0,
                 T_modes = [],
                 massf = [1.0,],
                 initial_fs_fun = None,
                 label="",
                 ncells = 10,
                 to_end_L=False,
                 to_end_R=False,
                 cluster_strength=0.0,
                 viscous_effects=0, # several options were available in L1d3
                 adiabatic=False,
                 hcells=[],
                 ):
        """
        Creates a gas slug with user-specified properties.

        Most parameters have default properties so that only the user
        needs to override the ones that they wish to set differently.

        If an initial_fs_fn is supplied, it is expected to be a user-defined function
        that accepts the x-position (in metres), and it returns a dictionary with
        appropriate data (p, T, T_modes, massf, vel) for that x-position.
        See the write_cell_data() function, below, to see what data is expected.
        """
        # Gas data related values
        self.gmodel_id = gmodel_id
        self.gmodel = config.gmodels[gmodel_id]
        self.gas = GasState(self.gmodel)
        nsp = self.gmodel.n_species
        nmodes = self.gmodel.n_modes
        if initial_fs_fun is None:
            # We assemble a single flow state from the pieces of flowstate data
            # that are (presumably) present.
            self.gas.p = p
            self.gas.massf = massf
            self.gas.T = T
            if len(T_modes) != nmodes:
                raise Exception("T_modes is inconsistent. gmodel.n_modes=%d T_modes=%s"
                                % (nmodes, T_modes));
            self.gas.T_modes = T_modes
            self.gas.update_thermo_from_pT()
            self.gas.update_sound_speed()
            self.gas.update_trans_coeffs()
            self.vel = vel
            self.initial_fs_fun = None
        else:
            # If we have something for initial_fs_fn, we will try to use it.
            if not callable(initial_fs_fun):
                raise RuntimeError("Was expecting a callable function for initial_fs_fun.")
            self.initial_fs_fun = initial_fs_fun
        #
        self.label = label
        #
        self.ncells = ncells
        self.to_end_L = to_end_L
        self.to_end_R = to_end_R
        self.cluster_strength = cluster_strength
        #
        self.viscous_effects = viscous_effects
        self.adiabatic = adiabatic
        if isinstance(hcells,int):
            self.hcells=[hcells,]
        elif isinstance(hcells,list):
            self.hcells=hcells
        else:
            print("Warning: hcells reset to empty list.")
            hcells = []
        #
        # Boundary object at each end of the slug will be
        # attached later when the gas-path is assembled.
        self.ecL = None
        self.ecR = None
        # The spatial limits of the gas slug will be determined later,
        # from the boundary-condition objects.
        self.xL = None
        self.xR = None
        #
        # The GasSlug objects need an identity that can be
        # transferred to the main simulation program.
        global slugList
        self.indx = len(slugList) # next available index
        slugList.append(self)
        return

    def write_config(self, fp):
        """
        Writes the configuration data in JSON format.
        """
        fp.write('"slug_%d": {\n' % self.indx)
        fp.write('  "label": %s,\n' % json.dumps(self.label))
        fp.write('  "gmodel_id": %d,\n' % self.gmodel_id)
        fp.write('  "ncells": %d,\n' % self.ncells)
        fp.write('  "viscous_effects": %d,\n' % self.viscous_effects)
        fp.write('  "adiabatic": %s,\n' % json.dumps(self.adiabatic))
        fp.write('  "ecL_id": %d,\n' % self.ecL.ecindx)
        fp.write('  "ecR_id": %d,\n' % self.ecR.ecindx)
        hncell = len(self.hcells)
        fp.write('  "nhcells": %d,\n' % len(self.hcells))
        fp.write('  "hcells": %s\n' % json.dumps(self.hcells))
        # Note no comma after last item in JSON dict.
        fp.write('},\n') # presume that this dict not the last
        return

    def construct_cells_and_faces(self):
        self.ifxs = roberts.distribute_points_1(self.xL, self.xR, self.ncells,
                                                self.to_end_L, self.to_end_R,
                                                self.cluster_strength)
        return

    def write_face_data(self, fp, tindx, write_header):
        """
        Write the initial state of all of the interfaces for the gas slug.
        """
        if write_header:
            fp.write("#   x   area\n")
        fp.write("# tindx %d\n" % tindx)
        for x in self.ifxs:
            d, area, K_over_L, Twall, vf, htcf = tube.eval(x)
            fp.write("%e %e\n" % (x, area))
        fp.write("# end\n")
        return

    def write_cell_data(self, fp, tindx, write_header):
        """
        Write the initial state of the cells within the gas slug in GNUPlot format.
        """
        nsp = self.gmodel.n_species
        nmodes = self.gmodel.n_modes
        if write_header:
            fp.write('# xmid  volume  vel  L_bar  rho  p  T  u  a')
            fp.write('  shear_stress  heat_flux')
            for i in range(nsp): fp.write('  massf[%d]' % i)
            if nsp > 1: fp.write('  dt_chem')
            for i in range(nmodes): fp.write('  T_modes[%d]  u_modes[%d]' % (i, i))
            if nmodes > 0: fp.write('  dt_therm')
            fp.write('\n')
        fp.write("# tindx %d\n" % tindx)
        L_bar = 0.0; dt_chem = -1.0; dt_therm = -1.0
        shear_stress=0.0; heat_flux = 0.0
        for j in range(self.ncells):
            xmid = 0.5*(self.ifxs[j+1] + self.ifxs[j])
            d, area, K_over_L, Twall, vf, htcf = tube.eval(xmid)
            volume = area * (self.ifxs[j+1] - self.ifxs[j])
            #
            if callable(self.initial_fs_fun):
                fs = self.initial_fs_fun(xmid)
                # We assemble a local flow state from the pieces of flow state data.
                # that are (presumably) present in the returned dictionary.
                self.gas.p = fs['p']
                self.gas.massf = fs['massf']
                self.gas.T = fs['T']
                if nmodes > 0:
                    if len(fs['T_modes']) != nmodes:
                        raise Exception("T_modes is inconsistent. gmodel.n_modes=%d T_modes=%s"
                                        % (nmodes, fs['T_modes']));
                    self.gas.T_modes = fs['T_modes']
                else:
                    self.gas.T_modes = []
                self.gas.update_thermo_from_pT()
                self.gas.update_sound_speed()
                self.gas.update_trans_coeffs()
                self.vel = fs['vel']
            #
            fp.write('%e %e %e %e' % (xmid, volume, self.vel, L_bar))
            fp.write(' %e %e %e %e' % (self.gas.rho, self.gas.p, self.gas.T, self.gas.u))
            fp.write(' %e %e %e' % (self.gas.a, shear_stress, heat_flux))
            for i in range(nsp): fp.write(' %e' % (self.gas.massf[i]))
            if nsp > 1: fp.write(' %e' % dt_chem)
            for i in range(nmodes):
                fp.write(' %e %e' % (self.gas.T_modes[i], self.gas.u_modes[i]))
            if nmodes > 0: fp.write(' %e' % dt_therm)
            fp.write('\n')
        fp.write("# end\n")
        return

    def write_initial_flow_data_for_any_cell(self, fp, t, x):
        """
        Write the flow data for a history cell or a history location.

        We assume that the data is constant for the whole slug.
        """
        L_bar = 0.0; shear_stress=0.0; heat_flux = 0.0
        fp.write('%e %e %e %e' % (t, x, self.vel, L_bar))
        fp.write(' %e %e %e %e' % (self.gas.rho, self.gas.p, self.gas.T, self.gas.u))
        fp.write(' %e %e %e' % (self.gas.a, shear_stress, heat_flux))
        massf = [0.0,]*config.overall_species_count
        gs_massf = self.gas.massf
        for i in range(len(gs_massf)):
            massf[config.overall_species_index[self.gmodel_id][i]] = gs_massf[i]
        for mf in massf: fp.write(' %e' % mf)
        Tmodes = [0.0,]*config.overall_modes_count
        gs_Tmodes = self.gas.T_modes
        for i in range(len(gs_Tmodes)):
            Tmodes[config.overall_modes_index[self.gmodel_id][i]] = gs_Tmodes[i]
        for Tmode in Tmodes: fp.write(' %e' % Tmode)
        fp.write('\n')
        return

    def write_history_loc_data(self, fp, t, x):
        """
        Write the initial state of the cell at the given x-location.

        It may be that nothing is written, if the x-location is not
        covered by a cell in the slug.
        """
        for j in range(self.ncells):
            if (x >= self.ifxs[j]) and (x <= self.ifxs[j+1]):
                self.write_initial_flow_data_for_any_cell(fp, t, x)
                break;
        return

    def write_history_cell_data(self, fp, t, ic):
        """
        Write the initial state of a particular cell.

        We are actually going to assume that all cells
        within the slug have the same data.
        """
        x = 0.5 * (self.ifxs[ic] + self.ifxs[ic+1])
        self.write_initial_flow_data_for_any_cell(fp, t, x)
        return
        
    @property
    def energy(self):
        """
        Returns total energy within slug of gas.
        """
        e = 0.0
        for j in range(self.ncells):
            xmid = 0.5*(self.ifxs[j+1] + self.ifxs[j])
            d, area, K_over_L, Twall, vf, htcf = tube.eval(xmid)
            volume = area * (self.ifxs[j+1] - self.ifxs[j])
            e += volume*self.gas.rho*(self.gas.internal_energy + 0.5*self.vel*self.vel)
        return e

#----------------------------------------------------------------------------

class Piston():
    """
    Contains the information for a piston.
    """

    __slots__ = 'indx', 'label', \
                'mass', 'diam', 'L', 'xL0', 'xR0', 'x0', 'vel0', \
                'massf', 'dmassfdt', 'massfmin', \
                'front_seal_f', 'front_seal_area', \
                'back_seal_f', 'back_seal_area', \
                'p_restrain', 'is_restrain', \
                'with_brakes', 'brakes_on', 'brakes_friction_force', \
                'x_buffer', 'on_buffer', 'ecL', 'ecR'

    def __init__(self, mass, diam, xL0, xR0, vel0,
                 massf=1.0, dmassfdt=0.0, massfmin=1.0e-3,
                 front_seal_f=0.0, front_seal_area=0.0,
                 back_seal_f=0.0, back_seal_area=0.0,
                 p_restrain=0.0, is_restrain=0,
                 with_brakes=False, brakes_on=0, brakes_friction_force=0.0,
                 x_buffer=10.e6, on_buffer = 0,
                 label=""):
        """
        Create a piston with specified properties.
        """
        self.mass = mass
        self.diam = diam
        if xR0 < xL0:
            # We would like the x-values to be increasing to the right
            # but we really don't care if the piston length is zero.
            xL0, xR0 = xR0, xL0
        self.xL0 = xL0
        self.xR0 = xR0
        self.L = xR0 - xL0
        self.x0 = 0.5*(xL0 + xR0)
        self.vel0 = vel0
        self.massf = massf
        self.dmassfdt = dmassfdt
        self.massfmin = massfmin
        self.front_seal_f = front_seal_f
        self.front_seal_area = front_seal_area
        self.back_seal_f = back_seal_f
        self.back_seal_area = back_seal_area
        self.p_restrain = p_restrain
        self.is_restrain = is_restrain
        self.with_brakes = with_brakes
        self.brakes_on = brakes_on
        self.brakes_friction_force = brakes_friction_force
        self.x_buffer = x_buffer
        self.on_buffer = on_buffer
        #
        # Connections to boundary conditions will be made later,
        # the gas path is assembled.
        self.ecL = None
        self.ecR = None
        #
        # The Piston objects need an identity that can be
        # transferred to the main simulation program.
        global pistonList
        self.indx = len(pistonList) # next available index
        if len(label) > 0:
            self.label = label
        else:
            # Construct a simple label.
            self.label = 'piston-' + str(self.indx)
        pistonList.append(self)
        return

    def write_config(self, fp):
        """
        Writes the configuration data in JSON format.
        """
        fp.write('"piston_%d": {\n' % self.indx)
        fp.write('  "label": %s,\n' % json.dumps(self.label))
        fp.write('  "mass": %e,\n' % self.mass)
        fp.write('  "diameter": %e,\n' % self.diam)
        fp.write('  "length": %e,\n' % self.L)
        fp.write('  "dmassfdt": %e,\n' % self.dmassfdt)
        fp.write('  "massfmin": %e,\n' % self.massfmin)
        fp.write('  "front_seal_f": %e,\n' % self.front_seal_f)
        fp.write('  "front_seal_area": %e,\n' % self.front_seal_area)
        fp.write('  "back_seal_f": %e,\n' % self.back_seal_f)
        fp.write('  "back_seal_area": %e,\n' % self.back_seal_area)
        fp.write('  "p_restrain": %e,\n' % self.p_restrain)
        fp.write('  "x_buffer": %e,\n' % self.x_buffer)
        fp.write('  "with_brakes": %s,\n' % json.dumps(self.with_brakes))
        fp.write('  "brakes_friction_force": %e,\n' % self.brakes_friction_force)
        # It may be that the piston has only one face connected.
        ecindx = -1
        if self.ecL: ecindx = self.ecL.ecindx
        fp.write('  "ecL_id": %d,\n' % ecindx)
        ecindx = -1
        if self.ecR: ecindx = self.ecR.ecindx
        fp.write('  "ecR_id": %d\n' % ecindx)
        # Note no comma after last item in JSON dict.
        fp.write('},\n') # presume that this dict not the last
        return

    def write_data(self, fp, tindx, write_header):
        """
        Write state data.
        """
        if write_header:
            fp.write("# tindx  x  vel  is_restrain  brakes_on  on_buffer  massf\n")
        fp.write("%d %e %e %d %d %d %e\n" %
                 (tindx, self.x0, self.vel0,
                  self.is_restrain, self.brakes_on, self.on_buffer,
                  self.massf))
        return

    @property
    def energy(self):
        """
        Returns kinetic energy.
        """
        return 0.5*self.mass*self.massf*self.vel0*self.vel0

#----------------------------------------------------------------------------

class Valve():
    """
    A valve may be used to retard the motion of the nearest interface
    that lies between gas cells in a slug.

    Note that we do not include valves into the gas-flow path.
    They are a bit like the tube object in that they influence or interfere
    with the gas flow from the edges.

    We may construct several of these objects somewhere in the input script
    and the information will be stored for later writing into the config file.
    """
    __slots__ = 'indx', 'label', 'x', 'times', 'fopen'

    def __init__(self, x, times=[0.0,], fopen=[1.0,], label=""):
        """
        x: location of the valve
        times: list of time values (in seconds)
        fopen: list of fraction-open values (in range 0.0 to 1.0)

        It may be convenient to generate these lists with a function
        inside your Python input script.
        """
        self.indx = len(valveList)
        self.x = x
        self.times = list(times)
        self.fopen = list(fopen)
        self.label = label
        valveList.append(self)
        return

    def write_config(self, fp):
        """
        Writes the configuration data in JSON format.
        """
        fp.write('"valve_%d": {\n' % self.indx)
        fp.write('  "label": %s,\n' % json.dumps(self.label))
        fp.write('  "x": %e,\n' % self.x)
        fp.write('  "times": %s,\n' % json.dumps(self.times))
        fp.write('  "fopen": %s\n' % json.dumps(self.fopen))
        # Note no comma after last item in JSON dict.
        fp.write('},\n') # presume that this dict not the last
        return

def add_valve(x, times=[0.0,], fopen=[1.0,], label=""):
    """
    Construct and store a Valve object.
    This function is provided just for neatness in the input script.
    """
    _ = Valve(x, times, fopen, label)
    return len(valveList)

#----------------------------------------------------------------------------

class EndCondition():
    """
    Base class for the particular end-(boundary)-conditions defined below.

    This class holds just the connection data.
    """
    __slots__ = 'ecindx', \
                'slugL', 'slugL_end', \
                'slugR', 'slugR_end', \
                'pistonL', 'pistonL_face', \
                'pistonR', 'pistonR_face'
    def __init__(self, slugL=None, slugL_end='R',
                 slugR=None, slugR_end='L',
                 pistonL=None, pistonL_face='R',
                 pistonR=None, pistonR_face='L'):
        # The following may be reassigned during assembly.
        self.slugL = slugL
        self.slugL_end = slugL_end
        self.slugR = slugR
        self.slugR_end = slugR_end
        self.pistonL = pistonL
        self.pistonL_face = pistonL_face
        self.pistonR = pistonR
        self.pistonR_face = pistonR_face
        #
        # The EndCondition objects need an identity that can be
        # transferred to the main simulation program.
        global ecList
        self.ecindx = len(ecList) # next available index
        ecList.append(self)
        return

    def json_str(self):
        myDict = {}
        indx = -1
        if self.slugL: indx = self.slugL.indx
        myDict["left-slug-id"] = indx
        myDict["left-slug-end-id"] = self.slugL_end
        indx = -1
        if self.slugR: indx = self.slugR.indx
        myDict["right-slug-id"] = indx
        myDict["right-slug-end-id"] = self.slugR_end
        indx = -1
        if self.pistonL: indx = self.pistonL.indx
        myDict["left-piston-id"] = indx
        myDict["left-piston-face-id"] = self.pistonL_face
        indx = -1
        if self.pistonR: indx = self.pistonR.indx
        myDict["right-piston-id"] = indx
        myDict["right-piston-face-id"] = self.pistonR_face
        return json.dumps(myDict)


class Diaphragm(EndCondition):
    """
    Contains the information for a diaphragm which controls the
    interaction of two GasSlug objects.
    """

    __slots__ = 'indx', \
                'x0', 'p_burst', 'state', \
                'dt_hold', 'dxL', 'dxR', 'label'

    def __init__(self, x0, p_burst, state=0, dt_hold=0.0,
                 dxL=0.0, dxR=0.0, label="",
                 slugL=None, slugL_end='R',
                 slugR=None, slugR_end='L'):
        """
        Creates a diaphragm with specified properties.

        state == 0, closed
        state == 1, triggered by over-pressure
        state == 2, open

        The connections to GasSlugs are made later via the function
        assemble_gas_path.
        """
        super().__init__(slugL=slugL, slugL_end=slugL_end,
                         slugR=slugR, slugR_end=slugR_end)
        self.x0 = x0
        self.p_burst = p_burst
        self.state = state
        self.dt_hold = dt_hold
        self.dxL = dxL
        self.dxR = dxR
        global diaphragmList
        self.indx = len(diaphragmList) # next available index
        if len(label) > 0:
            self.label = label
        else:
            self.label = "diaphragm-" + str(self.indx)
        diaphragmList.append(self)
        return

    def write_config(self, fp):
        """
        Writes the diaphragm information to the specified file.
        """
        fp.write('"end_condition_%d": {\n' % self.ecindx)
        fp.write('  "class": %s,\n' % json.dumps(self.__class__.__name__))
        fp.write('  "label": %s,\n' % json.dumps(self.label))
        fp.write('  "x0": %e,\n' % self.x0)
        fp.write('  "p_burst": %e,\n' % self.p_burst)
        fp.write('  "dt_hold": %e,\n' % self.dt_hold)
        fp.write('  "dxL": %e,\n' % self.dxL)
        fp.write('  "dxR": %e,\n' % self.dxR)
        fp.write('  "connections": %s\n' % self.json_str()) # no comma for last item
        fp.write('},\n')
        return

    def write_data(self, fp, tindx, write_header):
        """
        Write state data.
        """
        if write_header:
            fp.write("# tindx state \n")
        fp.write("%d %d\n" % (tindx, self.state))
        return


class GasInterface(EndCondition):
    """
    Contains the information for an interface between two slugs.

    The primary use of this class is to locate the ends of the connected GasSlugs.
    Implicitly, the logical connections are also made via the assemble_gas_path().
    """

    __slots__ = 'x0'

    def __init__(self, x0,
                 slugL=None, slugL_end='R',
                 slugR=None, slugR_end='L'):
        """
        Creates as interface between two gas slugs at specified location.

        :param x0: (float) Initial position, in metres.
        """
        super().__init__(slugL=slugL, slugL_end=slugL_end,
                         slugR=slugR, slugR_end=slugR_end)
        self.x0 = x0
        global interfaceList
        interfaceList.append(self)
        return

    def write_config(self, fp):
        """
        Write in JSON format.
        """
        fp.write('"end_condition_%d": {\n' % self.ecindx)
        fp.write('  "class": %s,\n' % json.dumps(self.__class__.__name__))
        fp.write('  "x0": %e,\n' % self.x0)
        fp.write('  "connections": %s\n' % self.json_str()) # no comma for last item
        fp.write('},\n')
        return


class FreeEnd(EndCondition):
    """
    Locates the end of a gas slug.
    """

    __slots__ = 'x0'

    def __init__(self, x0,
                 slugL=None, slugL_end='R',
                 slugR=None, slugR_end='L'):
        """
        Creates a GasSlug end-condition with a specified location.

        :param x0: (float) Initial position, in metres.
        """
        super().__init__(slugL=slugL, slugL_end=slugL_end,
                         slugR=slugR, slugR_end=slugR_end)
        self.x0 = x0
        global freeEndList
        freeEndList.append(self)
        return

    def write_config(self, fp):
        """
        Write in JSON format.
        """
        fp.write('"end_condition_%d": {\n' % self.ecindx)
        fp.write('  "class": %s,\n' % json.dumps(self.__class__.__name__))
        fp.write('  "x0": %e,\n' % self.x0)
        fp.write('  "connections": %s\n' % self.json_str()) # no comma for last item
        fp.write('},\n')
        return


class VelocityEnd(EndCondition):
    """
    Fixed-velocity end condition for a GasSlug.

    Velocity may be zero.
    """

    __slots__ = 'x0', 'vel'

    def __init__(self, x0, vel=0.0,
                 slugL=None, slugL_end='R',
                 slugR=None, slugR_end='L'):
        """
        Creates a GasSlug end-condition with a specified location and velocity.

        :param x0: (float) Initial position, in metres.
        :param v: (float) Velocity, in m/s, of the end-point of the GasSlug.
        """
        super().__init__(slugL=slugL, slugL_end=slugL_end,
                         slugR=slugR, slugR_end=slugR_end)
        self.x0 = x0
        self.vel = vel
        global velocityEndList
        velocityEndList.append(self)
        return

    def write_config(self, fp):
        """
        Write in JSON format.
        """
        fp.write('"end_condition_%d": {\n' % self.ecindx)
        fp.write('  "class": %s,\n' % json.dumps(self.__class__.__name__))
        fp.write('  "x0": %e,\n' % self.x0)
        fp.write('  "vel": %e,\n' % self.vel)
        fp.write('  "connections": %s\n' % self.json_str()) # no comma for last item
        fp.write('},\n')
        return


class PistonFace(EndCondition):
    """
    Just connection data between gas slug and piston.
    """

    __slots__ = 'x0'

    def __init__(self,
                 slugL=None, slugL_end='R',
                 slugR=None, slugR_end='L',
                 pistonL=None, pistonL_face='R',
                 pistonR=None, pistonR_face='L'):
        """
        Creates a GasSlug end-condition at a piston face.
        """
        if (not pistonL) and (not pistonR):
            raise Exception("Cannot have two pistons attached to a PistonFace.")
        if (not slugL) and (not slugR):
            raise Exception("Cannot have two gas slugs attached to a PistonFace.")
        if (pistonL):
            self.x0 = pistonL.xR0
        elif (pistonR):
            self.x0 = pistonR.xL0
        else:
            raise Exception("Could not find a position of a piston face.")
        super().__init__(slugL=slugL, slugL_end=slugL_end,
                         slugR=slugR, slugR_end=slugR_end,
                         pistonL=pistonL, pistonL_face=pistonL_face,
                         pistonR=pistonR, pistonR_face=pistonR_face)
        global pistonFaceList
        pistonFaceList.append(self)
        return

    def write_config(self, fp):
        """
        Write in JSON format.
        """
        fp.write('"end_condition_%d": {\n' % self.ecindx)
        fp.write('  "class": %s,\n' % json.dumps(self.__class__.__name__))
        fp.write('  "x0": %e,\n' % self.x0)
        fp.write('  "connections": %s\n' % self.json_str()) # no comma for last item
        fp.write('},\n')
        return

# --------------------------------------------------------------------

def assemble_gas_path(*components):
    """
    Assembles a gas path by making the logical connections between
    adjacent components.

    The components are assembled left-to-right,
    as they are supplied to this function.

    components: An arbitrary number of arguments representing
        individual components or lists of components.
        Each component may be a GasSlug, Piston, or any
        other gas-path object, however, it doesn't always make sense
        to connect arbitrary components.
        For example, connecting a GasSlug to a Piston is reasonable
        but connecting a Piston to a Diaphragm without an intervening
        GasSlug does not make sense in the context of this simulation
        program.
    """
    print("Assemble gas path:")
    clist = []
    for c in components:
        if isinstance(c,tuple) or isinstance(c,list):
            clist.extend(c)
        else:
            clist.append(c)
    for i in range(len(clist)-1):
        connect_pair(clist[i], clist[i+1])
    # Now that we have made all connections,
    # we can locate the ends of the gas slugs.
    for c in components:
        if isinstance(c, GasSlug):
            c.xL = c.ecL.x0
            c.xR = c.ecR.x0
    return


def connect_pair(cL, cR):
    """
    Make the logical connection between a pair of components.

    cL: is left object
    cR: is right object

    Usually called by assemble_gas_path.
    """
    print("connect_pair L:", cL.__class__.__name__,
          " R:", cR.__class__.__name__)

    if isinstance(cL,VelocityEnd) and isinstance(cR, GasSlug):
        cR.ecL = cL
        cL.slugR = cR
        cL.slugR_end = 'L'
        print("  velocity-end <--> gas-slug is done")
    elif isinstance(cL,FreeEnd) and isinstance(cR, GasSlug):
        cR.ecL = cL
        cL.slugR = cR
        cL.slugR_end = 'L'
        print("  free-end <--> gas-slug is done")
    elif isinstance(cL,GasInterface) and isinstance(cR, GasSlug):
        cR.ecL = cL
        cL.slugR = cR
        cL.slugR_end = 'L'
        print("  gas-interface <--> gas-slug is done")
    elif isinstance(cL,Piston) and isinstance(cR, GasSlug):
        print("  constructing a new PistonFace")
        pf = PistonFace(pistonL=cL, pistonL_face='R',
                        slugR=cR, slugR_end='L')
        cL.ecR = pf
        cR.ecL = pf
        print("  piston <--> piston-face <--> gas-slug connections done")
    elif isinstance(cL,PistonFace) and isinstance(cR, GasSlug):
        cL.slugR = cR
        cL.slugR_end = 'L'
        cR.ecL = cL
        print("  piston-face <--> gas-slug is done")
    elif isinstance(cL,Diaphragm) and isinstance(cR, GasSlug):
        cL.slugR = cR
        cL.slugR_end = 'L'
        cR.ecL = cL
        print("  diaphragm <--> gas-slug is done")
    elif isinstance(cL,GasSlug) and isinstance(cR, VelocityEnd):
        cL.ecR = cR
        cR.slugL = cL
        cR.slugL_end = 'R'
        print("  gas-slug <--> velocity-end is done")
    elif isinstance(cL,GasSlug) and isinstance(cR, FreeEnd):
        cL.ecR = cR
        cR.slugL = cL
        cR.slugL_end = 'R'
        print("  gas-slug <--> free-end is done")
    elif isinstance(cL,GasSlug) and isinstance(cR, GasInterface):
        cL.ecR = cR
        cR.slugL = cL
        cR.slugL_end = 'R'
        print("  gas-slug <--> gas-interface is done")
    elif isinstance(cL,GasSlug) and isinstance(cR, Piston):
        print("  constructing a new PistonFace")
        pf = PistonFace(pistonR=cR, pistonR_face='L',
                        slugL=cL, slugL_end='R')
        cL.ecR = pf
        cR.ecL = pf
        print("  gas-slug <--> piston-face <--> piston connections done")
    elif isinstance(cL,GasSlug) and isinstance(cR, PistonFace):
        cL.ecR = cR
        cR.slugL = cL
        cR.slugL_end = 'R'
        print("  gas-slug <--> piston-face is done")
    elif isinstance(cL,GasSlug) and isinstance(cR, Diaphragm):
        cL.ecR = cR
        cR.slugL = cL
        cR.slugL_end = 'R'
        print("  gas-slug <--> diaphragm is done")
    else:
        raise Exception("  Invalid pair to connect.")
    return


# --------------------------------------------------------------------

def write_initial_files():
    """
    Writes the files needed for the main simulation code.

    These files are found in the directory config.job_name.
    """
    global config, pistonList, diaphragmList, slugList
    print("Begin write initial files.")
    if not os.path.exists(config.job_name):
        os.mkdir(config.job_name)
    #
    fp = open(config.job_name+'/config.json', 'w')
    fp.write("{\n")
    config.write(fp)
    for slug in slugList: slug.write_config(fp)
    for piston in pistonList: piston.write_config(fp)
    for valve in valveList: valve.write_config(fp)
    for ec in ecList: ec.write_config(fp)
    # Previous entries all presume that they are not the last.
    fp.write('"dummy_item": 0\n')
    fp.write('}\n')
    fp.close()
    #
    tube.set_up_tables()
    fp = open(config.job_name+'/tube.data', 'w')
    tube.write(fp)
    fp.close()
    #
    for slug in slugList:
        slug.construct_cells_and_faces()
        fileName = config.job_name + ('/slug-%04d-faces.data' % slug.indx)
        fp = open(fileName, 'w')
        slug.write_face_data(fp, 0, True)
        fp.close()
        fileName = config.job_name + ('/slug-%04d-cells.data' %  slug.indx)
        fp = open(fileName, 'w')
        slug.write_cell_data(fp, 0, True)
        fp.close()
    #
    for piston in pistonList:
        fileName = config.job_name + ('/piston-%04d.data' % piston.indx)
        fp = open(fileName, 'w')
        piston.write_data(fp, 0, True)
        fp.close()
    #
    for diaphragm in diaphragmList:
        fileName = config.job_name + ('/diaphragm-%04d.data' % diaphragm.indx)
        fp = open(fileName, 'w')
        diaphragm.write_data(fp, 0, True)
        fp.close()
    #
    fileName = config.job_name + '/times.data'
    fp = open(fileName, 'w')
    fp.write('# tindx time\n')
    fp.write('%d %e\n' % (0, 0.0))
    fp.close()
    #
    history_header_line = '# 1:t  2:x  3:vel  4:L_bar  5:rho  6:p  7:T  8:u  9:a  10:shear_stress  11:heat_flux'
    column_count = 12
    for gi in range(len(config.gmodels)):
        species_names = config.gmodels[gi].species_names
        for i in range(config.gmodels[gi].n_species):
            history_header_line += ('  %d:gas-%d-massf-%s' % (column_count, gi, species_names[i]))
            column_count += 1
    for gi in range(len(config.gmodels)):
        for i in range(config.gmodels[gi].n_modes):
            history_header_line += ('  %d:gas-%d-Tmode[%d]' % (column_count, gi, i))
            column_count += 1
    history_header_line += '\n'
    #
    for slug in slugList:
        for ic in slug.hcells:
            fileName = config.job_name + ('/history-cell-%04d-in-slug-%04d.data' % (ic, slug.indx))
            fp = open(fileName, 'w')
            fp.write(history_header_line)
            slug.write_history_cell_data(fp, 0.0, ic)
            fp.close()
    #
    for ih in range(len(config.hloc_list)):
        fileName = config.job_name + ('/history-loc-%04d.data' % ih)
        fp = open(fileName, 'w')
        fp.write(history_header_line)
        for slug in slugList: slug.write_history_loc_data(fp, 0.0, config.hloc_list[ih])
        fp.close()
    #
    fileName = config.job_name + '/energies.data'
    fp = open(fileName, 'w')
    fp.write('# time')
    for i in range(len(slugList)): fp.write(' slug-%d' % i)
    for i in range(len(pistonList)): fp.write(' piston-%d' % i)
    fp.write(' total\n')
    fp.write('%e' % 0.0)
    e_total = 0.0
    for s in slugList:
        e = s.energy
        fp.write(' %e' % e)
        e_total += e
    for p in pistonList:
        e = p.energy
        fp.write(' %e' % e)
        e_total += e
    fp.write(' %e\n' % e_total)
    fp.close()
    #
    fileName = config.job_name + '/events.txt'
    fp = open(fileName, 'w')
    fp.write('# time-stamped event descriptions\n')
    fp.close()
    #
    print("End write initial files.")
    return


# --------------------------------------------------------------------

if __name__ == '__main__':
    print("Begin l1d4-prep...")

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
        print("Summary of components:")
        print("  gas slugs         :", len(slugList))
        print("  pistons           :", len(pistonList))
        print("  valves            :", len(valveList))
        print("  diaphragms        :", len(diaphragmList))
        print("  free-ends         :", len(freeEndList))
        print("  velocity-ends     :", len(velocityEndList))
        print("  gas-gas interfaces:", len(interfaceList))
        print("  piston faces      :", len(pistonFaceList))
        if len(slugList) < 1:
            print("Warning: no gas slugs defined; this is unusual.")
        write_initial_files()
    print("Done.")
    sys.exit(0)
