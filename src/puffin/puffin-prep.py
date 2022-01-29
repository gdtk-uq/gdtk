#! /usr/bin/env python3
"""
Python program to set up a simulation for the steady-state 2D Flow Solver.

It is intended for the user to define their particular flow path and
inflow in terms of the data objects defined in this module.
As part of its initialization, this program will read and execute
a user-specified input script that contains, in Python,
the code that defines the gas-path details.

Usage:
  $ puffin-prep --job=<jobName>

The simulation control data is then organised via the classes
GlobalConfig and StreamTube.
These classes provide places to store the configuration information
and their function/method names appear as commands in the user's
job description file.

Note that config is a reference to the global configuration information
describing the calculation.  There is one such variable set up by the
main program and the user's script should directly set the attributes
of this variable to adjust settings for the simulation.

Author:
  P.A. Jacobs

Versions:
  2022-01-21  First Python code
"""

# ----------------------------------------------------------------------
#
import sys
import os
import collections
from getopt import getopt
import math
import numpy as np
import json
from eilmer.gas import GasModel, GasState, ThermochemicalReactor

sys.path.append("") # so that we can find user's scripts in current directory

shortOptions = "hf:"
longOptions = ["help", "job="]

def printUsage():
    print("")
    print("Usage: puffin-prep" + \
          " [--help | -h]" + \
          " [--job=<jobName> | -f <jobName>]")
    print("")
    return

#----------------------------------------------------------------------
# This is where we store the core data for the simulation.

FluxCalc = collections.namedtuple('_', 'hanel riemann')(*range(2))

class GlobalConfig(object):
    """
    Contains the calculation control parameters.

    The user's script should not create one of these
    but should specify the simulation parameters by
    altering the attributes of the global object "config"
    that already exists by the time the user's script executes.
    """
    count = 0

    # We want to prevent the user's script from introducing new attributes
    # via typographical errors.
    __slots__ = 'job_name', 'title', \
                'gas_model_file', 'gmodel', \
                'reaction_file_1', 'reaction_file_2', 'reactor', 'reacting', 'T_frozen', \
                'axisymmetric', \
                'dx', 'print_count', 'plot_dx', 'max_x', 'max_step', \
                'cfl', 'max_step_relax', 'x_order', 't_order', 'flux_calc'

    def __init__(self):
        """Accepts user-specified data and sets defaults. Make one only."""
        if GlobalConfig.count >= 1:
            raise Exception("Already have a GlobalConfig object defined.")

        self.job_name = ""
        self.title = "Another Puffin calculation."
        self.gas_model_file = ""
        self.gmodel = None
        self.reaction_file_1 = ""
        self.reaction_file_2 = ""
        self.reactor = None
        self.reacting = False
        self.T_frozen = 300.0
        self.axisymmetric = False
        self.dx = 1.0e-3
        self.print_count = 50
        self.plot_dx = 1.0e-2
        self.max_x = 1.0
        self.max_step = 10
        self.cfl = 0.5
        self.max_step_relax = 100
        self.x_order = 2
        self.t_order = 2
        self.flux_calc = FluxCalc.hanel
        #
        GlobalConfig.count += 1
        return

    def write(self):
        """
        Writes the configuration data to the specified file in JSON format.

        fp: handle to an open file.
        """
        global config, streamTubeList
        #
        fp = open(config.job_name+'/config.json', 'w')
        fp.write('{\n')
        fp.write('  "title": "%s",\n' % self.title)
        fp.write('  "gas_model_file": "%s",\n' % self.gas_model_file)
        fp.write('  "reaction_file_1": "%s",\n' % self.reaction_file_1)
        fp.write('  "reaction_file_2": "%s",\n' % self.reaction_file_2)
        fp.write('  "reacting": %s,\n' % json.dumps(self.reacting))
        fp.write('  "T_frozen": %e,\n' % self.T_frozen)
        fp.write('  "axisymmetric": %s,\n' % json.dumps(self.axisymmetric))
        fp.write('  "dx": %e,\n' % self.dx)
        fp.write('  "print_count": %d,\n' % self.print_count)
        fp.write('  "plot_dx": %e,\n' % self.plot_dx)
        fp.write('  "max_x": %e,\n' % self.max_x)
        fp.write('  "max_step": %d,\n' % self.max_step)
        fp.write('  "cfl": %e,\n' % self.cfl)
        fp.write('  "max_step_relax": %d,\n' % self.max_step_relax)
        fp.write('  "x_order": %d,\n' % self.x_order)
        fp.write('  "t_order": %d,\n' % self.t_order)
        fp.write('  "flux_calc": %d,\n' % self.flux_calc)
        #
        fp.write('  "n_streams": %d,\n' % len(streamTubeList))
        for st in streamTubeList:
            # Assemble a dictionary defining the flowstate.
            fs = {'p':st.gas.p, 'T':st.gas.T, 'massf':st.gas.massf,
                  'velx':st.velx, 'vely':st.vely}
            fp.write('  "inflow_%d": %s,\n' % (st.indx, json.dumps(fs)))
            fp.write('  "ncells_%d": %d,\n' % (st.indx, st.ncells))
        # Note, no comma after last item inside JSON dict
        fp.write('  "dummy_item": 0\n')
        fp.write('}\n')
        fp.close()
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


# --------------------------------------------------------------------

streamTubeList = []


class StreamTube():
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
                'gas', 'velx', 'vely', 'ncells', \
                'y0', 'y1', 'bc0', 'bc1', \
                'xs', 'y0s', 'y1s', 'bc0s', 'bc1s', 'nbc'

    def __init__(self, gas=None, velx=1000.0, vely=0.0,
                 y0=None, y1=None, bc0=None, bc1=None,
                 ncells=10, nbc=101, label="",):
        """
        Creates an outline of a streamtube with user-specified boundaries.

        gas    : inflow gas state
        velx   : inflow x-velocity
        vely   : inflow y-velocity
        y0     : lower boundary y as a function of x
                 The boundaries extend from x=0.0 to x=config.max_x.
        y1     : upper boundary y as a function of x
        bc0    : lower boundary condition as an integer function of x
                 if bc0(x) == 0, this stream exchanges data with lower stream
                 if bc0(x) == 1, this stream exchanges data with upper stream
        bc1    : upper boundary condition as an integer function of x
        ncells : number of cells between lower and upper boundary
        nbc    : number of points in the boundary definitions
        label  : optional label for postprocessing
        """
        global config
        self.gas = GasState(config.gmodel)
        self.gas.copy_values(gas)
        self.velx = velx
        self.vely = vely
        if callable(y0):
            self.y0 = y0
        else:
            raise RuntimeError("Was expecting a function for y0")
        if callable(y1):
            self.y1 = y1
        else:
            raise RuntimeError("Was expecting a function for y1")
        if callable(bc0):
            self.bc0 = bc0
        else:
            raise RuntimeError("Was expecting a function for bc0")
        if callable(bc1):
            self.bc1 = bc1
        else:
            raise RuntimeError("Was expecting a function for bc1")
        self.ncells = ncells
        self.nbc = nbc
        #
        # The StreamTube objects need an identity that can be
        # transferred to the main calculation program.
        global streamTubeList
        self.indx = len(streamTubeList) # next available index
        streamTubeList.append(self)
        #
        return

    def construct_and_write_boundaries(self):
        """
        Constructs lower and upper boundary paths.

        Note that some of the information is in the global config.
        """
        global config
        xs = np.linspace(0.0, config.max_x, self.nbc)
        y0s = [self.y0(x) for x in xs]
        y1s = [self.y1(x) for x in xs]
        bc0s = [self.bc0(x) for x in xs]
        bc1s = [self.bc1(x) for x in xs]
        fileName = config.job_name + ('/streamtube-%d.data' % self.indx)
        with open(fileName, 'w') as fp:
            fp.write('# x  y0  y1  bc0  bc1\n')
            fp.write('# n_bc=%d\n' % self.nbc)
            fp.write('# dx_bc=%g\n' % (xs[1]-xs[0]))
            for i in range(len(xs)):
                fp.write('%g %g %g %g %g\n' % (xs[i], y0s[i], y1s[i], bc0s[i], bc1s[i]))
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
    config.write()
    for st in streamTubeList:
        st.construct_and_write_boundaries()
    print("End write initial files.")
    return

# --------------------------------------------------------------------

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
        config.job_name, ext = os.path.splitext(jobName)
        inputScriptName = config.job_name + ".py"
        print("Input script file name: %s" % inputScriptName)

        # The user-specified input comes in the form of Python code.
        # It is up to the user to be careful; there is no security.
        exec(compile(open(inputScriptName, "rb").read(), inputScriptName, 'exec'))
        if len(streamTubeList) < 1:
            print("Warning: no stream tubes defined; this is unusual.")
        write_initial_files()
    print("Done.")
    sys.exit(0)
