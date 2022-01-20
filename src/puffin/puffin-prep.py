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
                'dx', 'cfl', 'cfl_count', \
                'print_count', 'plot_count', \
                'max_x', 'max_step', 'x_order'

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
        self.dx = 1.0e-3
        self.cfl = 0.5
        self.cfl_count = 10;
        self.print_count = 50;
        self.plot_count = 10;
        self.max_x = 1.0
        self.max_step = 10
        self.x_order = 2
        #
        GlobalConfig.count += 1
        return

    def write(self, fp):
        """
        Writes the configuration data to the specified file in JSON format.

        fp: handle to an open file.
        """
        global config, streamTubeList
        #
        fp.write('"config": {\n')
        fp.write('  "title": "%s",\n' % self.title)
        fp.write('  "gas_model_file": %s,\n' % self.gas_model_file)
        fp.write('  "reaction_file_1": %s,\n' % self.reaction_file_1)
        fp.write('  "reaction_file_2": %s,\n' % self.reaction_file_2)
        fp.write('  "reacting": %s,\n' % json.dumps(self.reacting))
        fp.write('  "T_frozen": %e,\n' % self.T_frozen)
        fp.write('  "max_x": %e,\n' % self.max_x)
        fp.write('  "max_step": %d,\n' % self.max_step)
        fp.write('  "dx": %e,\n' % self.dx)
        fp.write('  "cfl": %e,\n' % self.cfl)
        fp.write('  "cfl_count": %d,\n' % self.cfl_count)
        fp.write('  "print_count": %d,\n' % self.print_count)
        fp.write('  "plot_count": %d,\n' % self.plot_count)
        fp.write('  "x_order": %d,\n' % self.x_order)
        #
        fp.write('  "nstreams": %d\n' % len(streamTubeList))
        # Note, no comma after last item inside JSON dict
        fp.write('},\n') # comma here because not the last item in file
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
    if not os.path.exists(fileName):
        raise Exception("Gas model file not found: " + fileName)
    gmodel = GasModel(fileName)
    gmodel_id = len(config.gmodels)
    config.gmodel = gmodel
    config.gas_model_file = fileName
    nsp = gmodel.n_species
    nmodes = gmodel.n_modes
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
    #
    return gmodel_id


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
                'gas', 'velx', 'vely', \
                'ncells', 'y0', 'y1'

    def __init__(self, gas=None, velx=1000.0, vely=0.0,
                 y0=None, y1=None, ncells=10, label="",):
        """
        Creates an outline of a streamtube with user-specified boundaries.

        gas : inflow gas state
        velx : inflow x-velocity
        vely : inflow y-velocity
        y0 : lower boundary y as a function of x
        y1 : upper boundary y as a function of x
        ncells : number of cells between lower and upper boundary
        label : optional label for postprocessing
        """
        global config
        #
        # The StreamTube objects need an identity that can be
        # transferred to the main calculation program.
        global streamTubeList
        self.indx = len(streamTubeList) # next available index
        streamTubeList.append(self)
        #
        return

    def write(fp):
        return

    def construct_boundaries():
        """
        Constructs inflow gas data, lower and upper boundary paths.
        """
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
    with open(config.job_name+'/config.json', 'w') as fp:
        fp.write("{\n")
        config.write(fp)
        # Previous entries all presume that they are not the last.
        fp.write('"dummy_item": 0\n')
        fp.write('}\n')
    #
    for st in streamTubeList:
        st.construct_boundaries()
        fileName = config.job_name + ('/streamtube-%04d.data' % st.indx)
        with open(fileName, 'w') as fp:
            st.write(fp)
    #
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
