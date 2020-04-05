#! /usr/bin/env python3
"""
Python program to set up a simulation for the Lagrangian 1D Flow Solver.

It is intended for the user to define their particular facility and
flow in terms of the data objects defined in this module.  As part of
its initialization, this program will execute a user-specified input 
script that contains, in Python, the user's code that defines both
facility geometry and gas-path details.

Usage::

    $ l1d4-prep -f <job>

The simulation control data is then organised via the classes:
GlobalData, GasSlug, Piston and Diaphragm.  These classes
provide places to store the configuration information and their
function/method names appear as commands in the user's
job description file.

When setting up a new simulation, first define the tube as a set
of (x,d) break-points and identify regions of head-loss and
regions where the wall-temperature varies from the nominal value.
Create the GasSlugs, Pistons, and Diaphragms that will make up the
gas path.  Note that places where two GasSlugs join will need a
GasInterface to be defined.  Once all of the components have been
created, assemble the gas path and then set any of the time-stepping
parameters for which you want values other than the default.

Here is an example script for the Sod shock-tube problem::

    # sod.py
    config.title = 'Sods ideal shock tube, 2020-04-04'
    my_gm = add_gas_model('ideal-air-gas-model.lua')
    
    # Define the tube walls.
    add_break_point(0.0, 0.01)
    add_break_point(3.0, 0.01)
    
    # Create the gas-path.
    left_wall = VelocityEnd(x0=0.0, vel=0.0)
    driver_gas = GasSlug(p=100.0e3, vel=0.0, T=348.4, gmodel_id=my_gm, nn=200)
    interface = GasInterface(x0=0.5)
    driven_gas = GasSlug(p=10.0e3, vel=0.0, T=278.7, gmodel_id=my_gm, nn=100)
    right_wall = VelocityEnd(x0=1.0, vel=0.0)
    assemble_gas_path(left_wall, driver_gas, interface, driven_gas, right_wall)
    
    # Set some time-stepping parameters
    config.dt_init = 1.0e-7
    config.max_time = 0.6e-3
    config.max_step = 5000
    add_dt_plot(0.0, 10.0e-6, 5.0e-6)
    add_history_loc(0.7)

This script should define the gas path::

  .       |+----- driver-gas -----+|+----- driven-gas -----+|
  .       |                        |                        |
  .       |                        |                        |
  .   left-wall                interface               right-wall
    
and can be invoked with the command::

    $ l1d4-prep -f sod

Upon getting to the end of the user's script, this program should then write 
(1) a complete simulation parameter file (config/sod.config) in JSON format
(2) A tube-definition file.
(3) State files for pistons, diaphragms and gas slugs. 

Note that Python is very picky about whitespace.  If you cut and paste the
example from above, make sure that the lines start in the first column and
that indentation is consistent with Python's syntax rules.

Globally-defined object
-----------------------

* gconfig: Contains the Global Configuration information describing the simulation.
  Note that there is one such variable set up by the main program and
  the user's script should directly set the attributes of this variable
  to adjust settings for the simulation.

.. Author: P.A. Jacobs

.. Versions: 
   2005-Jun
   2006-Jul-24 Ported to Rowan's new C++ chemistry.
   2010-Mar-May More reworking for L1d3 and 
                the latest-greatest thermochemistry
   2012-Sep-Oct Much cleaning up and Sphinx docs.
   2020-Apr-04 L1d4 flavour started
"""

# ----------------------------------------------------------------------
#
import sys
import os
from getopt import getopt
import math
import numpy as np
from eilmer.gas import GasModel, GasState, ThermochemicalReactor
from eilmer import roberts

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
    
    The following attributes are available:

    * title: Short title string for embedding in the parameter and solution files.
    
    * gas_file_names: list of file names for the detailed gas models.
      There may be one or more, but you have to specify one when you
      make each GasSlug.

    * reacting: (bool) If set to True, Rowan's finite-rate chemistry will
      be active.  (Default is False)

    * dt_init: (float) The size of the time-step that will be used for the
      first few simulation steps.
      After a few steps, the cfl condition takes over the determination
      of a suitable time-step.

    * max_time: (float) The simulation will stop if it reaches this time.
      It is most usual to use this critereon to stop the simulation.

    * max_step: The simulation will be stopped if it reaches
      this number of steps.
      This is mostly used to catch the problem of the calculation taking
      a very long time (measured by one's patience), possibly because
      the time-step size has decreased to an extremely small value.

    * cfl: (float) Largest allowable CFL number.
      The time step is adjusted to ensure that this value is not exceeded
      in any particular cell.
      A typical value of 0.25 seems to work well for simulations with
      sudden events such as diaphragm bursting, while a value as high as
      0.5 should be considered only for well-behaved flows.

    * t_order: (int) 
      1=Euler time-stepping. This is generally cheap and nasty.
      2=predictor-corrector time-stepping, nominally second order.
      This is the default setting.
      It is, however, twice as CPU intensive as Euler time-stepping.

    * x_order: (int) 
      1=use cell averages without high-order reconstruction.
      Use this only if the second-order calculation is showing problems.
      2=use limited reconstruction (nominally second order).
      This is the default selection. 

    * dt_plot_list: (list of tuples) 
      Specifies the frequency of writing complete solutions
      (for later plotting, maybe) and also for the writing of data at
      history locations.
      It may be convenient to have different frequencies of writing such
      output at different stages of the simulation.
      For example, free-piston driven shock tunnels have a fairly long
      period during which the piston travels the length of the compression
      tube and then a relatively short period, following diaphragm rupture,
      when all the interesting things happen.
      It is good to have low-frequency output during most of the compression
      process and higher-frequency output starting just before diaphragm
      rupture.
      Arranging good values may require some trial and error.
      Add entries to this list via the add_dt_plot function.

    * hloc_list: (list of floats)
      List of x-coordinates for the history locations.
      Add entries via the function add_history_loc.
    """
    count = 0

    # We want to prevent the user's script from introducing new attributes
    # via typographical errors.
    __slots__ = 'job_name', 'title', 'gas_model_files', 'gmodels', \
                'reaction_scheme_files', 'reacting', \
                'dt_init', 'cfl', 'dt_plot_list', \
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
        self.reaction_scheme_files = []
        self.reacting = False
        self.dt_init = 1.0e-6
        self.cfl = 0.5
        # If dt_plot_list is still an empty list when we write
        # the parameter file, just use the max_time value for
        # both dt_plot and dt_his.  Fill in later.
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
        fp.write("{\n")
        fp.write('title = %s\n' % self.title)
        fp.write('gas_model_files = %s\n' % self.gas_model_files)
        fp.write('reaction_scheme_files = %s\n' % self.reaction_scheme_files)
        fp.write('reacting = %s\n' % self.reacting)
        fp.write('max_time = %e\n' % self.max_time)
        fp.write('max_step = %d\n' % self.max_step)
        fp.write('dt_init = %e\n' % self.dt_init)
        fp.write('cfl = %e\n' % self.cfl)
        fp.write('x_order = %d\n' % self.x_order)
        fp.write('t_order = %d\n' % self.t_order)
        #
        if len(config.dt_plot_list) == 0:
            # Since the user did not specify any, default to the end.
            self.add_dt_plot(0.0, config.max_time, config.max_time)
        n_dt_plot = len(self.dt_plot_list)
        fp.write('n_dt_plot = %d\n' % n_dt_plot)
        fp.write('t_change =');
        for i in range(n_dt_plot):
            fp.write(' %e' % config.dt_plot_list[i][0])
        fp.write('\n')
        fp.write('dt_plot =');
        for i in range(n_dt_plot):
            fp.write(' %e' % config.dt_plot_list[i][1])
        fp.write('\n')
        fp.write('dt_his =');
        for i in range(n_dt_plot):
            fp.write(' %e' % config.dt_plot_list[i][2])
        fp.write('\n')
        #
        n_hloc = len(config.hloc_list)
        fp.write('hloc_n = %d\n' % n_hloc)
        fp.write('hloc_x =')
        for i in range(n_hloc):
            fp.write(' %e' % config.hloc_list[i])
        fp.write('\n')
        #
        fp.write('nslug = %d\n' % len(slugList))
        fp.write('npiston = %d\n' % len(pistonList))
        fp.write('ndiaphragm = %d\n' % len(diaphragmList))
        for p in pistonList: p.write_config(fp)
        for d in diaphragmList: d.write_config(fp)
        for s in slugList: s.write_config(fp)
        fp.write('}\n')
        return
    
# We will create just one GlobalConfig object that the user can alter.
config = GlobalConfig()

# --------------------------------------------------------------------
# The following functions are to provide convenient ways of setting
# some of the GlobalConfig elements.

def add_gas_model(fileName):
    """
    Initialize and add a GasModel in the GlobalConfig list of gas models.

    :param FileName: (string) Name of the detailed-gas-model file.
    """
    global config
    gmodel = GasModel(fileName)
    gmodel_id = len(config.gmodels)
    config.gmodels.append(gmodel)
    return gmodel_id


def add_dt_plot(t_change, dt_plot, dt_his):
    """
    Add a dt tuple to the dt_plot tuple list in GlobalConfig.

    :param t_change: (float) The time, in seconds, 
        at which this dt_plot and dt_his should take effect.
    :param dt_plot: (float) Time interval between writing whole solutions
        for later plotting.
    :param dt_his: (float) Time interval between writing data to history file.
    :returns: the new length of the list
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
    Add a location to the history-location list in L{GlobalConfig}.

    :param x: (float) x-coordinate, in metres, of the sample point.
    :returns: Number of sample points defined so far.
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
    
    The following attributes are available:

    * n: (int) The number of small segments that will be used to describe
      the tube's area distribution internal to the simulation.
      To enable a fast lookup process for the area calculation,
      the area variation between equally-spaced x-positions is taken
      to be linear.
      The default value is 4000 and probably won't need to be changed
      except for geometries with rapidly changing cross-sections.

    * xd_list: List of break-point tuples defining the tube wall.
      Add elements to the list via the function add_break_point.

    * T_nominal: (float) The nominal wall temperature (in degrees K)
      in the absence of a patch of differing temperature.

    * T_patch_list: (list of tuples)
      Regions of the tube wall that have temperature different to the 
      nominal value can be specified via the function add_T_patch.

    * loss_region_list: (list of tuples)
      List of head-loss regions, usually associated
      with sudden changes in tube cross-section and diaphragm stations.
      Add regions via the function add_loss_region.
    """
    count = 0

    __slots__ = 'n', 'x_list', 'd_list', 'T_nominal', 'T_patch_list', 'loss_region_list'
    
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
        #
        GlobalConfig.count += 1
        return

    def write(self, fp):
        """
        Writes the tube specification to the specified file, in small steps.
        """
        global tube
        nsegments = len(tube.x_list) - 1
        assert len(tube.d_list) == nsegments+1, "Mismatch in tube.x_list, tube.d_list lengths."
        if nsegments < 1:
            raise Exception("You did not specify at least two points in your tube.")
        xs = np.linspace(tube.x_list[0], tube.x_list[-1], num=tube.n+1)
        ds = []
        for x in xs:
            iseg = nsegments
            for i in range(nsegments):
                if x <= tube.x_list[i]:
                    iseg = i
                    break
            # Have selected the segment containing x.
            frac = (x - tube.x_list[iseg-1]) / (tube.x_list[iseg] - tube.x_list[iseg-1])
            d = (1.0-frac)*tube.d_list[iseg-1] + frac*tube.d_list[iseg]
            ds.append(d)
        #
        K_over_Ls = []
        for x in xs:
            value = 0.0
            for region in tube.loss_region_list:
                xL = region['xL']; xR = region['xR']
                if x >= xL and x <= xR: value = region['K']/(xR-xL)
            K_over_Ls.append(value)
        #
        Ts = []
        for x in xs:
            value = tube.T_nominal
            for region in tube.T_patch_list:
                xL = region['xL']; xR = region['xR']
                if x >= xL and x <= xR: value = region['T']
            Ts.append(value)
        #
        fp.write('# n= %d\n' % tube.n)
        fp.write('# 1:x,m  2:d,m  3:area,m^2  4:K_over_L,1/m  5:Twall,K\n')
        for i in range(len(xs)):
            fp.write('%e %e %e %e %e\n' %
                     (xs[i], ds[i], math.pi*(ds[i]**2)/4, K_over_Ls[i], Ts[i]))
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

    :param x: (float) x-coordinate, in metres, of the break point
    :param d: (float) diameter, in metres, of the tube wall at the break-point.
    :returns: Number of break points defined so far.
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

    :param xL: (float) Left-end location, in metres, of the loss region.
    :param xR: (float) Right-end location, in metres, of the loss region.
    :param K: (float) Head-loss coefficient.  A value of 0.25 seems to be good for a
        reasonably smooth contraction such as the T4 main diaphragm station.
    :returns: Number of loss regions defined so far.
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

    :param xL: (float) Left-end location, in metres, of the loss region.
    :param xR: (float) Right-end location, in metres, of the loss region.
    :param T: (float) Wall temperature in degrees K.
    :returns: Number of temperature patches defined so far.
    """
    if xR < xL:
        # Keep x-values in increasing order
        xL, xR = xR, xL
    if abs(xR - xL) < 1.0e-3:
        print("Warning: temperature patch is very short: (", xL, xR, ")")
    tube.T_patch_list.append({'xL':xL, 'xR':xR, 'T':T})
    return len(tube.T_patch_list)


#----------------------------------------------------------------------
# The following classes define the objects that will be assembled into
# the gas path.
# ---------------------------------------------------------------------
  
# We will accumulate references to defined objects.
slugList = []
pistonList = []
bcList = []
diaphragmList = []
interfaceList = []
freeEndList = []
velocityEndList = []

class GasSlug(object):
    """
    Contains the gas properties and discretisation for each gas slug.

    The user may create more than one gas slug to describe the initial
    gas properties throughout the facility.
    
    Note that a slug needs to have appropriate end-conditions.
    This is achieved by creating end-condition objects such as
    FreeEnd and VelocityEnd objects and then assembling 
    the gas-path via a call to the function assemble_gas_path.
    """

    __slots__ = 'indx', 'gas', 'gmodel', 'gmodel_id', \
                'vel', 'label', 'xL', 'xR', \
                'bcL', 'bcR', 'bcL_which_end', 'bcR_which_end', \
                'nn', 'to_end_L', 'to_end_R', 'cluster_strength', \
                'viscous_effects', 'adiabatic_flag', 'hcells', \
                'ifxs'
        
    def __init__(self,
                 gmodel_id = None,
                 p = 100.0e3,
                 vel = 0.0,
                 T = 300.0,
                 T_modes = [],
                 massf = [1.0,],
                 label="",
                 nn = 10,
                 to_end_L=False,
                 to_end_R=False,
                 cluster_strength=0.0,
                 viscous_effects=0, # several options, see Lp-file documentation
                 adiabatic_flag=0,
                 hcells=[],
                 ):
        """
        Creates a gas slug with user-specified properties.

        Most parameters have default properties so that only the user
        needs to override the ones that they wish to set differently.

        Note that the locations of the ends of the slug are communicated
        through end-condition objects that are attached during assembly
        of the gas path.

        :param gmodel_id: (int) index of the gas-model file name.
        :param p: (float) Pressure in Pa.
        :param vel: (float) Velocity in m/s.
        :param T: (float) Temperature in degrees K.
        :param T_modes: (list of float) Temperatures, in K, for the other energy modes.
        :param massf: Mass fractions supplied as a list of floats 
            or a dictionary of species names and floats. 
            The number of mass fraction values should match the number 
            of species expected by the selected gas model.
        :param label: Optional (string) label for the gas slug.
        :param nn: (int) Number of cells within the gas slug.
        :param to_end_L: (bool) Flag to indicate that cells should 
            be clustered to the left end.
        :param to_end_R: (bool) Flag to indicate that cells should
            be clustered to the right end.
        :param cluster_strength: (float) As this value approaches 1.0 from above,
            the clustering gets stronger.
            A value of zero indicates no clustering.
        :param viscous_effects: (int) A nonzero value activates the viscous effects.
            0 = inviscid equations only;
            1 = include viscous source terms F_wall, loss, q,
            friction factor for pipe flow;
        :param adiabatic_flag: (int) Flag to indicate that there should
            be no heat transfer at the tube wall.
        :param hcells: Either the index (int) of a single cell or 
            a list of indices of cells for which the data are 
            to be written every dt_his seconds, as set by add_dt_plot.
            Note that cells are indexed from 0 to nn-1.
        """
        global slugList
        self.indx = len(slugList) # next available index
        # Gas data related values
        self.gmodel_id = gmodel_id
        self.gmodel = config.gmodels[gmodel_id]
        self.gas = GasState(self.gmodel)
        nsp = self.gmodel.n_species
        nmodes = self.gmodel.n_modes
        self.gas.p = p
        self.gas.massf = massf
        self.gas.T = T
        self.gas.T_modes = T_modes
        self.gas.update_thermo_from_pT()
        self.gas.update_sound_speed()
        self.gas.update_trans_coeffs()
        self.vel = vel
        self.label = label
        #
        self.nn = nn
        self.to_end_L = to_end_L
        self.to_end_R = to_end_R
        self.cluster_strength = cluster_strength
        #
        self.viscous_effects = viscous_effects
        self.adiabatic_flag = adiabatic_flag
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
        self.bcL = None
        self.bcR = None
        # The spatial limits of the gas slug are determined from
        # the boundary-condition object.
        self.xL = None
        self.xR = None
        # We also want to know which end of the other object we
        # are attached to.
        self.bcL_which_end = 'R'
        self.bcR_which_end = 'L'
        #
        slugList.append(self)
        return
    
    def write_config(self, fp):
        """
        Writes the flow state information to the specified file.
        """
        fp.write("[slug-%d]\n" % self.indx)
        fp.write("    label = %s\n" % self.label)
        fp.write("    nn = %d\n" % self.nn)
        fp.write("    cluster_to_end_L = %d\n" % self.to_end_L)
        fp.write("    cluster_to_end_R = %d\n" % self.to_end_R)
        fp.write("    cluster_strength = %e\n" % self.cluster_strength)
        fp.write("    viscous_effects = %d\n" % self.viscous_effects)
        fp.write("    adiabatic_flag = %d\n" % self.adiabatic_flag)
        #
        fp.write("    BC_L = %s\n" % boundary_control_string(self.bcL, self.bcL_which_end))
        fp.write("    BC_R = %s\n" % boundary_control_string(self.bcR, self.bcR_which_end))
        #
        hncell = len(self.hcells)
        fp.write("    hncell = %d\n" % hncell)
        fp.write("    hxcell =")
        for i in range(hncell):
            fp.write(" %d" % self.hcells[i])
        fp.write("\n")
        #
        fp.write("    initial_xL = %e\n" % self.xL)
        fp.write("    initial_xR = %e\n" % self.xR)
        fp.write("    initial_p = %e\n" % self.gas.p)
        fp.write("    initial_vel = %e\n" % self.vel)
        fp.write("    initial_T = %e\n" % self.gas.T)
        nsp = self.gmodel.n_species
        fp.write("    massf =")
        for i in range(nsp):
            fp.write(" %e" % (self.gas.massf[i]))
        fp.write("\n")
        return

    def construct_cells_and_faces(self):
        self.ifxs = roberts.distribute_points_1(self.xL, self.xR, self.nn,
                                                self.to_end_L, self.to_end_R,
                                                self.cluster_strength)
        return
    
    def write_face_data(self, fp):
        fp.write("# tindx 0\n")
        fp.write("# end\n")
        return

    def write_cell_data(self, fp):
        fp.write("# tindx 0\n")
        fp.write("# end\n")
        return

def boundary_control_string(other_object, other_object_which_end):
    """
    Assembles a boundary-condition control string for the supplied object.

    Helper function for the GasSlug class.
    """
    if isinstance(other_object, FreeEnd):
        bcs = "F  free-end"
    elif isinstance(other_object, VelocityEnd):
        bcs = "V %e  specified-velocity-end: velocity" % other_object.vel
    elif isinstance(other_object, Piston):
        bcs = "P %d  piston: piston-id" % other_object.indx
    elif isinstance(other_object, GasSlug):
        bcs = "S %d %s  slug: slug-id, slug-end-id" % \
              (other_object.indx, other_object_which_end)
    elif isinstance(other_object, Diaphragm):
        # We need to get the details of the slug attached to
        # the other side of the diaphragm.
        if other_object_which_end == 'L':
            slug_id = other_object.slugR.indx
            slug_end_id = other_object.slugR_which_end
        elif other_object_which_end == 'R':
            slug_id = other_object.slugL.indx
            slug_end_id = other_object.slugL_which_end
        else:
            raise Exception("boundary_control_string() is confused")
        bcs = "SD %d %s %d  diaphragm+slug: slug-id, slug-end-id, diaphragm-id" % \
              (slug_id, slug_end_id, other_object.indx)
    elif isinstance(other_object, Valve):
        # We need to get the details of the slug attached to
        # the other side of the Valve.
        if other_object_which_end == 'L':
            slug_id = other_object.slugR.indx
            slug_end_id = other_object.slugR_which_end
        elif other_object_which_end == 'R':
            slug_id = other_object.slugL.indx
            slug_end_id = other_object.slugL_which_end
        else:
            raise Exception("boundary_control_string() is confused")
        bcs = "SV %d %s %d valve+slug: slug-id, slug-end-id, valve-id" % \
            (slug_id, slug_end_id, other_object.indx)

    return bcs

#----------------------------------------------------------------------------

class Piston(object):
    """
    Contains the information for a piston.

    * The left- and right-end positions of the piston are
      also used to locate the ends of adjoining GasSlugs.
    * The basic piston model has inertia but no friction.
      To make accurate simulations of a particular facility,
      it is usually important to have some account of
      the friction caused by gas-seals and guide-rings that
      may be present on the piston.
    * The f_decay parameter is used to model secondary diaphragms 
      in expansion tubes as pistons which lose their mass over time.
    """

    __slots__ = 'indx', 'label', \
                'm', 'd', 'L', 'xL0', 'xR0', 'x0', 'vel0', \
                'front_seal_f', 'front_seal_area', \
                'back_seal_f', 'back_seal_area', \
                'p_restrain', 'is_restrain', 'with_brakes', 'brakes_on', \
                'x_buffer', 'hit_buffer', \
                'slugL', 'slugL_which_end', \
                'slugR', 'slugR_which_end'
    
    def __init__(self, m, d, xL0, xR0, vel0,
                 front_seal_f=0.0, front_seal_area=0.0,
                 back_seal_f=0.0, back_seal_area=0.0,
                 p_restrain=0.0, is_restrain=0,    
                 with_brakes=0, brakes_on=0,      
                 x_buffer=10.e6, hit_buffer = 0,
                 label=""):
        """
        Create a piston with specified properties.
            
        :param m: (float) Mass in kg.
        :param d: (float) Face diameter, metres.
        :param xL0: (float) Initial position of left-end, metres.
            The initial position of the piston centroid is set
            midway between xL0 and xR0 while piston length is the
            difference (xR0 - xL0).
        :param xR0: (float) Initial position of right-end, metres.
        :param vel0: (float) Initial velocity (of the centroid), m/s.
        :param front_seal_f: (float) friction coefficient.
            Typical value might be 0.2.
        :param front_seal_area: (float) Seal area over which the front-side 
            pressure acts.
            This is the effective area over which the compressed gas pressed the 
            front-side seal against the tube wall.
            Friction force is this area multiplied by downstream-pressure by
            friction coefficient.
        :param back_seal_f: (float) friction coefficient. 
            A typical value might be 0.2.
        :param back_seal_area: (float) Seal area over which the back-side 
            pressure acts.
            Friction force is this area multiplied by downstream-pressure by
            friction coefficient.  This is for gun tunnel pistons that have
            flexible skirts that are pressed onto the tube wall by the pushing gas.
        :param p_restrain: (float) Pressure at which restraint will release.
            Some machines, such as two-stage light-gas guns, will
            hold the projectile in place with some form of mechanical
            restraint until the pressure behind the piston reaches
            a critical value.  The piston is then allowed to slide.
        :param is_restrain: (int) Status flag for restraint.
            0=free-to-move, 1=restrained, 2=predefined trajectory read from external file
        :param with_brakes: (int) Flag to indicate the presence of brakes.
            0=no-brakes, 1=piston-does-have-brakes.
            Such brakes, as on the T4 shock tunnel, allow forward
            motion of the piston but prevent backward motion by
            locking the piston against the tube wall.
        :param brakes_on: (int) Flag to indicate the state of the brakes.
            0=off, 1=on.
        :param x_buffer: (float) Position of the stopping buffer in metres.
            This is the location of the piston centroid at which the piston
            would strike the buffer (or brake, in HEG terminology).
            Note that it is different to the location of the front of
            the piston at strike.
        :param hit_buffer: (int) Flag to indicate state of buffer interaction.
            A value of 0 indicates that the piston has not (yet) hit the
            buffer.
            A value of 1 indicates that it has.
            Details of the time and velocity of the strike are recorded in
            the event file.
        :param label: (string) A bit of text for corresponding line in the Lp file.
        """
        global pistonList
        self.indx = len(pistonList) # next available index
        if len(label) > 0:
            self.label = label
        else:
            # Construct a simple label.
            self.label = 'piston-' + str(self.indx)
        self.m = m
        self.d = d
        if xR0 < xL0:
            # We would like the x-values to be increasing to the right
            # but we really don't care if the piston length is zero.
            xL0, xR0 = xR0, xL0
        self.xL0 = xL0
        self.xR0 = xR0
        self.L = xR0 - xL0
        self.x0 = 0.5*(xL0 + xR0)
        self.vel0 = vel0
        self.front_seal_f = front_seal_f
        self.front_seal_area = front_seal_area
        self.back_seal_f = back_seal_f
        self.back_seal_area = back_seal_area
        self.p_restrain = p_restrain
        self.is_restrain = is_restrain
        self.with_brakes = with_brakes
        self.brakes_on = brakes_on
        self.x_buffer = x_buffer
        self.hit_buffer = hit_buffer
        #
        # The following will be assigned during assembly.
        self.slugL = None
        self.slugL_which_end = 'R'
        self.slugR = None
        self.slugR_which_end = 'L'
        #
        pistonList.append(self)
        return

    def write_config(self, fp):
        """
        Writes the piston information to the specified file.
        """
        fp.write("[piston-%d]\n" % self.indx)
        fp.write("    label = %s\n" % self.label)
        fp.write("    front_seal_f = %e\n" % self.front_seal_f)
        fp.write("    front_seal_area = %e\n" % self.front_seal_area)
        fp.write("    back_seal_f = %e\n" % self.back_seal_f)
        fp.write("    back_seal_area = %e\n" % self.back_seal_area)
        fp.write("    mass = %e\n" % self.m)
        fp.write("    diameter = %e\n" % self.d)
        fp.write("    length = %e\n" % self.L)
        fp.write("    p_restrain = %e\n" % self.p_restrain)
        fp.write("    is_restrain = %d\n" % self.is_restrain)
        fp.write("    x_buffer = %e\n" % self.x_buffer)
        fp.write("    hit_buffer = %d\n" % self.hit_buffer)
        fp.write("    with_brakes = %d\n" % self.with_brakes)
        fp.write("    brakes_on = %d\n" % self.brakes_on)
        if self.slugL != None:
            indx = self.slugL.indx
        else:
            indx = -1
        fp.write("    left-slug-id = %d\n" % indx)
        fp.write("    left-slug-end-id = %s\n" % self.slugL_which_end)
        if self.slugR != None:
            indx = self.slugR.indx
        else:
            indx = -1
        fp.write("    right-slug-id = %d\n" % indx)
        fp.write("    right-slug-end-id = %s\n" % self.slugR_which_end)
        fp.write("    x0 = %e\n" % self.x0)
        fp.write("    vel0 = %e\n" % self.vel0)
        return

    def write_data(self, fp):
        assert False, "TODO"
        return
    
#----------------------------------------------------------------------------

class Diaphragm(object):
    """
    Contains the information for a diaphragm which controls the
    interaction of two GasSlugs.
    """

    __slots__ = 'indx', 'bcindx', 'x0', 'p_burst', 'is_burst', \
                'slugL', 'slugR', \
                'slugL_which_end', 'slugR_which_end', \
                'dt_hold', \
                'dxL', 'dxR', 'label'
    
    def __init__(self,
                 x0, p_burst, is_burst=0, dt_hold=0.0,
                 dxL=0.0, dxR=0.0, label=""):
        """
        Creates a diaphragm with specified properties.

        The connections to GasSlugs are made later via the function
        assemble_gas_path.

        :param x0: (float) x-position in the tube, metres.
            This value is used to determine the end-points of the GasSlugs.
        :param p_burst: (float) Pressure, in Pa, at which rupture is triggered.
        :param is_burst: (int) Flag to indicate the state of diaphragm.
            A value of 0 indicates that the diaphragm is intact while
            a value of 1 indicates that the diaphragm is ruptured and the
            GasSlugs are interacting.
        :param dt_hold: (float) Time delay, in seconds, from rupture trigger
            to actual rupture.
        :param dxL: (float) The distance over which p is averaged on left of
            the diaphragm.  The pressure difference between the left-
            and right-sided of the diaphragm is used to trigger rupture.
            The default value of 0.0 will cause the pressure in the
            gas cell immediately adjacent to the diaphragm to be used.
        :param dxR: (float) The distance, in metres, over which p is averaged
            on right-side of the diaphragm.
        :param label: A (string) label that will appear in the parameter file
            for this diaphragm.
        """
        global diaphragmList, bcList
        self.indx = len(diaphragmList) # next available index
        self.bcindx = len(bcList)
        if len(label) > 0:
            self.label = label
        else:
            self.label = "diaphragm-" + str(self.indx)
        self.x0 = x0
        self.p_burst = p_burst
        self.is_burst = is_burst
        self.dt_hold = dt_hold
        self.dxL = dxL
        self.dxR = dxR
        #
        # The following will be assigned in assembly.
        self.slugL = None
        self.slugL_which_end = 'R'
        self.slugR = None
        self.slugR_which_end = 'L'
        #
        diaphragmList.append(self)
        bcList.append(self)
        return

    def write_config(self, fp):
        """
        Writes the diaphragm information to the specified file.
        """
        fp.write("[diaphragm-%d]\n" % self.indx)
        fp.write("    label = %s\n" % self.label)
        fp.write("    is_burst = %d\n" % self.is_burst)
        fp.write("    p_burst = %e\n" % self.p_burst)
        fp.write("    dt_hold = %e\n" % self.dt_hold)
        if self.slugL != None:
            indx = self.slugL.indx
        else:
            indx = -1
        fp.write("    left-slug-id = %d\n" % indx)
        fp.write("    left-slug-end-id = %s\n" % self.slugL_which_end)
        fp.write("    dxL = %e\n" % self.dxL)
        if self.slugR != None:
            indx = self.slugR.indx
        else:
            indx = -1
        fp.write("    right-slug-id = %d\n" % indx)
        fp.write("    right-slug-end-id = %s\n" % self.slugR_which_end)
        fp.write("    dxR = %e\n" % self.dxR)
        return

    def write_data(self, fp):
        assert False, "TODO"
        return

#-------------------------------------------------------------------------------

class GasInterface(object):
    """
    Contains the information for an interface between two slugs.

    The primary use of this class is to locate the ends of
    the connected GasSlugs.
    Implicitly, the logical connections are also made via the
    function assemble_gas_path.
    """

    __slots__ = 'bcindx', 'x0', 'slugL', 'slugL_which_end', \
                'slugR', 'slugR_which_end'
    
    def __init__(self, x0):
        """
        Creates as interface between two L{GasSlug}s at specified location.

        :param x0: (float) Initial position, in metres.
        """
        global interfaceList, bcList
        self.bcindx = len(bcList)
        self.x0 = x0
        self.slugL = None
        self.slugL_which_end = 'R'
        self.slugR = None
        self.slugR_which_end = 'L'
        interfaceList.append(self)
        bcList.append(self)
        return

    def write_config(self, fp):
        assert False, "TODO"
        return
    
#----------------------------------------------------------------------------

class FreeEnd(object):
    """
    Contains the information for a free-end condition.
    """

    __slots__ = 'bcindx', 'x0'
    
    def __init__(self, x0):
        """
        Creates a GasSlug end-condition with a specified location.

        :param x0: (float) Initial position, in metres.
        """
        global freeEndList, bcList
        self.bcindx = len(bcList)
        self.x0 = x0
        freeEndList.append(self)
        bcList.append(self)
        return

    def write_config(self, fp):
        assert False, "TODO"
        return

# --------------------------------------------------------------------

class VelocityEnd(object):
    """
    Contains the information for a fixed-velocity end condition
    for a GasSlug.
    """

    __slots__ = 'bcindx', 'x0', 'vel'
    
    def __init__(self, x0, vel=0.0):
        """
        Creates a GasSlug end-condition with a specified location
        and velocity.

        :param x0: (float) Initial position, in metres.
        :param v: (float) Velocity, in m/s, of the end-point of the GasSlug.
        """
        global velocityEndList, bcList
        self.bcindx = len(bcList)
        self.x0 = x0
        self.vel = vel
        velocityEndList.append(self)
        bcList.append(self)
        return

    def write_config(self, fp):
        assert False, "TODO"
        return

# --------------------------------------------------------------------

def assemble_gas_path(*components):
    """
    Assembles a gas path by making the logical connections between
    adjacent components.

    The components are assembled left-to-right,
    as they are supplied to this function.

    :param components: An arbitrary number of arguments representing
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
    #
    # We now need to go through the component list and,
    # for any GasInterface components, we need to connect
    # the slugs on either side.  Once this is done,
    # the GasInterface objects have done their job.
    for i in range(len(clist)):
        if isinstance(clist[i], GasInterface):
            connect_slugs(clist[i-1], clist[i+1])
    return


def connect_slugs(cL, cR):
    """
    Make the logical connection between a pair of gas slugs.
    
    :param cL: is left slug
    :param cR: is right slug

    Usually called by assemble_gas_path.
    """
    print("connect_slugs()")
    if isinstance(cL, GasSlug) and isinstance(cR, GasSlug):
        print("   Info: make slug <--> slug connection")
        cL.bcR = cR
        cL.bcR_which_end = 'L'
        cR.bcL = cL
        cR.bcL_which_end = 'R'
    else:
        raise Exception("Error: Both objects must be GasSlugs.")
    return


def connect_pair(cL, cR):
    """
    Make the logical connection between a pair of components.
    
    :param cL: is left object
    :param cR: is right object

    Usually called by assemble_gas_path.
    """
    print("connect_pair()")
    print("    left component", cL)
    print("    right component", cR)

    if isinstance(cL,VelocityEnd) and isinstance(cR, GasSlug):
        cR.bcL = cL
        cR.xL = cL.x0
        print("    velocity-end <--> gas-slug is done")
    elif isinstance(cL,FreeEnd) and isinstance(cR, GasSlug):
        cR.bcL = cL
        cR.xL = cL.x0
        print("    free-end <--> gas-slug is done")
    elif isinstance(cL,GasInterface) and isinstance(cR, GasSlug):
        cR.bcL = cL
        cR.xL = cL.x0
        print("    gas-interface <--> gas-slug is done")
    elif isinstance(cL,Piston) and isinstance(cR, GasSlug):
        cL.slugR = cR
        cL.slugR_which_end = 'L'
        cR.bcL = cL
        cR.bcL_which_end = 'R'
        cR.xL = cL.xR0
        print("    piston <--> gas-slug is done")
    elif isinstance(cL,Diaphragm) and isinstance(cR, GasSlug):
        cL.slugR = cR
        cL.slugR_which_end = 'L'
        cR.bcL = cL
        cR.bcL_which_end = 'R'
        cR.xL = cL.x0
        print("    diaphragm <--> gas-slug is done")
    elif isinstance(cL,GasSlug) and isinstance(cR, VelocityEnd):
        cL.bcR = cR
        cL.xR = cR.x0
        print("    gas-slug <--> velocity-end is done")
    elif isinstance(cL,GasSlug) and isinstance(cR, FreeEnd):
        cL.bcR = cR
        cL.xR = cR.x0
        print("    gas-slug <--> free-end is done")
    elif isinstance(cL,GasSlug) and isinstance(cR, GasInterface):
        cL.bcR = cR
        cL.xR = cR.x0
        print("    gas-slug <--> gas-interface is done")
    elif isinstance(cL,GasSlug) and isinstance(cR, Piston):
        cL.bcR = cR
        cL.bcR_which_end = 'L'
        cL.xR = cR.xL0
        cR.slugL = cL
        cR.slugL_which_end = 'R'
        print("    gas-slug <--> piston is done")
    elif isinstance(cL,GasSlug) and isinstance(cR, Diaphragm):
        cL.bcR = cR
        cL.bcR_which_end = 'L'
        cL.xR = cR.x0
        cR.slugL = cL
        cR.slugL_which_end = 'R'
        print("    gas-slug <--> diaphragm is done")
    else:
        raise Exception("    Invalid pair to connect.")
    return


# --------------------------------------------------------------------

def write_initial_files():
    """
    Writes the files needed for the main simulation code.

    These files are found in the directory config.job_name.
    """
    global config
    print("Begin write initial files.")
    if not os.path.exists(config.job_name):
        os.mkdir(config.job_name)
    #
    fp = open(config.job_name+'/config.json', 'w')
    config.write(fp)
    fp.close()
    #
    fp = open(config.job_name+'/tube.data', 'w')
    tube.write(fp)
    fp.close()
    #
    for slug in slugList:
        slug.construct_cells_and_faces()
        fileName = config.job_name+'/slug-%04d-faces.data'.format(slug.indx) 
        fp = open(fileName, 'w')
        slug.write_face_data(fp)
        fp.close()
        fileName = config.job_name+'/slug-%04d-cells.data'.format(slug.indx) 
        fp = open(fileName, 'w')
        slug.write_cell_data(fp)
        fp.close()
    #
    for piston in pistonList:
        fileName = config.job_name+'/piston-%04d.data'.format(piston.indx) 
        fp = open(fileName, 'w')
        piston.write_data(fp)
        fp.close()
    #
    for diaphragm in diaphragmList:
        fileName = config.job_name+'/diaphragm-%04d.data'.format(diaphragm.indx) 
        fp = open(fileName, 'w')
        diaphragm.write_data(fp)
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
        print("    gas slugs         :", len(slugList))
        print("    pistons           :", len(pistonList))
        print("    diaphragms        :", len(diaphragmList))
        print("    free-ends         :", len(freeEndList))
        print("    velocity-ends     :", len(velocityEndList))
        print("    gas-gas interfaces:", len(interfaceList))
        if len(slugList) < 1:
            print("Warning: no gas slugs defined; this is unusual.")
        write_initial_files()
    print("Done.")
    sys.exit(0)


