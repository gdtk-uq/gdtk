#! /usr/bin/env python3
# optimize.py
# Automatic code for finding the optimum nozzle contour for a Mach 4 nozzle in the T4 shock tunnel.
#
# originally written by PJ 2018-03-04 for ENGG7601 assignment
#
# Objective function used here is that of Craddock (2000), and Chan (2013)
# which aims to reduce the flow non-uniformity at the exit plane of the nozzle.
# In particular, the functions for flow uniformity and the objective function
# were adapted from Wilson's design files.
# cfcfd3/examples/eilmer3/2D/nenzfr-optimisation-UQ-X3R-M7-nozzle/optimise_X3R_Mach_7_nozzle.py
#
# Authors:
# Peter J., Daniel Smith, Wilson Chan
#
# 25-June-2020

import sys, os
import numpy
import math
from scipy.optimize import minimize
import string, shlex, subprocess

# Global parameters for the design.
jobname = "t4m4b"
# Design targets
M_target = 4.0        # Mach number
dM_target = 0.01      # Variation in Mach number
# Desired flow angle is zero, of course.
dtheta_target = 0.02  # Variation in outflow angle (in degrees)
# Weighting parameters.
phi_theta = 1.0 / math.tan(math.radians(dtheta_target))
phi_M = 1.0 / dM_target

numpy.set_printoptions(linewidth=100)

def run_command(cmd):
    """
    Run the command as a subprocess, saving the standard streams to text files.
    Input is a tuple, as indicated below.
    """
    cmdText, tagText = cmd
    if (type(cmdText) is list):
        args = cmdText
    else:
        args = shlex.split(cmdText)
    if not os.path.exists("logs"): os.mkdir("logs")
    with open(f"logs/stdout-{tagText}.txt","w") as out, \
         open(f"logs/stderr-{tagText}.txt","w") as err:
        completedCmd = subprocess.run(args, stdout=out, stderr=err)
    return completedCmd.returncode

def run_simulation(param_dict):
    """
    Run a flow calculation with particular shape parameters
    provided by the optimizer.
    """
    fp = open("nozzle.template.lua", 'r')
    text = fp.read()
    fp.close()
    template = string.Template(text)
    text = template.substitute(param_dict)
    fp = open(f"{jobname}.lua", 'w')
    fp.write(text)
    fp.close()
    return_code = run_command((f'e4shared-debug --prep --job={jobname}', 'prep'))
    if return_code != 0:
        raise RuntimeError("Preparation for simulation failed.")
    return_code = run_command((f'mpirun -np 2 e4mpi --run --job={jobname} --verbosity=1', 'run'))
    if return_code != 0:
        raise RuntimeError("MPI run for simulation failed.")
    return_code = run_command(
        (f'e4shared --post --job={jobname} --tindx-plot=last --add-vars="mach,pitot"' +
         f' --slice-list="60:62,$,:,0" --output-file="{jobname}-exit.data"',
         'post')
    )
    if return_code != 0:
        raise RuntimeError("Postprocessing for simulation failed.")
    return

def flow_uniformity():
    """
    Reads in the most recent processed data called jobname-exit.data and
    returns measures of the badness of the Mach number and flow angle profiles.
    """
    # Load in the most recent results
    data = numpy.loadtxt(f"{jobname}-exit.data", skiprows=1)
    # Secondary functions that contribute to the objective function.
    f_theta = 0.0; f_M = 0.0 # Initialise both functions to zero first.
    N = 0 # Initialise the counter for the number of cells in the core flow.
    for i in range(len(data[:,0])):
        # Definition used by Chris Craddock to estimate the boundary layer edge.
        M_i = data[i,18] # mach number
        y_i = data[i,1] # y-position
        vely_i = data[i,6]
        velx_i = data[i,5]
        if i == 0:
            dMdy = 0.0  # Set to some number so that the first point
                        # is not set as the boundary layer edge.
        else:
            M_im1 = data[i-1,18]
            y_im1 = data[i-1,1]
            dMdy = (M_i - M_im1) / (y_i - y_im1)
        # If dMdy >= -20.0, then we are in the core flow.
        if dMdy >= -20.0:
            f_theta += (vely_i / velx_i)**2
            f_M += (M_i - M_target)**2
            N += 1
    #
    # Weight the secondary functions by weighting parameters.
    f_theta = phi_theta**2 / N * f_theta
    f_M = phi_M**2 / N * f_M
    return f_theta, f_M

def objective(params):
    """
    Given a list of parameter values, run a simulation and compute the flow uniformity
    """
    pdict = {"dy1":params[0], "dy2":params[1], "dy3":params[2], "dy4":params[3],
             "dy5":params[4], "dy6":params[5], "dy7":params[6]}
    print(40*"-")
    print("Commencing simulation for x = ", params)
    run_simulation(pdict)
    print("Simulation complete.")
    f_theta, f_M = flow_uniformity()
    obj_funct = (f_theta + f_M)**2
    print("Objective function evaluated to ", obj_funct)
    return obj_funct


def main():
    """
    This script was built in stages.
    The if-statements are for testing the functions as the script
    was being developed. They might be still useful for exploring.
    """
    if 0:
        print("Let's run a simulation.")
        pdict = {"dy1":0.0, "dy2":0.0, "dy3":0.0, "dy4":0.0, "dy5":0.0, "dy6":0.0, "dy7":0.0}
        run_simulation(pdict)
    if 0:
        print("Compute the flow uniformity from previously run simulation.")
        f_theta, f_M = flow_uniformity()
        print("Mach uniformity = ", f_M)
        print("angular uniformity = ", f_theta)
    if 0:
        print("Evaluate objective function.")
        params = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        objv = objective(params)
        print("objective value=", objv)
    if 1:
        print("Let the optimizer take control and run the numerical experiment.")
        x0 = numpy.array([[-1.01359136e-03, 1.40857503e-05, -9.27802658e-04, 3.82603999e-04,
                           2.39843018e-03, -2.91766389e-04, -3.60642387e-04]])
        result = minimize(objective, x0, method='Nelder-Mead',
                          options={'disp':True, 'maxiter':5000})
        print('optimized result:')
        print('  x=', result.x)
        print('  fun=', result.fun)
        print('  success=', result.success)
        print('  message=', result.message)

        # Save the new set of bezier points for later.
        bezCtrlPts_orig = numpy.loadtxt(f"Bezier-control-pts-{jobname}-initial.data", skiprows=1)
        bezCtrlPts_opt = bezCtrlPts_orig
        for i in range(len(result.x)):
            bezCtrlPts_opt[i+2,1] = bezCtrlPts_opt[i+2,1] + result.x[i]
        numpy.savetxt("Bezier-control-pts-{jobname}.opt.data", bezCtrlPts_opt,
                      header='#     x, m       y, m')
    return

if __name__ == "__main__":
    main()
