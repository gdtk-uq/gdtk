#! /usr/bin/env python3
# optimize.py
# Automate the running of the flow simulations while searching
# for the optimum parameters (angles) that define the nozzle shape.
#
# PJ, 2018-03-04, take bits from nenzfr

import sys, os
DGDINST = os.path.expandvars("$HOME/dgdinst")
sys.path.append(DGDINST)

import shlex, subprocess, string, math
from scipy.optimize import minimize

def run_command(cmdText):
    """
    Run the command as a subprocess.
    """
    # Flush before using subprocess to ensure output
    # text is in the right order.
    sys.stdout.flush()    
    if (type(cmdText) is list):
        args = cmdText
    else:
        args = shlex.split(cmdText)
    print("About to run cmd:", " ".join(args))
    return subprocess.check_call(args)

def prepare_input_script(substituteDict):
    """
    Prepare the actual input file for Eilmer4 from a template
    which has most of the Lua input script in place and just
    a few place-holders that need to be substituted for actual
    values.
    """
    fp = open("nozzle.template.lua", 'r')
    text = fp.read()
    fp.close()
    template = string.Template(text)
    text = template.substitute(substituteDict)
    fp = open("nozzle.lua", 'w')
    fp.write(text)
    fp.close()
    return

def run_simulation(param_dict):
    run_command('prep-gas ideal-air.inp ideal-air-gas-model.lua')
    prepare_input_script(param_dict)
    run_command('e4shared --prep --job=nozzle')
    run_command('e4shared --run --job=nozzle --verbosity=1')
    return

def post_simulation_files():
    run_command('e4shared --post --job=nozzle --tindx-plot=all'+
                ' --vtk-xml --add-vars="mach,pitot,total-p,total-h"')
    run_command('e4shared --post --job=nozzle --slice-list="1,0,:,0"'+
                ' --output-file="nozzle-throat.data"')
    run_command('e4shared --post --job=nozzle --slice-list="1,$,:,0"'+
                ' --output-file="nozzle-exit.data"')
    return

def neg_thrust(tindx):
    """
    Read the loads file and compute the x-component of thrust.

    Input:
    tindx : integer specifying which loads file
    """
    fileName = 'loads/t%04d-loads.dat' % tindx
    print("Estimating thrust from loads file ", fileName)
    f = open(fileName, 'r')
    thrustx = 0.0
    for line in f.readlines():
        items = line.strip().split()
        if items[0] == '#': continue
        if len(items) < 11: break
        # r = float(items[1]) # 2:pos.y
        dA = float(items[3]) # 4:area per radian for axisymmetric simulation
        nx = float(items[10]) # 11:n.x
        p = float(items[9]) # 10:sigma
        thrustx = thrustx + 2*math.pi*dA*nx*p
    return thrustx
    
def objective(params, *args):
    """
    Given a list of parameter values, run a simulation and compute the thrust.
    Since the thrust is in the negative x-direction, large negative values
    are good.  So, the minimizer will drive toward good values.
    """
    pdict = {"theta_init":params[0], "alpha":params[1],
             "beta":params[2], "theta_cone":params[3]}
    run_simulation(pdict)
    # Note that, if we run several simulations, there will be several
    # loads files, indexed in order of creation.
    # We want the most recent, if we are in the optimization process.
    f = open('loads/nozzle-loads.times', 'r')
    tindx = 0
    for line in f.readlines():
        items = line.strip().split()
        if items[0] == '#': continue
        if len(items) < 2: break
        tindx = int(items[0])
    thrustx = neg_thrust(tindx)
    print("pdict=", pdict, "thrust=", thrustx)
    return thrustx

def main():
    """
    This script was built in stages.
    The if-statements are for testing the functions as the script
    was being developed. They might be still useful for exploring.
    """
    if 0:
        print("Let's run a simulation.")
        pdict = {"theta_init":30.0, "alpha":0.0, "beta":0.0, "theta_cone":30.0}
        run_simulation(pdict)
    if 0:
        print("Compute thrust from previously run simulation.")
        print("thrust=", neg_thrust(0))
    if 0:
        print("Evaluate objective function.")
        params = [30.0, 0.0, 0.0, 30.0] # [theta_init, alpha, beta, theta_cone]
        objv = objective(params)
        print("objective value=", objv)
    if 1:
        print("Let the optimizer take control and run the numerical experiment.")
        x0 = [30.0, 0.0, 0.0, 30.0] # [theta_init, alpha, beta, theta_cone]
        result = minimize(objective, x0, method='Nelder-Mead',
                          options={'disp':True, 'maxfev':50})
        print('optimized result:')
        print('    x=', result.x)
        print('    fun=', result.fun)
        print('    success=', result.success)
        print('    message=', result.message)
    return

if __name__ == "__main__":
    main()
