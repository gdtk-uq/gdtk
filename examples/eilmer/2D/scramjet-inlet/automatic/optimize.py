#! /usr/bin/env python3
# optimize.py
# Automate the running of the flow simulations while searching
# for the optimum parameters that define the scramjet inlet.
#
# PJ, 2020-08-03, adapt from the nozzle-optimize-bell example.

import sys, os
DGDINST = os.path.expandvars("$HOME/dgdinst")
sys.path.append(DGDINST)

import shlex, subprocess, string, math
from scipy.optimize import minimize
import numpy

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
    fp = open("inlet.template.lua", 'r')
    text = fp.read()
    fp.close()
    template = string.Template(text)
    text = template.substitute(substituteDict)
    fp = open("inlet.lua", 'w')
    fp.write(text)
    fp.close()
    return

def run_simulation(param_dict):
    run_command('prep-gas ideal-air.inp ideal-air-gas-model.lua')
    prepare_input_script(param_dict)
    run_command('e4shared --prep --job=inlet')
    run_command('e4shared --run --job=inlet --verbosity=1')
    run_command('e4shared --post --job=inlet --tindx-plot=last'+
                ' --add-vars="mach,total-p"'+
                ' --slice-list="1,$,:,0"'+
                ' --output-file=profile.data')
    return

def objective(params, *args):
    """
    Given a list of parameter values, run a simulation and compute
    a measure of the quality of the compressed-flow after the inlet.
    Remember that the minimizer will try to reduce the magnitude
    of this measure.
    """
    pdict = {"THETA":params[0], "H0":params[1], "H1":params[2], "H2":params[3]}
    try:
        run_simulation(pdict)
        #
        p_inf = 10000.0
        p_total_inf = 152.0 * p_inf # from NACA 1135 table for M=4
        PR_target = 12.0 # from Mike Smart's 1999 paper
        #
        data = numpy.loadtxt("profile.data", skiprows=1)
        ys = data[:,1] # remember that Python array indices start at 0
        n = len(ys)
        p_static = sum(data[:,8])/n
        p_total = sum(data[:,19])/n
        mach = sum(data[:,18])/n
        print("p_static=", p_static, "PR=", p_static/p_inf)
        print("p_total=", p_total, "normalized=", p_total/p_total_inf)
        print("mach=", mach)
        #
        low_M_penalty = 100.0 if mach < 1.5 else 0.0
        #
        obj = abs(p_static/p_inf - PR_target) + (1.0 - p_total/p_total_inf) + low_M_penalty
    except:
        obj = 10000.0 # a big penalty for a failed simulation
    print("pdict=", pdict, "obj=", obj)
    return obj

def main():
    """
    This script was built in stages.
    The if-statements are for testing the functions as the script
    was being developed. They might be still useful for exploring.
    """
    if 0:
        print("Let's run a simulation.")
        pdict = {"THETA":10.0, "H0":0.05, "H1":0.05, "H2":0.15}
        run_simulation(pdict)
    if 0:
        print("Evaluate objective function.")
        params = [10.0, 0.05, 0.05, 0.15] # [THETA, H0, H1, H2]
        objv = objective(params)
        print("objective value=", objv)
    if 1:
        print("Let the optimizer take control and run the numerical experiment.")
        x0 = [10.0, 0.05, 0.05, 0.15] # [THETA, H0, H1, H2]
        result = minimize(objective, x0, method='Nelder-Mead',
                          options={'disp':True, 'maxfev':50})
        print('optimized result:')
        print('  x=', result.x)
        print('  fun=', result.fun)
        print('  success=', result.success)
        print('  message=', result.message)
    return

if __name__ == "__main__":
    main()
