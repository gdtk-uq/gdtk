#! /usr/bin/env python3
# optimize.py
# Automate the running of the flow simulations while searching for
# the optimum parameters that define the truncated Busemann inlet.
#
# PJ, 2020-09-19, adapt from the 2D/scramjet example, single worker thread.

import sys, os
DGDINST = os.path.expandvars("$HOME/dgdinst")
sys.path.append(DGDINST)

import shlex, subprocess, string, math
from eilmer.nelmin import minimize
import numpy
import logging
logging.basicConfig(
    filename="optimize.log.txt",
    format="%(levelname)-10s %(created)f %(threadName)s %(message)s",
    level=logging.NOTSET
)

def run_command(cmdText, fileName=None):
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
    logging.info("About to run cmd: "+" ".join(args))
    fo = open(fileName, 'w') if fileName else None
    return subprocess.check_call(args, stdout=fo)

def prepare_input_script(substituteDict):
    """
    Prepare the actual input file for Eilmer4 from a template
    which has most of the Lua input script in place and just
    a few place-holders that need to be substituted for actual
    values.
    """
    fp = open("tbo.template.lua", 'r')
    text = fp.read()
    fp.close()
    template = string.Template(text)
    text = template.substitute(substituteDict)
    fp = open("tbo.lua", 'w')
    fp.write(text)
    fp.close()
    return

def run_simulation(param_dict):
    # run_command('prep-gas ideal-air.inp ideal-air-gas-model.lua')
    prepare_input_script(param_dict)
    run_command('e4shared --prep --job=tbo', 'prep.log.txt')
    run_command('mpirun -np 6 --oversubscribe -- e4mpi --run --job=tbo',
                'run.log.txt')
    run_command('e4shared --post --job=tbo --tindx-plot=last'+
                ' --add-vars="mach,total-p"'+
                ' --slice-list="5,$,:,0"'+
                ' --output-file=profile.data',
                'post.log.txt')
    return

def objective(params, *args):
    """
    Given a list of parameter values, run a simulation and compute
    a measure of the quality of the compressed-flow after the inlet.
    Remember that the minimizer will try to reduce the magnitude
    of this measure.
    """
    pdict = {
        "Y0":params[0], "Y1":params[1], "Y2":params[2],
        "Y3":params[3], "Y4":params[4], "Y5":params[5],
        "Y6":params[6], "Y7":params[7], "Y8":params[8]
    }
    try:
        run_simulation(pdict)
        #
        p_inf = 1000.0
        p_total_inf = 4613.0 * p_inf # from NACA 1135 table for M=7.12
        PR_target = 340.1 # from Busemann ideal solution
        #
        data = numpy.loadtxt("profile.data", skiprows=1)
        ys = data[:,1] # remember that Python array indices start at 0
        n = len(ys)
        p_static = sum(data[:,8])/n
        p_total = sum(data[:,19])/n
        mach = sum(data[:,18])/n
        logging.info(f"p_static={p_static}, PR={p_static/p_inf}")
        logging.info(f"p_total={p_total}, normalized={p_total/p_total_inf}")
        logging.info(f"mach={mach}")
        #
        low_M_penalty = 10.0 if mach < 1.5 else 0.0
        #
        obj = 3.0*abs(p_static/p_inf/PR_target - 1.0) + \
            (1.0 - p_total/p_total_inf) + low_M_penalty
    except:
        obj = 10000.0 # a big penalty for a failed simulation
    logging.info(f"pdict={pdict}, obj={obj}")
    return obj

def main():
    """
    This script was built in stages.
    The if-statements are for testing the functions as the script
    was being developed. They might be still useful for exploring.
    """
    if 0:
        print("Let's run a simulation.")
        # Initial truncated Busemann points.
        pdict = {"Y0":0.881750, "Y1":0.826332, "Y2":0.823080,
                 "Y3":0.674828, "Y4":0.744892, "Y5":0.508068,
                 "Y6":0.505361, "Y7":0.287438, "Y8":0.146429}
        run_simulation(pdict)
    if 0:
        print("Evaluate objective function.")
        # Initial truncated Busemann points.
        # params = [0.881750, 0.826332, 0.823080,
        #           0.674828, 0.744892, 0.508068,
        #           0.505361, 0.287438, 0.146429]
        # After second run of optimizer with 105+108 fn evaluations.
        # fun= 0.3966710038268789
        params = [0.8951219647936348, 0.8458731967690415, 0.8096059871085017,
                  0.6492786657904442, 0.7438374618335349, 0.471683442280917,
                  0.5128294314724264, 0.31951280691860723, 0.1457829739923719]
        objv = objective(params)
        print("objective value=", objv)
    if 1:
        print("Let the optimizer take control and run the numerical experiment.")
        # Initial truncated Busemann points.
        # x0 = [0.881750, 0.826332, 0.823080,
        #       0.674828, 0.744892, 0.508068,
        #       0.505361, 0.287438, 0.146429]
        # After frist run of optimizer with 105 fn evaluations.
        # fun= 0.44387911743189534
        x0 = [0.8997840012046376, 0.844200184878106, 0.8299524988526994,
              0.6531033284047321, 0.7494029828177131, 0.4819390874141007,
              0.5132636478978485, 0.29685473312349775, 0.1460182292909658]
        dx = [0.02, 0.02, 0.02,
              0.02, 0.02, 0.02,
              0.02, 0.02, 0.02]
        result = minimize(objective, x0, dx, options={'P':3, 'maxfe':50})
        print('optimized result:')
        print('  x=', result.x)
        print('  fun=', result.fun)
        print('  success=', result.success)
        print('  nfe=', result.nfe)
        print('  nrestarts=', result.nrestarts)
    return

if __name__ == "__main__":
    main()
