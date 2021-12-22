#! /usr/bin/env python3
# optimize.py
# Automate the running of the flow simulations while searching
# for the optimum parameters (angles) that define the nozzle shape.
#
# To monitor the progress of the optimizer, you can run the command:
# tail -f progress.txt
# to see the result of each objective evaluation.
#
# PJ, 2018-03-04, take bits from nenzfr
#     2021-12-22, update to accommodate loads-file changes, use nelmin
#                 and run jobs in their own directories

import sys, os
DGDINST = os.path.expandvars("$HOME/dgdinst")
sys.path.append(DGDINST)

import shutil, shlex, subprocess, string, math
from eilmer.nelmin import minimize, NelderMeadMinimizer
import time

start_time = time.time()
progress_file = open("progress.txt", 'w')
progress_file.write("# job wall_clock params[0] params[1] params[2] params[3] thrustx\n")
progress_file.flush()
obj_eval_number = 0

def run_command(cmdText, jobDir, logFile=None):
    """
    Run the command as a subprocess in directory jobDir.
    If a logFile name is provided, capture stdout+stderr and write to that file,
    else just let those streams go to the console.
    """
    # Flush before using subprocess to ensure output
    # text is in the right order.
    sys.stdout.flush()
    if (type(cmdText) is list):
        args = cmdText
    else:
        args = shlex.split(cmdText)
    # print("jobDir:", jobDir, "cmdText:", " ".join(args)) # debug
    captureFlag = True if logFile else False
    result = subprocess.run(args, cwd=jobDir, check=True, capture_output=captureFlag)
    if captureFlag:
        f = open(logFile, 'w')
        f.write(result.stdout.decode())
        f.write(result.stderr.decode())
        f.close()
    return result.returncode

def prepare_input_script(substituteDict, jobDir):
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
    fp = open(jobDir+"/nozzle.lua", 'w')
    fp.write(text)
    fp.close()
    return

def run_simulation(param_dict, jobDir):
    """
    Prepare and run a simulation in it own directory.
    We do this so that several simulations may run concurrently
    and so that we can easily find the data files later.
    """
    if not os.path.exists(jobDir): os.mkdir(jobDir)
    shutil.copy('ideal-air.inp', jobDir)
    logFile = jobDir+'.log'
    run_command('prep-gas ideal-air.inp ideal-air-gas-model.lua', jobDir, logFile)
    prepare_input_script(param_dict, jobDir)
    run_command('e4shared --prep --job=nozzle', jobDir, logFile)
    run_command('e4shared --run --job=nozzle --max-cpus=4 --verbosity=1', jobDir, logFile)
    return

def post_simulation_files(jobDir):
    """
    Postprocess a simulation is its own directory.
    Although not used by the optimizer, this function may be handy
    when exploring the behaviour of the optimization procedure.
    """
    run_command('e4shared --post --job=nozzle --tindx-plot=all'+
                ' --vtk-xml --add-vars="mach,pitot,total-p,total-h"', jobDir)
    run_command('e4shared --post --job=nozzle --slice-list="1,0,:,0"'+
                ' --output-file="nozzle-throat.data"', jobDir)
    run_command('e4shared --post --job=nozzle --slice-list="1,$,:,0"'+
                ' --output-file="nozzle-exit.data"', jobDir)
    return

def neg_thrust(tindx, jobDir):
    """
    Read the loads file and return the x-component of thrust.

    Input:
    tindx : integer specifying which loads file to inspect
    jobDir : the dedicated directory for the simulation files.
    """
    fileName = jobDir+'/loads/t%04d/b0001.t%04d.loads.dat' % (tindx, tindx)
    print("Estimating thrust from loads file ", fileName) # debug
    f = open(fileName, 'r')
    thrustx = 0.0
    for line in f.readlines():
        items = line.strip().split()
        if items[0] == '#': continue
        if len(items) < 10: break
        dA = float(items[6]) # 7:area per radian for axisymmetric simulation
        nx = float(items[3]) # 4:n.x
        p = float(items[9]) # 10:p
        thrustx = thrustx + 2*math.pi*dA*nx*p
    print("thrustx=", thrustx) # debug
    return thrustx

def objective(params, *args):
    """
    Given a list of parameter values, run a simulation and compute the thrust.
    Since the thrust is in the negative x-direction, large negative values are good.
    The minimizer will drive toward good values.

    [TODO] make this function thread-safe.
    Presently it is reliable for n_workers=1, only.
    """
    global start_time, progress_file, obj_eval_number
    pdict = {"theta_init":params[0], "alpha":params[1],
             "beta":params[2], "theta_cone":params[3]}
    obj_eval_number += 1
    print("Start job number:", obj_eval_number)
    jobDir = 'job-%04d' % obj_eval_number
    run_simulation(pdict, jobDir)
    # Note that, if we run several simulations, there will be several
    # loads files, indexed in order of creation.
    # We want the most recent, if we are in the optimization process.
    f = open(jobDir+'/loads/nozzle-loads.times', 'r')
    tindx = 0
    for line in f.readlines():
        items = line.strip().split()
        if items[0] == '#': continue
        if len(items) < 2: break
        tindx = int(items[0])
    thrustx = neg_thrust(tindx, jobDir)
    progress_file.write("%d %.1f %.4f %.4f %.4f %.4f %.2f\n" %
                        (obj_eval_number, time.time()-start_time, params[0],
                         params[1], params[2], params[3], thrustx))
    progress_file.flush() # so that we can see the results as simulations run
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
        run_simulation(pdict, "job-0000")
    if 0:
        print("Compute thrust from previously run simulation.")
        print("thrust=", neg_thrust(0, 'job-0000'))
    if 0:
        print("Evaluate objective function.")
        params = [30.0, 0.0, 0.0, 30.0] # [theta_init, alpha, beta, theta_cone]
        obj_eval_number = 1
        objv = objective(params)
        print("objective value=", objv)
    if 1:
        print("Let the optimizer take control and run the numerical experiment.")
        x0 = [30.0, 0.0, 0.0, 30.0] # [theta_init, alpha, beta, theta_cone] degrees
        result = minimize(objective, x0, [2.0, 2.0, 2.0, 2.0],
                          options={'tol':1.0e-4, 'P':2, 'maxfe':60, 'n_workers':1})
        print('optimized result:')
        print('  x=', result.x)
        print('  fx=', result.fun)
        print('  convergence-flag=', result.success)
        print('  number-of-fn-evaluations=', result.nfe)
        print('  number-of-restarts=', result.nrestarts)
        print('  vertices=', [str(v) for v in result.vertices])
    #
    print("overall calculation time:", time.time()-start_time)
    return

# Let's actually do some work...
main()
progress_file.close()
