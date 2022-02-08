#! /usr/bin/env python3
# optimize.py
# Automate the running of the flow simulations while searching
# for the optimum parameters (angles and displacements) that define the nozzle shape.
#
# To monitor the progress of the optimizer, you can run the command:
# tail -f progress.txt
# to see the result of each objective evaluation.
#
# PJ, 2018-03-04, Eilmer4 implementation, take bits from nenzfr
#     2021-12-22, update to accommodate loads-file changes, use nelmin
#                 and run jobs in their own directories
#     2022-02-09, Puffin implmentation for CfH seminar.

import sys, os
DGDINST = os.path.expandvars("$HOME/dgdinst")
sys.path.append(DGDINST)

import shutil, shlex, subprocess, queue
import string, math
import time
from eilmer.nelmin import minimize, NelderMeadMinimizer

start_time = time.time()
progress_file = open("progress.txt", 'w')
progress_file.write("# job wall_clock theta_init theta_cone deltas_0 deltas_1 deltas_2 thrustx\n")
progress_file.flush()
# Each of the Puffin simulations will be run in its own directory, identified by job number,
# so that it can be run independently of all other simulations.
# We will use a queue to regulate the specification of this job number.
obj_eval_queue = queue.Queue()
obj_eval_queue.put(0)

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
    fp = open("nozzle.template.py", 'r')
    text = fp.read()
    fp.close()
    template = string.Template(text)
    text = template.substitute(substituteDict)
    fp = open(jobDir+"/nozzle.py", 'w')
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
    run_command('puffin-prep --job=nozzle', jobDir, logFile)
    run_command('puffin --job=nozzle --verbosity=1', jobDir, logFile)
    return

def neg_thrust(jobDir):
    """
    Read the loads file and return the x-component of thrust.

    Input:
    jobDir : the dedicated directory for the simulation files.
    """
    run_command('puffin-post --job=nozzle --output=stream --cell-index=$ --stream-index=0', jobDir)
    fileName = jobDir+'/nozzle/nozzle-streamline-0-cell-39.data'
    print("Estimating thrust from streamline file ", fileName) # debug
    f = open(fileName, 'r')
    ys = []; ps = []
    for line in f.readlines():
        items = line.strip().split()
        if items[0] == '#': continue
        if len(items) < 10: break
        ys.append(float(items[1])) # 2:y
        ps.append(float(items[6])) # 7:p
    thrustx = 0.0
    for i in range(len(ps)-1):
        y_avg = 0.5*(ys[i+1]+ys[i])
        dy = ys[i+1]-ys[i]
        p_avg = 0.5*(ps[i+1]+ps[i])
        dA = 2*math.pi*y_avg*dy
        thrustx = thrustx - dA*p_avg
    print("thrustx=", thrustx) # debug
    return thrustx

def objective(params, *args):
    """
    Given a list of parameter values, run a simulation and compute the thrust.
    Since the thrust is in the negative x-direction, large negative values are good.
    The minimizer will drive toward good values.
    """
    global start_time # Immutable.
    global progress_file # Shared output stream, hopefully not too messed up.
    global obj_eval_queue # Thread-safe.
    job_number = obj_eval_queue.get()
    job_number += 1
    obj_eval_queue.put(job_number)
    print("Start job number:", job_number)
    pdict = {"theta_init":params[0], "theta_cone":params[1],
             "deltas_0":params[2], "deltas_1":params[3], "deltas_2":params[4]}
    jobDir = 'job-%04d' % job_number
    run_simulation(pdict, jobDir)
    thrustx = neg_thrust(jobDir)
    progress_file.write("%d %.1f %.4f %.4f %.4f %.4f %.4f %.2f\n" %
                        (job_number, time.time()-start_time, params[0],
                         params[1], params[2], params[3], params[4], thrustx))
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
        pdict = {"theta_init":30.0, "deltas_0":0.0, "deltas_1":0.0, "deltas_2":0.0, "theta_cone":30.0}
        run_simulation(pdict, "job-0000")
    if 0:
        print("Compute thrust from previously run simulation.")
        print("thrust=", neg_thrust('job-0000'))
    if 0:
        print("Evaluate objective function.")
        params = [30.0, 30.0, 0.0, 0.0, 0.0] # [theta_init, theta_cone, deltas_0, _1, _2]
        obj_eval_number = 1
        objv = objective(params)
        print("objective value=", objv)
    if 1:
        print("Let the optimizer take control and run the numerical experiment.")
        x0 = [30.0, 30.0, 0.0, 0.0, 0.0] # [theta_init, theta_cone, deltas_0, _1, _2] angles in degrees
        result = minimize(objective, x0, [2.0, 2.0, 0.01, 0.01, 0.01],
                          options={'tol':1.0e-4, 'P':2, 'maxfe':60, 'n_workers':2})
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
