"""
Extract timing data from a set of of e4mpi simulations

@author: Nick Gibbons
"""

from numpy import array, linspace, argmax
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from subprocess import check_output, DEVNULL
from glob import glob
from re import findall, split
from platform import platform, node

import matplotlib
matplotlib.use("Agg")
plt.rcParams.update({'font.size': 14})
plt.rcParams['svg.fonttype'] = 'none'

def number(thing):
    try:
        outnumber = int(thing)
    except(ValueError):
        try:
            outnumber = float(thing)
        except(ValueError):
            outnumber = complex(thing.replace('i','j')).real
    return outnumber

def read_stdoutfile(filename):
    with open(filename) as fp:
        lines  = fp.read()

    mdata = {}
    ntasks = int(findall("MPI-parallel, number of tasks ([0-9]+)", lines)[0])
    mdata['ntasks'] = ntasks

    steps = [line for line in lines.split('\n') if line.startswith("Step= ")]
    fields = ["Step", "t", "dt", "cfl", "WC", "WCtFT", "WCtMS"]
    delims = '|'.join([i+'=' for i in fields])

    data = [split(delims, line)[1:] for line in steps]
    if len(data[-1])==2: data.pop(-1)
    data = [list(map(number,line)) for line in data]
    cols = [array(i) for i in zip(*data)]
    data = dict(zip(fields,cols))
    return mdata, data

def linear_fit(x, y):
    sx = x.sum()
    sy = y.sum()
    sx2= (x*x).sum()
    sxy= (x*y).sum()
    M = y.size
    m = (sxy - sx*sy/M)/(sx2 - sx*sx/M)
    c = (sy - sx*m)/M
    return m, c

def execution_time(T1, p, s, ncpus):
    Tn_est = s*T1 + p/ncpus*T1
    return Tn_est

def objective_function(x, ncpus, Tn):
    T1, p = x
    s = 1.0 - p
    Tn_est = execution_time(T1, p, s, ncpus)
    error = ((Tn_est - Tn)**2).sum()
    return error


platform = platform()
cpuname = check_output("lscpu | grep 'Model name:'", shell=True).decode("utf-8").split(':')[1].strip()
hostname = node()
mpi = check_output("ompi_info --version", shell=True).decode("utf-8").split('\n')[0].strip()
ldc2 = check_output("ldc2 --version", shell=True).decode("utf-8").split('\n')[0].strip()
e4dist = check_output("e4-nk-dist --help", shell=True, stderr=DEVNULL).decode("utf-8").split('\n')[1].split(':')[1].strip()

print("Begin processing NK solver parallel scaling results...")
print("        CPU: {}".format(cpuname))
print("   Platform: {}".format(platform))
print("   Hostname: {}".format(hostname))
print("   Compiler: {}".format(ldc2))
print(" Eilmer Rev: {}".format(e4dist))
print("MPI Version: {}".format(mpi))
print("\nTests:")
print("-------------------------------------------------------------------------------")


dirs = glob('test.*')
dirs.sort()
times = []
nps = []
cells = []
nspeeds = []

# First read the diagnostic files and estimate the time per iteration of each one
for dir in dirs:
    file = '{}/log.txt'.format(dir)
    with open('{}/makefile'.format(dir)) as fp:
        makefile = fp.read()
    np = int(findall('np=(\d+)', makefile)[0])

    N = 10 # Ignore the first N entries of the log file
    m, f = read_stdoutfile(file)
    d = {k:v[10:].copy() for k,v in f.items()}

    m,c = linear_fit(d['Step'], d['WC'])

    # Use grep to check how many cells in each simulation there are
    ncellscmd = "grep NELEM {}/su2grid/*.su2 | awk '{{ sum+=$2}} END {{print sum}}'".format(dir)
    ncells = check_output(ncellscmd, shell=True)
    ncells = ncells.decode("utf-8")
    ncells = int(ncells.strip())

    speed = 1.0/m
    nspeed = speed*ncells/np

    print("    {:s}: cpus={:4d} time/iter={:7.3f} (sec) nspeed={:6.1f} ncells={:d}".format(
           dir, np, m, nspeed, ncells))

    times.append(m)
    nps.append(np)
    cells.append(ncells)
    nspeeds.append(nspeed)

    #fig = plt.figure()
    #ax = plt.gca()
    #ax.plot(f['Step'], f['WC'], 'k.')
    #ax.plot(d['Step'], m*d['Step']+c, 'k-')
    #ax.set_xlabel('Step')
    #ax.set_xlabel('Time (seconds)')
    #ax.grid()
    #ax.set_title('Timing Results for directory {}'.format(dir))
    #plt.show()

print("")
testid = hash(tuple(nspeeds))
testid = hex(testid)

times = array(times)
nps = array(nps)
cells = array(cells)
nspeeds = array(nspeeds)

best = argmax(nspeeds)
bestspeed = nspeeds[best]
bestdir = dirs[best]
print("Most efficient simulation was {} with {:6.3f} cell-steps/sec/cpu".format(bestdir, bestspeed))

guess = array([times[0]*nps[0], 0.99])
f = lambda x : objective_function(x, nps, times)

# Now estimate the parallel scaling fraction
result = minimize(f, guess)
T1, p = result.x
s = 1.0-p
print("Computed parallel fraction of p={:5.5f}".format(p))
print("Test ID is {:s}".format(testid))

n = linspace(1,nps.max())
speed = T1/times
times_model = execution_time(T1, p, s, n)
speed_model = T1/times_model
speed_ideal = n 

fig = plt.figure(figsize=(8,6))
ax = plt.gca()
ax.plot(nps, speed, 'ko')
ax.plot(n, speed_model, 'k-')
ax.plot(n, speed_ideal, 'k--')
ax.set_xlabel('nps')
ax.set_xlabel('Speedup')
ax.set_title("e4mpi Solver Scaling: {}".format(testid))
ax.grid()
plt.savefig("results.svg")
plt.close()
#plt.show()

