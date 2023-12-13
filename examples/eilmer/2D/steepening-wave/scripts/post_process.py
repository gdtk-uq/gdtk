"""
Check the order of accuracy using the steepening wave problem.

@author: Nick Gibbons
"""

from numpy import sqrt, sum, absolute, mean, diff, log10, linspace
#from pylab import plot,show,loglog,legend,xlabel,ylabel,grid
from collections import defaultdict
from libpost import *
from analytic_solution import *
from sys import argv

def root_mean_square(x):
    return sqrt(sum(x*x)/x.size)

    
def linear_fit(x, y):
    sx = x.sum()
    sy = y.sum()
    sx2= (x*x).sum()
    sxy= (x*y).sum()
    M = y.size
    m = (sxy - sx*sy/M)/(sx2 - sx*sx/M)
    c = (sy - sx*m)/M
    return m, c

dirs = sorted(argv[1:])
if len(dirs)==0: raise Exception("No simulations given")

cfds = []
sols = []

# Loop over each simulation directory and compute the analytic solution
# at its specific set of x locations.
for dir in dirs:
    f = get_fluid_block_array(dir)
    cfd = {k:v[1,:].copy() for k,v in f.items()}

    x = cfd['pos.x']

    M_inf = -1.0
    p_inf = 1e5
    rho_inf = 1.0
    wave = SteepeningWave(x, rho_inf, p_inf, M_inf)
    t = 0.5*wave.shock_formation_time()
    p, rho, u = wave.gen_profiles(t)

    sol = {'pos.x':x, 'p':p, 'rho': rho, 'vel.x':u}

    sols.append(sol)
    cfds.append(cfd)

# Loop over the variables we care about and compile some error
# estimates for each one, into a list that is n simulations long
variables = ['p', 'rho', 'vel.x']
L0s = defaultdict(list)
L2s = defaultdict(list)

for sol,cfd in zip(sols, cfds):
    for var in variables:
        s = sol[var]
        c = cfd[var]
        error = s - c
        L0s[var].append(absolute(error).max())
        L2s[var].append(root_mean_square(error))
        #if var=='rho':
        #    plot(cfd['pos.x'], sol['rho'], 'k-')
        #    plot(sol['pos.x'], cfd['rho'], 'ro')
        #    show()

# We want these lists to be arrays later, so convert them here
L0s = {k:array(v) for k,v in L0s.items()}
L2s = {k:array(v) for k,v in L2s.items()}
dxs= array([mean(diff(cfd['pos.x'])) for cfd in cfds])

# Extract the order of convergence by fitting a slope in log/log space.
print("Order of Convergence: ")
for var in variables:
    m,c = linear_fit(log10(dxs),log10(L0s[var]))
    L0s[var+'m'] = m
    L0s[var+'c'] = c
    print(" {:>5s} L0={:3.2f}".format(var, m))

    m,c = linear_fit(log10(dxs),log10(L2s[var]))
    L2s[var+'m'] = m
    L2s[var+'c'] = c
    print(" {:>5s} L2={:3.2f}".format(var, m))

# Plot the results on the same graph
#log10dxs = log10(dxs)
#colours = ['red', 'blue', 'green']
#for i,var in enumerate(variables):
#    colour = colours[i]
#    plot(log10(dxs),log10(L0s[var]), linestyle="None", marker='o', color=colour, label=var)
#    plot(log10(dxs),log10(L2s[var]), linestyle="None", marker='+', color=colour)
#
#    fit_log10dxs = linspace(int(log10dxs.min())-1, int(log10dxs.max()))
#    plot(fit_log10dxs, fit_log10dxs*L0s[var+'m'] + L0s[var+'c'],  linestyle='--', color=colour)
#    plot(fit_log10dxs, fit_log10dxs*L2s[var+'m'] + L2s[var+'c'],  linestyle='-', color=colour)
#
#grid()
#xlabel('log10(dxs)')
#ylabel('log10(error)')
#legend()
#show()
