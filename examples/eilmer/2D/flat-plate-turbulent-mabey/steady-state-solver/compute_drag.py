#! /usr/bin/env python
# script for calculating drag force acting on a flat plate
# KD, 2019-12-17

from math import *
import os

# read last tindx
timesFile = "loads/mabey-loads.times"
tindx = 0
with open(timesFile, 'r') as f:
    data = f.readlines()
    for line in data:
        dat = line.split()
        if (dat[0] != "#"):
            tindx = dat[0]

print("computing drag at tindx: {}".format(tindx))
            
# data structures
loadsFile = "loads/t" + tindx + "/surface.data"
p     =   [] # pressure (Pa) 
A     =   [] # surface element area (m^2)
nx    =   [] # x-component of normal vector
tau_w =   [] # wall shear stress (Pa)
lx    =   [] # x-direction cosine for shear stress
outsign = [] # interface outsign
n     = 0  # number of surface elements

# collate all surface files into one file
cmd = "cat loads/t" + tindx + "/*dat > loads/t" + tindx + "/surface.data"
os.system(cmd)

# read surface data file
with open(loadsFile, 'r') as f:
    data = f.readlines()
    idx = 0
    for line in data:
        dat = line.split()
        if (dat[0] != "#"):
            dat = line.split()
            p.append(float(dat[11]))
            A.append(float(dat[3]))
            nx.append(float(dat[12]))
            tau_w.append(float(dat[7]))
            lx.append(float(dat[8]))
            outsign.append(float(dat[19]))
            n += 1

# compute force along flat plate
force_i = 0.0 # N
force_v = 0.0 # N
for i in range(0,n):
    force_i += p[i]*A[i]*(nx[i]*outsign[i])     # inviscid contribution (should be 0 N)
    force_v += tau_w[i]*A[i]*(lx[i]*outsign[i]) # viscous contribution
force_t = force_i + force_v # total drag force acting on the plate

# output force
print("F_i: {} N".format(abs(force_i)))
print("F_v: {} N".format(abs(force_v)))
print("F_t: {} N".format(abs(force_t)))

# remove surface data file
cmd = "rm loads/t" + tindx + "/surface.data"
os.system(cmd)
