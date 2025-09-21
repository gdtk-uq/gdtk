#! /usr/bin/env python
#
# author: Kyle A. Damm
# date:   02-10-2025
#

import sys
import os
import shutil
import math
import matplotlib.pyplot as plt

# Freestream conditions from Mabey report
rho0  = 0.177131    # kg/m^3
vel0  = 712.9       # m/s
T0    = 62.157      # K
mach0 = 4.5099      # -

def final_loads_dir(dir_path):
    final_dir = next(os.walk(dir_path))[1]
    return dir_path+ '/' + str(max(final_dir))

def gather_surface_data(dir_path, container, headers):
    loads_files = os.listdir(dir_path)
    for filename in loads_files:
        file_path = dir_path + '/' + filename
        print("Found file: ", file_path)
        read_loads_file(file_path, container, headers)

def read_loads_file(filename, container, headers):
    with open(filename, 'r') as f:
        data = f.readlines()
        count = 0
        for line in data:
            dat = line.split()
            # pull out the headers once
            if count == 0 and len(headers) == 0:
                for idx in range(0,len(dat)):
                    var_name = dat[idx]
                    headers.append(var_name)
            elif count > 0:
                arr = []
                for entry in dat:
                    arr.append(float(entry))
                container.append(arr)
            else:
                print("nothing to do on this line")
            count += 1    

def extract_line(output, p0, p1, N):
    cmd = 'lmr extract-line -l "{},{},{},{},{},{},{}" -f -o="{}" --add-vars=mach'.format(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], N, output)
    os.system(cmd)

def read_validation_file(filename, xidx, yidx):
    xvals = []
    yvals = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):  # skip empty lines and comments
                continue
            dat = line.split()
            xvals.append(float(dat[xidx]))
            yvals.append(float(dat[yidx]))
    return xvals, yvals

def compute_skin_friction(loads, headers):
    posx = []; tauw = [];
    xidx = headers.index('pos.x');
    posx = [x[xidx] for x in loads];
    yidx = headers.index('tau_wall.x')
    tauw = [x[yidx] for x in loads]
    posx, tauw = zip(*sorted(zip(posx, tauw)))
    cf = [2 * t / (rho0 * vel0**2) for t in tauw]
    return posx, cf

def sample_boundary_layer(p0, p1, N, outputFileName, var, var0):
    data = []; headers = []
    extract_line(outputFileName, p0, p1, N)
    read_loads_file(outputFileName, data, headers)
    data = sorted(data, key=lambda x: x[1]) # sort data by pos.y
    posy = []; xval = [];
    xidx = headers.index(var);
    xval = [x[xidx] for x in data];
    xval = [x/var0 for x in xval]
    yidx = headers.index('pos.y')
    posy = [x[yidx] for x in data]
    posy = [p0[1]-y for y in posy]
    posy, xval = zip(*sorted(zip(posy, xval)))
    return posy,xval

def plot(xcfd, ycfd, xexp, yexp, xmax, ymax, xerr, yerr, xlabel, ylabel):
    plt.errorbar(xexp, yexp, xerr=xerr, yerr=yerr, fmt='o', color='black', ecolor='black', elinewidth=1, capsize=3, label="Mabey exp.")
    plt.plot(xcfd, ycfd, color='black', label="Eilmer")
    plt.xlim(0.0, xmax)
    plt.ylim(0.0, ymax)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True, color='black', linestyle='--', linewidth=0.5)
    plt.show()

if __name__=='__main__':

    # define paths
    path2loads = "./lmrsim/loads"
    path2exp   = "./validation_data"

    # plot skin-friction along length of plate
    data = []; headers = []
    path = final_loads_dir(path2loads)
    gather_surface_data(path, data, headers)
    data = sorted(data, key=lambda x: x[0]) # sort data by pos.x
    posx,cf = compute_skin_friction(data, headers)
    xexp, yexp = read_validation_file(path2exp+"/skin_friction.dat", 0, 1)
    plot(posx, cf, xexp, yexp, 1.4, 0.002, 0.0, 0.0001, 'xloc [m]', '$C_f$')

    # plot boundary layer slice (slice from p0 to p1 with N samples)
    p0 = [0.368, 0.56, 0.0]
    p1 = [0.358, 0.5, 0.0]
    N  = 1000
    outputFileName = 'boundary_layer_slice.dat'
    # Mach number
    posy, xval = sample_boundary_layer(p0, p1, N, outputFileName, 'M_local', mach0)
    xexp, yexp = read_validation_file(path2exp+"/bl_profile.dat", 5, 1)
    plot(xval, posy, xexp, yexp, 1.5, 0.02, 0.0, 0.0, 'M/M0', 'yloc [m]')
    # velocity
    posy, xval = sample_boundary_layer(p0, p1, N, outputFileName, 'vel.x', vel0)
    xexp, yexp = read_validation_file(path2exp+"/bl_profile.dat", 6, 1)
    plot(xval, posy, xexp, yexp, 1.5, 0.02, 0.0, 0.0, 'U/U0', 'yloc [m]')
    # temperature
    posy, xval = sample_boundary_layer(p0, p1, N, outputFileName, 'T', T0)
    xexp, yexp = read_validation_file(path2exp+"/bl_profile.dat", 7, 1)
    plot(xval, posy, xexp, yexp, 5, 0.02, 0.0, 0.0, 'T/T0', 'yloc [m]')
