#! /usr/bin/env python
#
# Script to plot the normalized heat transfer from the final loads file from an Eilmer simulation.
#
# @author: Kyle A. Damm (2024-08-23)
#

import sys
import os
import math
import matplotlib.pyplot as plt

def final_loads_dir(dir_path):
    final_dir = next(os.walk(dir_path))[1]
    return dir_path+ '/' + str(max(final_dir))

def gather_surface_data(dir_path, container, headers):
    loads_files = os.listdir(dir_path)
    for filename in loads_files:
        file_path = dir_path + '/' + filename
        print(file_path)
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

def read_reference_file(filename, container):
    with open(filename, 'r') as f:
        data = f.readlines()
        for line in data:
            dat = line.split()
            if len(dat) > 0 and dat[0] != "#":
                container[0].append(float(dat[0]))
                container[1].append(float(dat[1]))

def evaluate_angle_from_stagnation_point(container, headers):
    xtmp = [x[headers.index("pos.x")] for x in container]
    ytmp = [x[headers.index("pos.y")] for x in container]
    xval = []
    for i in range(len(xtmp)):
        angle = math.atan(abs(ytmp[i]/xtmp[i])) * 180.0/math.pi
        xval.append(angle)
    return xval

def access_heat_transfer(container, headers):
    yidx = headers.index("q_total")
    yval = [abs(x[yidx]) for x in container]
    ymax = max(yval)
    for i in range(0,len(yval)):
        yval[i] /= ymax
    return yval

def plot_normalized_heat_transfer(container1, container2, headers, exp_data, theory_data):

    # compute and store the angle from stagation point in degrees
    xval1 = evaluate_angle_from_stagnation_point(container1, headers)
    xval2 = evaluate_angle_from_stagnation_point(container2, headers)

    # pull out the heat transfer
    yval1 = access_heat_transfer(container1, headers)
    yval2 = access_heat_transfer(container2, headers)

    # plot the data
    plt.plot(xval1, yval1, color='blue', label="GMSH grid")
    plt.plot(xval2, yval2, color='red', label="Native grid")
    plt.plot(theory_data[0], theory_data[1], color='black', label="Kemp theory")
    plt.scatter(exp_data[0], exp_data[1], color='black', label="Exp. data")
    plt.xlabel("Angle from stagnation point, degrees")
    plt.ylabel("Normalized heat flux")
    plt.grid(True, color='black', linestyle='--', linewidth=0.5)
    plt.legend()
    plt.show()

if __name__=='__main__':

    # path to Eilmer simulation loads directory
    gmsh_path = final_loads_dir("gmsh-grid/lmrsim/loads")
    native_path = final_loads_dir("native-grid/lmrsim/loads")

    # gather reference data
    kemp_exp = [[],[]]
    read_reference_file("./reference-data/kemp_experiment.dat", kemp_exp)
    kemp_theory = [[],[]]
    read_reference_file("./reference-data/kemp_theory.dat", kemp_theory)

    # gather Eilmer data
    native_data = []; gmsh_data = []; headers = []
    # gmsh grid
    gather_surface_data(gmsh_path, gmsh_data, headers)
    gmsh_data = sorted(gmsh_data, key=lambda x: x[0]) # sort data by xloc
    # native grid
    gather_surface_data(native_path, native_data, headers)
    native_data = sorted(native_data, key=lambda x: x[0]) # sort data by xloc

    # generate plot
    plot_normalized_heat_transfer(gmsh_data, native_data, headers, kemp_exp, kemp_theory)
