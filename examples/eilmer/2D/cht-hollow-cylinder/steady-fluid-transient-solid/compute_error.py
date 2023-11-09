#! /usr/bin/env python
#
# This script reads in two Eilmer loads files and computes the RMS error between the surface temperature profiles.
#
# author: Kyle A. Damm
# date:   2023-11-09

import os as os
import numpy as np

def final_loads_dir(dir_path):
    dir_list = os.listdir(dir_path)
    dir_list.sort(reverse=True)
    final_dir = ''
    for subdir in dir_list:
        subdir_path = dir_path + subdir
        if (os.path.isdir(subdir_path)):
            final_dir = subdir_path
            break
    return final_dir

def gather_surface_data(dir_path, fname):
    cwd = os.getcwd() # current working directory
    os.chdir(dir_path)
    cmd = 'cat *dat > ' + fname
    os.system(cmd) # bundle block files together
    os.chdir(cwd)
    cmd = 'mv ' + dir_path + '/' + fname + ' .'
    os.system(cmd) # move file to current working directory

def read_file(filename, container, headers):
    with open(filename, 'r') as f:
        data = f.readlines()
        line_num = 0
        for line in data:
            dat = line.split()
            # pull out the headers once
            if (dat[0] == '#' and line_num >= 1 and len(headers) == 0):
                for idx in range(1,len(dat)):
                    var_name = dat[idx].split(':')[1]
                    headers.append(var_name)
            elif (dat[0] != '#'):
                arr = []
                for entry in dat:
                    arr.append(float(entry))
                container.append(arr)
            line_num += 1

def compare_temperature(data1, data2, headers1, headers2):
    idx = headers1.index("T");
    val1 = np.array([x[idx] for x in data1])
    idx = headers2.index("T");
    val2 = np.array([x[idx] for x in data2]);
    assert len(val1) == len(val2), "The Eilmer solution and reference solution files have different lengths."
    diff = val1 - val2
    n = len(val1)
    error = np.sqrt(np.dot(diff,diff)/n)
    return error

if __name__=='__main__':
    import sys
    sol_fname = "eilmer_solution.dat"
    ref_fname = "reference_solution.dat"
    path = final_loads_dir("./loads/")
    gather_surface_data(path, sol_fname)
    sol_data = []; sol_headers = []
    read_file(sol_fname, sol_data, sol_headers)
    ref_data = []; ref_headers = []
    read_file(ref_fname, ref_data, ref_headers)
    error = compare_temperature(sol_data, ref_data, sol_headers, ref_headers)
    print("Surface Temperature RMS:", error)
    
