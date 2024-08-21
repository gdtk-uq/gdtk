#! /usr/bin/env python
#
# Script to plot the density profile from an Eilmer snapshot along a line.
#
# @author: Kyle A. Damm (2024-08-21)
#

import os
from numpy import asarray
import matplotlib.pyplot as plt

def post_process(output):
    # line definition
    p0 = [0.0, 0.3, 0.0]
    p1 = [2.0, 0.3, 0.0]
    N = 1000
    
    # post-process Eilmer solution
    cmd = 'lmr extract-line -l "{},{},{},{},{},{},{}" -f -o="{}"'.format(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], N, output)
    os.system(cmd)

def read_file(filename):
    with open(filename) as fp:
        lines = [line.strip() for line in fp]

    keys = []
    values = []

    for idx, line in enumerate(lines):
        if (idx == 0):
            keys.append(line)
        else:
            values.append(line)

    keys = keys[-1].split(' ')
    keys = [k.split(':')[-1].strip() for k in keys]
    values = [asarray(list(line.split()),dtype=float) for line in values]
    cols = list(zip(*values))
    data = dict(zip(keys,cols))
    return data

def generate_plot(data):

    fig, ax = plt.subplots()
    fig.subplots_adjust(right=0.75)

    p1, = ax.plot(data['pos.x'], data['rho'], color='black', linewidth=2)

    ax.set_xlim(0, 2.0)
    ax.set_ylim(0.7, 1.65)

    ax.set_xlabel("x")
    ax.set_ylabel("Density")

    ax.grid()
    plt.show()

if __name__=='__main__':

    # output file name
    output = 'solution_over_line.dat'

    # generate post-process data
    post_process(output)
    
    # data along line
    data = read_file(output)
    generate_plot(data)
