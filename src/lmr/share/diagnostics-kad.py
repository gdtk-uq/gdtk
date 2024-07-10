#! /usr/bin/env python
#
# A script that plots the contents of lmrsim/diagnostics/nk-diagnostics
#
# RBO 2024-07-04: animate convergence live
#
# author: Kyle A. Damm
# date: 2024-03-26

from numpy import asarray
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def import_data(filename):
    # read in data from file
    with open(filename) as fp:
        lines = [line.strip() for line in fp]

    # separate the header and the actual data
    keys = []
    values = []
    for n,line in enumerate(lines):
        if n == 0:
            keys = line.split()
        else:
            values.append(asarray(line.split(),dtype=float))

    # the data is currently stored per line, unzip the data to store it per variable
    cols = list(zip(*values))
    # create a dictionary with the keys and the unzipped data
    data = dict(zip(keys,cols))

    return data

def fig_anim(i):
    data = import_data(diag_file)
    p1.set_data(data["step"], data['relaxation-factor'])
    p2.set_data(data["step"], data['linear-solve-residual-rel-achieved'])
    p3.set_data(data["step"], data['n-iters'])
    p4.set_data(data["step"], data['CFL'])
    p5.set_data(data["step"], data['global-residual-rel'])
    p6.set_data(data["step"], data['mass-balance'])

    ax.set_xlim(0, int(data["step"][-1]*2))

if __name__=='__main__':
    import sys
    fig = plt.figure(figsize=(8, 4.5))
    plt.rcParams.update({'font.size': 10})

    ax = plt.gca()
    twin1 = ax.twinx()
    twin2 = ax.twinx()

    fig.subplots_adjust(right=0.75)

    # Offset the right spine of twin2
    twin2.spines['right'].set_position(("axes", 1.125))

    # plot container contents
    diag_file = "lmrsim/diagnostics/nk-diagnostics"
    data = import_data(diag_file)
    global p1, p2, p3, p4, p5, p6
    p1, = ax.plot(data["step"], data['relaxation-factor'], color='red', linewidth=3, alpha = 0.5, label='Relaxation Factor')
    p2, = ax.semilogy(data["step"], data['linear-solve-residual-rel-achieved'], color='blue', alpha = 0.5, linewidth=3, label='Linear Solve Relative Residual')
    p3, = twin1.plot(data["step"], data['n-iters'], color='orange', linewidth=2, alpha = 0.5, label='Krylov Vectors')
    p4, = twin2.semilogy(data["step"], data['CFL'], color='green', linewidth=2, alpha = 0.5, label='CFL')
    p5, = ax.semilogy(data["step"], data['global-residual-rel'], color='black', linewidth=3, label='Global Relative Residual')
    p6, = ax.semilogy(data["step"], data['mass-balance'], color='purple', linewidth=2, alpha=0.5, label='Mass Balance')

    # set axis limits
    ax.set_xlim(0, int(data["step"][-1]*2))
    ax.set_ylim(1e-14, 100)
    twin1.set_ylim(0, 100)
    twin2.set_ylim(1e-2, 1e6)
    twin2.set_yscale("log")

    # set axis labels
    #ax.set_xlabel(xLabel)
    ax.set_ylabel("Global Relative Residual")
    twin1.set_ylabel("Krylov Vectors")
    twin2.set_ylabel("CFL")

    # colour the axes
    twin1.yaxis.label.set_color(p3.get_color())
    twin2.yaxis.label.set_color(p4.get_color())

    # colour the ticks
    tkw = dict(size=4, width=1.5)
    twin1.tick_params(axis='y', colors=p3.get_color(), **tkw)
    twin2.tick_params(axis='y', colors=p4.get_color(), **tkw)

    # some final formatting
    ax.legend(handles=[p1, p2, p3, p4, p5, p6], loc='upper right', frameon=True, framealpha=1, edgecolor="black")
    ax.grid()
    if (sys.argv[1] == "live"):
        anim = FuncAnimation(fig, fig_anim, interval=1000, cache_frame_data=False)
        plt.show()
    else:
        # assume two arguments passed
        plt.savefig(sys.argv[2], dpi=600)

