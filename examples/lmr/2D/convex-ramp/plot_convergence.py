#! /usr/bin/env python
#
# A script that plots the contents of lmrsim/diagnostics/nk-diagnostics
#
# author: Kyle A. Damm
# date: 2024-03-26
#
from numpy import asarray
import matplotlib.pyplot as plt

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

def generate_plot(container, xVar, xLabel):

    # set some plot settings
    plt.rcParams['figure.figsize'] = [16, 8]
    plt.rcParams.update({'font.size': 16})

    fig, ax = plt.subplots()
    fig.subplots_adjust(right=0.75)

    twin1 = ax.twinx()
    twin2 = ax.twinx()

    # offset the right spine of twin2
    twin2.spines['right'].set_position(("axes", 1.125))

    # plot container contents
    p1, = ax.semilogy(container[xVar], container['global-residual-rel'], color='black', linewidth=3, label='Global Relative Residual')
    p2, = twin1.plot(container[xVar], container['n-iters'], color='blue', linewidth=2, label='Krylov Vectors')
    p3, = twin2.semilogy(container[xVar], container['CFL'], color='green', linewidth=2, label='CFL')

    # set axis limits
    ax.set_xlim(0, int(container[xVar][-1]+50))
    ax.set_ylim(1e-14, 100)
    twin1.set_ylim(0, 100)
    twin2.set_ylim(1e-2, 1e6)

    # set axis labels
    ax.set_xlabel(xLabel)
    ax.set_ylabel("Global Relative Residual")
    twin1.set_ylabel("Krylov Vectors")
    twin2.set_ylabel("CFL")

    # colour the axes
    twin1.yaxis.label.set_color(p2.get_color())
    twin2.yaxis.label.set_color(p3.get_color())

    # colour the ticks
    tkw = dict(size=4, width=1.5)
    twin1.tick_params(axis='y', colors=p2.get_color(), **tkw)
    twin2.tick_params(axis='y', colors=p3.get_color(), **tkw)

    # illustrate the phases by fill colour
    phase1 = ax.fill_between(container[xVar], container['global-residual-rel'], where=asarray(container['phase']) <= 1, facecolor='red', alpha=.2, label="Phase #1")
    phase2 = ax.fill_between(container[xVar], container['global-residual-rel'], where=asarray(container['phase']) >= 1, facecolor='red', alpha=.4, label="Phase #2")

    # some final formatting
    ax.legend(handles=[phase1, phase2], loc='upper right', frameon=False)
    ax.grid()

    fig.savefig('convergence_vs_'+xVar+'.png')
    #plt.show()

if __name__=='__main__':

    fileDir = './lmrsim/diagnostics/'
    fileName = 'nk-diagnostics'
    data = import_data(fileDir+fileName)
    generate_plot(data, 'step', 'Nonlinear Step')
    generate_plot(data, 'wall-clock', 'Wall-clock (s)')
