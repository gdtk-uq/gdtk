#! /usr/bin/env python
#
# Script that plots the contents of the residual file e4-nk.diagnostics.dat
#
# Kyle A. Damm [06-06-2022]

from numpy import asarray
import matplotlib.pyplot as plt
import time

def read_diagnostics_file(filename):
    with open(filename) as fp:
        lines = [line.strip() for line in fp]

    keys = []
    values = []

    for line in lines:
        if line.startswith('#'):
            keys.append(line)
        else:
            values.append(line)

    keys = [k.split(':')[-1].strip() for k in keys]
    values = [asarray(list(line.split()),dtype=float) for line in values]
    cols = list(zip(*values))
    data = dict(zip(keys,cols))
    return data

def generate_plot():

    # turn on interactive mode
    plt.ion()

    # read in data to generate initial plot
    data = read_diagnostics_file('e4-nk.diagnostics.dat')

    fig, ax = plt.subplots()
    fig.subplots_adjust(right=0.75)

    twin1 = ax.twinx()
    twin2 = ax.twinx()

    # Offset the right spine of twin2
    twin2.spines['right'].set_position(("axes", 1.1))

    p1, = ax.semilogy(data['step'], data['global-residual-rel'], color='black', linewidth=3, label='global-residual-rel')
    p2, = ax.semilogy(data['step'], data['linear-solve-residual'], color='blue', label='linear-solve-residual')
    p3, = ax.semilogy(data['step'], data['mass-balance'], color='purple', label='mass-balance')
    p4, = ax.semilogy(data['step'], data['omega'], color='red', label='omega')
    p5, = twin1.plot(data['step'], data['nIters'], color='orange', label='Krylov vectors')
    p6, = twin2.semilogy(data['step'], data['CFL'], color='green', label='CFL')

    ax.set_xlim(0, int(data['step'][-1]+50))
    ax.set_ylim(1e-14, 100)
    twin1.set_ylim(0, 100)
    twin2.set_ylim(1e-2, 1e6)

    ax.set_xlabel("Nonlinear Step")
    ax.set_ylabel("Relative Residual")
    twin1.set_ylabel("Krylov Vectors")
    twin2.set_ylabel("CFL")

    # Colour the axes
    twin1.yaxis.label.set_color(p5.get_color())
    twin2.yaxis.label.set_color(p6.get_color())

    # Colour the ticks
    tkw = dict(size=4, width=1.5)
    twin1.tick_params(axis='y', colors=p5.get_color(), **tkw)
    twin2.tick_params(axis='y', colors=p6.get_color(), **tkw)

    ax.legend(handles=[p1, p2, p3, p4, p5, p6], loc='upper right')
    #ax.legend(loc='center left', bbox_to_anchor=(1.2, 0.5))

    ax.grid()
    plt.show()

    # while loop to continuously plot the data
    while True:
        data = read_diagnostics_file('e4-nk.diagnostics.dat')
        ax.set_xlim(0, int(data['step'][-1]+50))
        p1.set(xdata = data['step'], ydata = data['global-residual-rel'])
        p2.set(xdata = data['step'], ydata = data['linear-solve-residual'])
        p3.set(xdata = data['step'], ydata = data['mass-balance'])
        p4.set(xdata = data['step'], ydata = data['omega'])
        p5.set(xdata = data['step'], ydata = data['nIters'])
        p6.set(xdata = data['step'], ydata = data['CFL'])
        plt.draw()
        #plt.savefig('residuals.png')
        plt.pause(1)

if __name__=='__main__':
    generate_plot()

