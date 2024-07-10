#!/usr/bin/python3
"""
Plot the eilmer 4 newton-krylov solver residuals, interactively.

@author: Nick Gibbons
"""

from numpy import log10, array
from random import seed, sample
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from time import sleep
from datetime import datetime

RESIDUALS_FILE = "lmrsim/diagnostics/nk-diagnostics"

def number(thing):
    try:
        outnumber = int(thing)
    except(ValueError):
        try:
            outnumber = float(thing)
        except(ValueError):
            outnumber = complex(thing.replace('i','j')).real
    return outnumber

def read_diagnostics_file(filename):
    with open(filename) as fp:
        lines = [line.strip() for line in fp]

    keys = []
    values = []

    # first line contins column headers
    keys = lines[0].split()
    # remaining lines have data
    
    for line in lines[1:]:
        values.append(line)

    #keys = [k.split(':')[-1].strip() for k in keys]
    values = [list(map(number, line.split())) for line in values]
    cols = [array(col) for col in zip(*values)]
    assert len(cols)==len(keys)
    data = dict(zip(keys,cols))
    return data

def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

today = datetime.today().strftime('%Y-%m-%d')
seed(today)
rgbs = {key:hex_to_rgb(value) for key,value in mcolors.CSS4_COLORS.items()}
luminances = {k:0.2126*v[0] + 0.7152*v[1] + 0.0722*v[2] for k,v in rgbs.items()}
COLOURS = [name for name,luminance in luminances.items() if luminance<220]

class ResidualPlot(object):
    keys = [['linear-solve-residual-achieved', 'mass-balance', 'global-residual-rel'],
            ['n-iters', 'n-restarts', 'relaxation-factor']]
    yscales = ['log', 'linear']

    def __init__(self, filename):
        self.filename = filename
        data = read_diagnostics_file(self.filename)

        self.keys = []
        for group in self.__class__.keys:
            selfgroup = []
            for key in group:
                if key in data.keys():
                    selfgroup.append(key)
            self.keys.append(selfgroup)
        
        nplots = len(self.keys)
        self.fig = plt.figure(figsize=(4.5*nplots, 4.5))
        self.axes = self.fig.subplots(1,nplots, squeeze=False)[0,:]
        colours = sample(COLOURS, sum(len(keylist) for keylist in self.keys))
        self.colours = {key:colours[i] for i,key in enumerate([key for keylist in self.keys for key in keylist])}
        print("Colour palette for {}".format(today))
        print('\n'.join(['{: >21s} : {}'.format(k,v) for k,v in self.colours.items()]))
        
        self.lines = {}
        for i in range(nplots):
            ax = self.axes[i]
            for key in self.keys[i]:
                line = ax.plot(data['step'], data[key], color=self.colours[key], label=key)
                self.lines[key] = line[0]
            lines, labels = ax.get_legend_handles_labels()
            ax.set_yscale(self.yscales[i])
            ax.grid()
            ax.legend(lines, labels, framealpha=1.0, loc=0)
        
        self.fig.tight_layout()
        plt.show()

class LiveResidualPlot(ResidualPlot):
    def __init__(self, filename):
        super().__init__(filename)

    def update_lines(self):
        data = read_diagnostics_file(self.filename)
        for key in self.lines.keys():
            self.lines[key].set_xdata(data['step'])
            self.lines[key].set_ydata(data[key])
        for axes in self.axes: axes.relim()
        for axes in self.axes: axes.autoscale_view()
        self.fig.canvas.draw()
        #axes.legend()


if __name__=='__main__':
    import sys
    
    plt.ion()
    resplot = LiveResidualPlot(RESIDUALS_FILE)
    resplot.fig.tight_layout()
    resplot.fig.canvas.draw()

    if sys.argv[1] == "live":
        while True:
            print(".",end='',flush=True)
            resplot.update_lines()
            plt.pause(0.1)
            sleep(10.0)
    else:
        plt.savefig(sys.argv[2], dpi=600)

