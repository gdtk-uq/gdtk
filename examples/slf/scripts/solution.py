"""
Read and plot an SLF solution.

@author: Nick Gibbons
"""

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from random import seed,choice
from sys import argv

from gdtk.slf import read_solution_file
from gdtk.gas import GasModel

plt.rcParams['font.family'] = 'serif'
plt.rcParams['svg.fonttype'] = 'none'

seed(100)
colours = [choice(list(mcolors.CSS4_COLORS.keys())) for i in range(20)]

if __name__=='__main__':
    file = argv[1]
    gm = GasModel("gm.lua")
    data = read_solution_file(file)

    fig = plt.figure(figsize=(10,4.5))
    ax0,ax1 = fig.subplots(1,2)

    lines = []
    for i,name in enumerate(gm.species_names):
        line = ax1.semilogy(data['Z'], data['Y'][:,i], color=colours[i], label=name, linewidth=1.5)[0]
        lines.append(line)
        offset = 6
        y = data['Y'][-1,i]
        if y<1e-10: continue
        ax1.annotate('<'+name, xy=(1,y), xytext=(offset,0), color=line.get_color(), 
                    xycoords = ax1.get_yaxis_transform(), textcoords="offset points",
                    size=14, va="center")

    #fig.subplots_adjust(bottom=0.3, wspace=0.2)
    #ax0.plot(init['Z'], init['T'], 'b-')
    ax0.plot(data['Z'], data['T'], 'r-', linewidth=1.5)
    ax0.set_ylabel('Temperature (K)')
    ax0.set_xlabel('Mixture Fraction')
    ax1.set_ylabel('Mass Fraction')
    ax1.set_xlabel('Mixture Fraction')
    ax1.set_ylim((1e-10, 1e0))
    #ax1.legend(handles = lines , labels=gm.species_names,loc='upper center', 
    #            bbox_to_anchor=(-0.2, -0.2),fancybox=False, shadow=False, ncol=7)
    ax0.grid()
    ax1.grid()
    plt.tight_layout()
    plt.show()

