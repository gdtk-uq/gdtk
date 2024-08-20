"""
Script for plotting boundary layer slices from Eilmer.

@author: Nick Gibbons (20/08/24)
"""

from numpy import array, argsort, concatenate, mean, argmin, log, diff, log10
from sys import argv
from glob import glob
import matplotlib.pyplot as plt
from wall_heat_transfer import read_loads_file

if __name__=='__main__':
    directory_names = argv[1:]

    datas = []

    for directory_name in directory_names:
        print("Reading name: ", directory_name)
        data = read_loads_file(directory_name+'/line.txt')
        datas.append(data)

    fig = plt.figure(figsize=(11,4.5))
    ax,ax2 = fig.subplots(1,2) 
    #colours = ['blue','red','darkgreen','magenta', 'goldenrod', 'teal', 'wheat']
    colours = ['blue', 'olive']

    for j, dirs in enumerate(directory_names):
        data = datas[j]
        colour = colours[j]
        
        ax.plot(data['vel.x'], data['pos.y']*1000, color=colour, marker='o',mfc='none', linewidth=1.5, linestyle='-', label=dirs)
        ax2.plot(data['T'], data['pos.y']*1000, color=colour,  marker='o', mfc='none', linewidth=1.5, linestyle='-')

    ax.legend(framealpha=1.0)
    ax.set_ylabel('Y Position Above Plate (mm)')
    ax.set_xlabel('X-Velocity (m/s)')
    ax2.set_xlabel('Temperature (K)')
    ax.grid()
    ax2.grid()
    plt.tight_layout()
    plt.savefig('wht.svg')
    #plt.close()
    plt.show()

