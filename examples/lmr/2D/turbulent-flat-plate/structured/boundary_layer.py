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
    #directory_names = argv[1:]
    directory_names = ['80um','01um']
    colours = ['olive','blue']

    datas = []

    for directory_name in directory_names:
        print("Reading name: ", directory_name)
        data = read_loads_file(directory_name+'/line.txt')
        datas.append(data)

    fig = plt.figure(figsize=(11,4.5))
    axes = fig.subplots(1,2) 
    lims = [(-305,2200),(-320,2200)]
    lines = []

    for i,data in enumerate(datas):
        ax = axes[i]
        colour = colours[i]
        dirs = directory_names[i]
        
        axy = ax.twiny()
        lines.extend(ax.plot(data['vel.x'], data['pos.y']*1000, color=colour, marker='.',mfc='none', linewidth=1.0, linestyle='--', label=dirs+'- vel'))
        lines.extend(axy.plot(data['T'], data['pos.y']*1000, color=colour,  marker='o', mfc='none', linewidth=1.5, linestyle='-', label=dirs+'- T'))
        ax.plot([0.0,0.0], [data['pos.y'].min()*1000, data['pos.y'].max()*1000], 'k--')

        ax.set_ylabel('Y Position Above Plate (mm)')
        ax.set_xlabel('X-Velocity (m/s)')
        axy.set_xlabel('Temperature (K)')
        ax.set_xlim(lims[i])
        #ax.grid()


    leg = ax.legend(lines, ['', '', 'Velocity', 'Temperature'],
                    title="80um    01um", columnspacing=1.5,
                    ncol=2, framealpha=1.0, loc="upper center")
    leg._legend_box.align = "left"
    plt.tight_layout()
    plt.savefig('bl.svg')
    #plt.close()
    plt.show()

