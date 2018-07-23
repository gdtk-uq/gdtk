#!/usr/bin/env python
"""
Plot surface Mach number of Kiock planar turbine cascade.

Peter Blyton
    11 July 2011: Original implementation, rotate from x coord to normalised chord position.

"""
from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

# currently an issue with matplotlib and latex fonts, removed for time being.
#rc("text", usetex=True) # LaTeX on plot, requires "dvipng" package in Ubuntu
#rc("font", family="serif")

STAGGER_ANGLE = (33.3/180.0)*np.pi

def plot_surface_Mach():
    sim_data = np.loadtxt("simulation_results.dat") # simulation data
    RG_data = np.loadtxt("kiock_experimental_RG.dat") # data from Rhode-St.-Genese wind tunnel
    GO_data = np.loadtxt("kiock_experimental_GO.dat") # data from Goettingen wind tunnel
    BS_data = np.loadtxt("kiock_experimental_BS.dat") # data from Braunschweig wind tunnel
    OX_data = np.loadtxt("kiock_experimental_OX.dat") # data from Oxford wind tunnel
    
    # rotate the x and y coords into normalised chord position reference frame
    chord_pos = sim_data[:, 0]*np.cos(STAGGER_ANGLE) - sim_data[:, 1]*np.sin(STAGGER_ANGLE)
    
    # calculate surface Mach number with isentropic flow relation
    # surface_Mach  = np.sqrt((2.0/(gma - 1.0))*((p_tot/sim_data[:, 8])**((gma - 1.0)/gma) - 1.0))
    surface_Mach = sim_data[:,20] #M_local is computed by post processor and stored in 20th column
    
    plt.plot(chord_pos, surface_Mach, "x")
    plt.plot(RG_data[:, 0], RG_data[:, -1], "o w")
    plt.plot(GO_data[:, 0], GO_data[:, -1], "s w")
    plt.plot(BS_data[:, 0], BS_data[:, -1], "^ w")
    plt.plot(OX_data[:, 0], OX_data[:, -1], "v w")
    plt.xlabel("Chord position (x/C)")
    plt.ylabel("Surface Mach number")
    plt.axis(xmin=-0.2, xmax=1.2, ymin=-0.2, ymax=1.2)
    plt.suptitle("Kiock planar turbine cascade")
    plt.legend(("CO2 SW", "RG", "GO", "BS", "OX"), loc="upper left")
    plt.savefig("surface_Mach.pdf")
    plt.show()
    
    return 0

if __name__ == "__main__":
    plot_surface_Mach()
    print("Surface Mach number plotted.")

