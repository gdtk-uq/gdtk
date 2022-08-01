#! /usr/bin/env python3
# plot_results.py
# Peter J, Daniel Smith and Wilson Chan.

import numpy as np
import matplotlib.pyplot as plt
from math import atan, degrees
from gdtk.geom.vector3 import Vector3
from gdtk.geom.path import Bezier

if __name__ == "__main__":
    print("Plot the contours")
    plt.figure()
    profile_axes = plt.subplot(1,1,1,
                               xlim=(0.0, 0.50), ylim=(0.0, 0.05),
                               ylabel="radial position, m",
                               xlabel="axial position, m",
                               title="Nozzle profiles")
    control_points = np.loadtxt("Bezier-control-pts-t4m4b-initial.data", skiprows=1)
    bez_init = Bezier([Vector3(p[0], p[1]) for p in control_points])
    ts = np.linspace(0,1,100)
    ps = [bez_init(t) for t in ts]
    profile_axes.plot([p.x for p in ps], [p.y for p in ps], label="init")
    #
    control_points = np.loadtxt("Bezier-control-pts-t4m4b.opt.data", skiprows=1)
    bez_opt = Bezier([Vector3(p[0], p[1]) for p in control_points])
    ps = [bez_opt(t) for t in ts]
    profile_axes.plot([p.x for p in ps], [p.y for p in ps], label="opt")
    plt.legend()
    #
    print("Plot the outflow")
    plt.figure()
    mach_axes = plt.subplot(2,1,1,
                            xlim=(0.0, 0.05), ylim=(3.95, 4.05),
                            ylabel="Mach number",
                            title="Optimized nozzle exit flow")
    theta_axes = plt.subplot(2,1,2,
                             xlim=(0.0, 0.05), ylim=(-0.05, +0.05),
                             xlabel="radial position, m",
                             ylabel="Flow angle, degrees")
    exit_data = np.loadtxt("t4m4b-exit-initial.data", skiprows=1)
    Ms = exit_data[:,18]
    ys = exit_data[:,1]
    velxs = exit_data[:,5]
    velys = exit_data[:,6]
    thetas = [degrees(atan(vy/vx)) for vx,vy in zip(velxs,velys)]
    mach_axes.plot(ys, Ms, label="init")
    theta_axes.plot(ys, thetas, label="init")
    #
    exit_data = np.loadtxt("t4m4b-exit.data", skiprows=1)
    Ms = exit_data[:,18]
    ys = exit_data[:,1]
    velxs = exit_data[:,5]
    velys = exit_data[:,6]
    thetas = [degrees(atan(vy/vx)) for vx,vy in zip(velxs,velys)]
    mach_axes.plot(ys, Ms, label="opt")
    theta_axes.plot(ys, thetas, label="opt")
    mach_axes.legend()
    theta_axes.legend()
    plt.show()
