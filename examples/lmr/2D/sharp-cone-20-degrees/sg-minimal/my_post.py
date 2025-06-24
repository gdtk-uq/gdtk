#! /usr/bin/env python3
#
# my_post.py
# Custom postprocessing for an Eilmer5(lmr) simulation.
#
# Peter J. 2025-05-14

import os
from gdtk.lmr import LmrConfig, SimInfo
import matplotlib.pyplot as plt

lmrcfg = LmrConfig()
# print("content=", lmr.data)
print("sim_dir=", lmrcfg["simulation-directory"],
      "vtk_dir=", lmrcfg["vtk-output-directory"])

sim = SimInfo(lmrcfg)
# print("sim.sim_cfg=", sim.sim_cfg)
print("grids_info=", sim.grids)
print("gridarrays_info=", sim.gridarrays)
print("blocks_info=", sim.blocks)
print("times=", sim.times)
print("snapshots=", sim.snapshots)
print("variables=", sim.fluid_variables)
print("grid tags=", [sim.grids[i].tag for i in range(len(sim.grids))])
print("gas_model=", sim.gas_model)

# Read grids from the grid-preparation stage.
grids = sim.read_grids()
for i in range(len(grids)):
    grids[i].write_to_vtk_file("grid-%d.vtk" % (i))

# Read a couple of snapshots
snap_first = sim.read_snapshot('0000')
snap_last = sim.read_snapshot(sim.snapshots[-1])

# Add a derived variable.
snap_last.add_vars(['mach','pitot','total_p','total_h','total_T'])

xs, ys, ps = snap_last.get_slice(var='p', j=10)
xs, ys, Ms = snap_last.get_slice(var='mach', j=10)
xs, ys, pitot = snap_last.get_slice(var='pitot', j=10)
xs, ys, total_h = snap_last.get_slice(var='total_h', j=10)
# print("x,p=", xs, ps)
fig, ax = plt.subplots(4,1)
ax[0].plot(xs, ps/1000)
ax[1].plot(xs, Ms)
ax[2].plot(xs, pitot/1000)
ax[3].plot(xs, total_h/1.0e6)
ax[3].set_xlabel('x, m')
ax[0].set_ylabel('pressure, kPa')
ax[1].set_ylabel('Mach number')
ax[2].set_ylabel('Pitot pressure, kPa')
ax[3].set_ylabel('Total enthalpy, MJ/kg')
plt.show()

# Load the PyVista data via VTK files.
pvdata = sim.load_pvd_into_pyvista()
print("pvdata=", pvdata)

