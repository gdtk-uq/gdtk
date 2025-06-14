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
# print("snap_last=", snap_last)
# Sample the data along a grid-line for block 1,
# which happens to be over the ramp.
xs, ys, ps = snap_last.get_slice(var='p', j=10)
# xs = snap_last.fields[1]['pos.x'][:,10]
# ps = snap_last.fields[1]['p'][:,10]
# print("x,p=", xs, ps)
fig, ax = plt.subplots()
ax.plot(xs, ps/1000)
ax.set_xlabel('x, m')
ax.set_ylabel('pressure, kPa')
plt.show()

# Load the PyVista data via VTK files.
pvdata = sim.load_pvd_into_pyvista()
print("pvdata=", pvdata)

