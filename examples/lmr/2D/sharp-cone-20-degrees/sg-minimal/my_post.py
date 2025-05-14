#! /usr/bin/env python3
#
# my_post.py
# Custom postprocessing for an Eilmer5(lmr) simulation.
#
# Peter J. 2025-05-14

import os
from gdtk.lmr import LmrConfig, SimInfo

lmrcfg = LmrConfig()
# print("content=", lmr.data)
print("sim_dir=", lmrcfg["simulation-directory"],
      "vtk_dir=", lmrcfg["vtk-output-directory"])

sim = SimInfo(lmrcfg)
print("blks=", sim.blocks)
print("times=", sim.times)
print("snapshots=", sim.snapshots)
print("variables=", sim.fluid_variables)
print("grid tags=", [sim.grids[i]['tag'] for i in range(len(sim.grids))])
