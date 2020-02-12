#! /usr/bin/env python3
# t1d.py
# Transient quasi-one-dimensional gas flow simulation.
#
# Usage:
# $ prep-gas ideal-air.inp ideal-air-gas-model.lua
# $ python3 t1d.py
#
# PJ, 2020-02-12
# 
import math
from eilmer.gas import GasModel, GasState, GasFlow
from t1d_classes import Cell, Face

gmodel = GasModel('ideal-air-gas-model.lua')
flow = GasFlow(gmodel)

print("Set up a collection of cells representing the gas within the duct.")
#     Face  Cell  Face  Cell  Face  Cell  Face
#   ----+-----------+-----------+-----------+----
#       |           |           |           |
#      i-1   i-1    i     i    i+1   i+1   i+2
#       |           |           |           |
#   ----+-----------+-----------+-----------+----
n = 10
cells = [Cell(gmodel) for i in range(n)]
faces = [Face(gmodel) for i in range(n+1)]

print("Quasi-one-dimensional duct is defined by its area at position x.")
def Area(x): return 1.0 # metres^2
dx = 0.01 # metres
for i in range(n+1):
    x = dx * i
    faces[i].x = x
    faces[i].area = Area(x)

print("Set up initial gas states.")
# Here we set up something like the Sod shock tube,
# with high pressure gas in the left half of the duct
# and low pressure gas in the right half.
for i in range(n):
    cells[i].x = 0.5*(faces[i].x + faces[i+1].x)
    cells[i].vol = 0.5*(faces[i].area + faces[i+1].area) * dx
    gas = cells[i].gas
    gas.p, gas.T = (1.0e5, 348.4) if i < n/2 else (1.0e4, 278.8)
    gas.update_thermo_from_pT()
    cells[i].encode_conserved()

for i in range(n):
    print(i, cells[i])

print("Start time stepping.")
# Leave the cells at the ends as ghost-cells, with fixed properties,
# and update all of the interior cells over a small time step.
dt = 1.0e-6 # seconds
n_steps = 10
gsLstar = GasState(gmodel)
gsRstar = GasState(gmodel)
for j in range(n_steps):
    # Phase 1: compute flow states at interior faces.
    for i in range(1,n):
        pstar, wstar, wL, wR, faces[i].vel = \
            flow.osher_riemann(cells[i-1].gas, cells[i].gas,
                               cells[i-1].vel, cells[i].vel,
                               gsLstar, gsRstar, faces[i].gas)
    # Phase 2: Update interior cells.
    for i in range(1,n-1):
        ci = cells[i]; fi = faces[i]; fip1 = faces[i+1]
        vol = ci.vol; areai = fi.area; areaip1 = fip1.area
        rhoi = fi.gas.rho; rhoip1 = fip1.gas.rho
        veli = fi.vel; velip1 = fip1.vel
        mass_dt = (rhoi*veli*areai - rhoip1*velip1*areaip1)/vol
        pi = fi.gas.p; pip1 = fip1.gas.p
        momentum_dt = (pi*areai - pip1*areaip1 + ci.gas.p*(areaip1-areai) +
                       rhoi*(veli**2)*areai - rhoip1*(velip1**2)*areaip1)/vol
        hi = pi/rhoi + fi.gas.u + 0.5*(veli**2)
        hip1 = pip1/rhoip1 + fip1.gas.u + 0.5*(velip1**2)
        energy_dt = (rhoi*hi*areai - rhoip1*hip1*areaip1)/vol
        ci.mass += dt * mass_dt
        ci.momentum += dt * momentum_dt
        ci.energy += dt * energy_dt
        ci.decode_conserved()
    # end of time step

for i in range(n):
    print(i, cells[i])
   
