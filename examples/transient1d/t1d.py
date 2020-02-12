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
n = 200
cells = [Cell(gmodel) for i in range(n)]
faces = [Face(gmodel) for i in range(n+1)]

print("Quasi-one-dimensional duct is defined by its area at position x.")
def Area(x): return 1.0 # metres^2
dx = 1.0/n # metres
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

def write_flow_data(fp, t, cdata):
    "Write flow data for this instant in a format suitable for gnuplot."
    for i in range(len(cdata)):
        ci = cdata[i]
        fp.write("%d\t%g\t%g\t%g\t%g\t%g\t%g\n" %
                 (i, ci.x, t, ci.gas.p, ci.gas.T, ci.gas.rho, ci.vel))
    fp.write("\n")
    return

t = 0.0 # time, in seconds
fp = open('xt.data', 'w')
write_flow_data(fp, t, cells)

print("Start time stepping.")
# Leave the cells at the ends as ghost-cells, with fixed properties,
# and update all of the interior cells over a small time step.
dt = 0.5e-6 # time step, in seconds
n_steps = 1200
plot_steps = 20
gsLstar = GasState(gmodel)
gsRstar = GasState(gmodel)
for j in range(1, n_steps+1):
    # Phase 1: compute flow states at interior faces.
    for i in range(1,n):
        pstar, wstar, wL, wR, faces[i].vel = \
            flow.osher_riemann(cells[i-1].gas, cells[i].gas,
                               cells[i-1].vel, cells[i].vel,
                               gsLstar, gsRstar, faces[i].gas)
    # Phase 2: Update interior cells.
    for i in range(1,n-1):
        # Local names for quantities.
        ci = cells[i]; fi = faces[i]; fip1 = faces[i+1]
        vol = ci.vol; Ai = fi.area; Aip1 = fip1.area
        rhoi = fi.gas.rho; rhoip1 = fip1.gas.rho
        vi = fi.vel; vip1 = fip1.vel
        pi = fi.gas.p; pip1 = fip1.gas.p
        hi = pi/rhoi + fi.gas.u + 0.5*(vi**2)
        hip1 = pip1/rhoip1 + fip1.gas.u + 0.5*(vip1**2)
        # Mass fluxes at the faces.
        mfluxi = rhoi*vi*Ai
        mfluxip1 = rhoip1*vip1*Aip1
        # Rate of accumulation within the cell.
        mass_dt = (mfluxi - mfluxip1)/vol
        momentum_dt = (pi*Ai - pip1*Aip1 + ci.gas.p*(Aip1-Ai) +
                       mfluxi*vi - mfluxip1*vip1)/vol
        energy_dt = (mfluxi*hi - mfluxip1*hip1)/vol
        # Euler update for conserved quantities.
        ci.mass += dt * mass_dt
        ci.momentum += dt * momentum_dt
        ci.energy += dt * energy_dt
        ci.decode_conserved()
    t += dt
    if (j % plot_steps) == 0: write_flow_data(fp, t, cells)
    print("end of time step %d, t=%g seconds" % (j, t))
fp.close()
fp = open('final.data', 'w')
write_flow_data(fp, t, cells)
fp.close()
   
