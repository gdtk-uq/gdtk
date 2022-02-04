# Salas' axisymmetric nozzle example.
# PJ 2022-02-03

init_gas_model("ideal-air-gas-model.lua")
gas1 = GasState(config.gmodel)
gas1.p = 100.0e3
gas1.T = 300.0
gas1.update_thermo_from_pT()
gas1.update_sound_speed()
M1 = 2.0
V1 = M1 * gas1.a
print("V1=", V1)

config.max_step_relax = 40

import csv
xs = []; ys = []
with open('salas-nozzle-path.tsv', 'r') as tsvf:
    tsv_data = csv.reader(tsvf, delimiter='\t')
    for row in tsv_data:
        if row[0] == 'x': continue
        xs.append(float(row[0]))
        ys.append(float(row[1]))

from eilmer.spline import CubicSpline
upper_y = CubicSpline(xs, ys)
def lower_y(x): return 0.0
def lower_bc(x): return 0
def upper_bc(x): return 0

config.max_x = xs[-1]
config.dx = config.max_x/1000

st1 = StreamTube(gas=gas1, velx=V1, vely=0.0,
                 y0=lower_y, y1=upper_y, bc0=lower_bc, bc1=upper_bc,
                 ncells=80)
