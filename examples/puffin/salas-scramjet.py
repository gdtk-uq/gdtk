# Salas' scramjet example from Figure 11.
# PJ 2022-02-03

init_gas_model("ideal-air-gas-model.lua")
gas1 = GasState(config.gmodel)
gas1.p = 100.0e3
gas1.T = 300.0
gas1.update_thermo_from_pT()
gas1.update_sound_speed()
M1 = 3.0; V1 = M1 * gas1.a
print("V1=", V1)

from math import sin, cos, radians
V1x = V1 * cos(radians(8.0))
V1y = V1 * sin(radians(8.0))

config.max_step_relax = 40

from eilmer.geom.xpath import XPath
path_0 = XPath().moveto(0,0).lineto(1,0).lineto(2,0)
path_0.lineto(3,0).lineto(4,0.02).lineto(5,0.1)
path_0.lineto(6,0.244).lineto(7,0.425).lineto(8,0.625).lineto(9,0.825)
path_1a = XPath().moveto(0,1).lineto(3,0.75).lineto(5,1).lineto(9,1.4)
path_1b = XPath().moveto(0,1).lineto(3,1.25).lineto(5,1).lineto(9,1.4)
path_2a = XPath().moveto(0,2).lineto(0.5,2).lineto(2.5,1.8).lineto(4,2).lineto(9,2)
path_2b = XPath().moveto(0,2).lineto(0.5,2).lineto(2.5,2.2).lineto(4,2).lineto(9,2)
path_3a = XPath().moveto(0,3).lineto(3,2.75).lineto(5,3).lineto(9,2.6)
path_3b = XPath().moveto(0,3).lineto(3,3.25).lineto(5,3).lineto(9,2.6)
path_4 = XPath().moveto(0,4).lineto(1,4).lineto(2,4)
path_4.lineto(3,4).lineto(4,4-0.02).lineto(5,4-0.1)
path_4.lineto(6,4-0.244).lineto(7,4-0.425).lineto(8,4-0.625).lineto(9,4-0.825)
def bc_0(x): return 0
def bc_1(x):
    if x < 5.0: return 0
    return 1
def bc_2(x):
    if x < 0.5: return 1
    if x < 4.0: return 0
    return 1
def bc_3(x):
    if x < 5.0: return 0
    return 1
def bc_4(x): return 0

config.max_x = path_0.xs[-1]
config.dx = config.max_x/1000

st0 = StreamTube(gas=gas1, velx=V1x, vely=V1y,
                 y0=path_0, y1=path_1a, bc0=bc_0, bc1=bc_1, ncells=60)
st1 = StreamTube(gas=gas1, velx=V1x, vely=V1y,
                 y0=path_1b, y1=path_2a, bc0=bc_1, bc1=bc_2, ncells=60)
st2 = StreamTube(gas=gas1, velx=V1x, vely=V1y,
                 y0=path_2b, y1=path_3a, bc0=bc_2, bc1=bc_3, ncells=60)
st3 = StreamTube(gas=gas1, velx=V1x, vely=V1y,
                 y0=path_3b, y1=path_4, bc0=bc_3, bc1=bc_4, ncells=60)
