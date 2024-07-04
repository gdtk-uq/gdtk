"""
A Puffin input file for simulating a Busemann diffuser.

This example is taken from Section 5 of Moelder (2003).
It makes use of the busemann.py module in the gdtk package
to generate the Busemann contour.

Reference:

Moelder, S. (2003)
A Benchmark for Internal Flow CFD Codes
Computational Fluid Dynamics Journal, 12(2):47, pp.408--414

.. Author: RJG

.. Version: 2024-07-04
"""
from gdtk.busemann import BusemannDiffuser
from gdtk.numeric.spline import CubicSpline
from gdtk.geom.xpath import XPath
from gdtk.geom.path import Polyline
from math import asin

config.axisymmetric = True

init_gas_model("ideal-air-gas-model.lua")
gas1 = GasState(config.gmodel)
gamma = gas1.gamma
q1 = 50e3 # Pa
M1 = 5.76788
gas1.p = 2*q1/((gamma - 1)*M1*M1)
gas1.T = 250.0
gas1.update_thermo_from_pT()
gas1.update_sound_speed()
V1 = M1 * gas1.a

config.max_step_relax = 40

# Set up Busemann diffuser contour
M2 = 3.0
k = 1.4
theta_23 = asin(k/M2)
bd = BusemannDiffuser(M2, theta_23)
npts = 50
bd.generate_contour(1.0, 0.001)
contour = bd.contour_as_spline(npts)
thrt_start = contour.x[-1]
thrt_height = contour(thrt_start)
thrt_end = thrt_start + 0.2
L = thrt_end - contour.x[0]
# need to shift surface so that it starts at x = 0
d_start = contour.x[0]
def upper_y(x):
    xdash = x + d_start
    if xdash < thrt_start:
        return contour(xdash)
    else:
        return thrt_height
    
def upper_bc(x): return 0

def lower_y(x): return 0.0
def lower_bc(x): return 0

config.max_x = thrt_end - d_start
config.dx = L/1000

st1 = StreamTube(gas=gas1, velx=V1, vely=0.0,
                 y0=lower_y, y1=upper_y,
                 bc0=lower_bc, bc1=upper_bc,
                 ncells=80)









