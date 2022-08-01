# anderson_11p1.py
# Script for Exercise 11.1 in J.D.Anderson Modern Compressible Flow
#
# Usage:
# python3 anderson_11p1.py
#
# PJ 2022-01-11
#
import math
import gdtk.imoc.kernel as kernel
import gdtk.imoc.unit_process as unit
from gdtk.ideal_gas_flow import PM1, PM2

print("Exercise 11.1: Minimum-length 2D nozzle")
kernel.axisymmetric = False
kernel.g = 1.4

# Put in a wall along the x-axis.
def f0(x): return 0.0
wall0 = kernel.Wall(f0, 0.0, 1.0)

print("Create a centered expansion fan at (0.0 0.1)")
fan = []
for i in range(7):
    nu_value = math.radians(0.375 + i*3.0)
    new_node = kernel.Node(x=0.0, y=0.1, nu=nu_value, mach=PM2(nu_value), theta=nu_value)
    kernel.register_node_in_mesh(new_node)
    fan.append(new_node.indx)

print("Compute the first wall node radiating from the fan.")
old_first = unit.cminus_wall(wall0, fan[0])

print("Compute the rest of the fan radiating down to the wall.")
for i in range(1,7):
    new_nodes = unit.march_along_cminus(old_first, fan[i], 'down')
    axis_node = unit.cminus_wall(wall0, new_nodes[-1])
    # Note that for the next line we start marching with the second node on the recent line.
    old_first = new_nodes[1]

node = kernel.nodes[axis_node]
x_cone = node.x; M_cone = node.mach
print(f"Mach cone starts at x={x_cone} with M={M_cone}")

print("Put down a number of nodes along the x-axis with constant M.")
print("Work back upstream from each of these nodes.")
dL = 0.05
mach_angle = math.asin(1.0/M_cone)
dx = dL * math.cos(mach_angle); dy = dL * math.sin(mach_angle)

# When starting to march upstream on the last Cminus line of the fan,
# we want to start with the last interior point that was generated
# while marching down that characteristic line.
old = new_nodes[-1]

# Starting at final axis_node of the fan, we want to construct the leading edge
# of the uniform region going downstream.
old_edge = axis_node
for i in range(1,11):
    new_edge_node = kernel.Node(x=x_cone+i*dx, y=i*dy, nu=PM1(M_cone), mach=M_cone, theta=0.0)
    kernel.register_node_in_mesh(new_edge_node)
    new_edge_node.cplus_up = old_edge
    kernel.nodes[old_edge].cplus_down = new_edge_node.indx
    new_nodes = unit.march_along_cminus(old, new_edge_node.indx, 'up')
    old = new_nodes[1] # The node just off the axis.
    old_edge = new_edge_node.indx

print("Start at the last node on the fan and step along a streamline.")
# We continue until either
# (1) the interpolation fails, or
# (2) we cross the characteristic defining the start of the uniform flow region.
indx = fan[-1]
kernel.register_streamline_start(indx)
while True:
    indx = unit.step_stream_node(indx, 0.05)
    if indx is None:
        break
    else:
        node = kernel.nodes[indx]
        y_cone = (node.x - x_cone)*dy/dx
        if node.y <= y_cone: break
