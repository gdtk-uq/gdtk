# TUSQ-Mach4-MOC.py
# A TUSQ Mach 4 nozzle design for the RDE work
#
# The isentropic nozzle design is undertaken using the method of Sivell
# ("Aerodynamic Design of Axisymmetric Hypersonic Wind-Tunnel Nozzles",
# James C. Sivells, J. Spacecraft, Vol. 7, No. 11, November, 1970).
#
# This method calculates an axial Mach number distribution along the nozzle
# which forms the basis for the MOC calculation. This is a little bit
# different to other common uses of MOC as it starts from the axis and works
# 'up' to the nozzle wall (and beyond it). A sample solution for the centreline
# Mach distribution is provided in the dat file. The other important number is
# the node number at which the design Mach has been reached. In this case it is
# node 27.
#
# The inputs to the analytical calculation of the centreline profile can be used
# as a check for the final result. The centreline profile was calculated for a
# Mach 4 nozzle with a throat radius of 20mm, an exit radius of 65mm and has a
# computed length of 480mm.
#
# Fabian Zander
# Tuesday 14th April, 2020
# Covid Isolation, Boshammer St, Rangeville
#
import gdtk.imoc.kernel as kernel
import gdtk.imoc.unit_process as unit
from gdtk.ideal_gas_flow import PM1
import numpy as np
import matplotlib.pyplot as mplt
# Import the centreline data, two columns with x and Mach
centreline_data = np.loadtxt('machRDEnozzleMachProfile.dat')
cl_x = centreline_data[:,0]
cl_mach = centreline_data[:,1]
# Set up the calculation parameters
kernel.axisymmetric = True
kernel.g = 1.4
# First create and add all the centreline nodes
for i in range(len(cl_x)):
    n = kernel.Node(x=cl_x[i], y=0.0, nu=PM1(cl_mach[i]), mach=cl_mach[i], theta=0.0)
    kernel.register_node_in_mesh(n)
# Now add the interior nodes. As we are working up from the centreline at each
# step up we lose one node (two nodes 'merge' into one). The max height of
# 150mm was iteratively chosen to ensure the computed mesh covered the entire
# region of interest
top_node_count = len(cl_mach)
while kernel.nodes[-1].y < 0.15 and top_node_count >= 2:
    for i in range(len(kernel.nodes) - top_node_count, len(kernel.nodes) - 1):
        unit.interior(i+1, i)
    top_node_count -= 1
# Now, we know our Mach cone starts at node 27, we use this to calculate the
# node which represents our nozzle exit. From this we can then use streamlines
# to define our wall profile
nozz_exit_node = 27
# We are going to find the first node above the design radius of 65mm
while kernel.nodes[nozz_exit_node].y < 0.065:
    nozz_exit_node = kernel.nodes[nozz_exit_node].cplus_down
# We can do a bit of a check here by looking at the Mach number and theta of
# the computed node, these should be Mach 4.0 and 0.0
print(f"The nozzle exit node (node {nozz_exit_node}) has a Mach number of "
      + f"{kernel.nodes[nozz_exit_node].mach} and a theta of "
      + f"{kernel.nodes[nozz_exit_node].theta:.5f}")
# Now, we are going to create our nozzle wall by following the streamline
# upstream from the nozzle exit point, until we reach the throat.
nozz_wall_nodes = [nozz_exit_node] # Start at nozzle exit node.
dL = -8.0e-3 # Negative increment steps upstream.
kdt = kernel.create_kd_tree() # Use kdtree for a fast search.
while kernel.nodes[nozz_wall_nodes[-1]].x > 0.01:
    new_idx = unit.step_stream_node(nozz_wall_nodes[-1], dL, kdtree=kdt)
    if new_idx is None: break
    nozz_wall_nodes.append(new_idx)
nozz_wall_nodes.reverse()
kernel.register_streamline_start(nozz_wall_nodes[0])
# Now, we can do a check of the throat radius and see that this matches what
# we are expecting from our analytical solution, it should be ~20mm.
# The value will probably be slightly higher due to the way we defined the end
# of the nozzle wall solution above
print(f"The throat radius is approximately "
      + f"{kernel.nodes[nozz_wall_nodes[0]].y*1000:.2f} mm.")

# We just want the walls to appear as distinct items in the GUI.
# The axis is the lower wall.
def f0(x): return 0.0
# And the create the upper (nozzle) wall as a closure around
# the arrays of xs and ys.
xs = np.array([kernel.nodes[i].x for i in nozz_wall_nodes])
ys = np.array([kernel.nodes[i].y for i in nozz_wall_nodes])
def f1(x, wall_x=xs, wall_y=ys):
    import numpy as np
    return np.interp(x, wall_x, wall_y)
kernel.walls = [kernel.Wall(f0, 0.0, 0.5), kernel.Wall(f1, 0.0, 0.5)]
