# kernel.py
"""
Kernel data and functions for the Isentropic Method-Of-Characteristics.
This is a regrowth of the old IMOC code that was implemented in C+Tcl/Tk.

Author(s):
  Peter J.
  Centre for Hypersonics,
  School of Mechanical Engineering, U of Q.

Version:
  2019-12-28: Let's start coding in Python and see how it develops...
"""

import numpy as np
from scipy import spatial

debugLevel = 0

g = 1.4  # Ratio of specific heats for gas
p0 = 1.0 # Nondimensional total pressure for flow
T0 = 1.0 # Nondimensional total temperature for flow
axisymmetric = False

nodes = [] # Storage for the Node objects.
mesh_indices = [] # Indices of the nodes that have been added to the characteristics mesh.

walls = [] # Storage for the user-defined functions that specify the walls.

class Node(object):
    __slots__ = ['indx', 'x', 'y',
                 'theta', 'nu', 'mach',
                 'cplus_up', 'cplus_down',
                 'cminus_up', 'cminus_down',
                 'czero_up', 'czero_down']

    def __init__(self, indx=None, x=0.0, y=0.0,
                 theta=0.0, nu=0.0, mach=0.0,
                 cplus_up=-1, cplus_down=-1,
                 cminus_up=-1, cminus_down=-1,
                 czero_up=-1, czero_down=-1):
        """
        Initialize a Node.
        The default state will be empty and unconnected.

        We store the new Node in the nodes list at indx.
        If indx is not supplied, the new Node is appended to nodes.

        Note that we depend upon the index within the nodes storage list
        being the same as the values of the indx attribute.
        """
        if indx is None: indx = len(nodes)
        self.indx = indx
        self.x = x; self.y = y
        self.theta = theta; self.nu = nu; self.mach = mach
        self.cplus_up = cplus_up; self.cplus_down = cplus_down
        self.cminus_up = cminus_up; self.cminus_down = cminus_down
        self.czero_up = czero_up; self.czero_down = czero_down
        if indx < len(nodes):
            nodes[indx] = self
        else:
            nodes.append(self)
        return

    def __repr__(self):
        """
        A string representation of a Node.
        """
        strng = "Node(indx=%d, x=%g, y=%g" % (self.indx, self.x, self.y)
        strng += ", theta=%g, nu=%g, mach=%g" % (self.theta, self.nu, self.mach)
        strng += ", cplus_up=%d, cplus_down=%d" % (self.cplus_up, self.cplus_down)
        strng += ", cminus_up=%d, cminus_down=%d" % (self.cminus_up, self.cminus_down)
        strng += ", czero_up=%d, czero_down=%d)" % (self.czero_up, self.czero_down)
        return strng

def create_kd_tree():
    kdtree = spatial.KDTree(np.array([(node.x, node.y) for node in nodes]))
    return kdtree

def find_nodes_near(x, y, tol=0.0, max_count=30, kdtree=None):
    """
    An attempt to replicate PJs C code to find the nodes near a given
    x, y position. If tol is <= 0.0 then only a single node is returned,
    otherwise an array of nodes, up to max_count long, will be returned.
    INPUTS:
            x, y      - position of interest
            tol       - radius of interest
            max_count - maximum number of nodes to collect
            kdtree    - if a kdtree object has been created this can be
                        used for much, much faster searching
    OUTPUT:
            idx_near  - a list of node indices nearby
    """
    idx_near = []
    if kdtree is None:
        if tol <= 0.0:
            # Find the nearest node
            idx_near.append(-1)
            dist_near = 1.0e6 # Arbitrarily large
            for idx, node in enumerate(nodes):
                node_dist = np.sqrt((x - node.x)**2 + (y - node.y)**2)
                if node_dist < dist_near:
                    dist_near = node_dist
                    idx_near[-1] = idx
        else:
            # Collect an array of the closest nodes
            # NOTE: There doesn't seem to be any mechanisms of ensuring they are
            # the closest nodes, just the first ones within the radius of interest
            for idx, node in enumerate(nodes):
                node_dist = np.sqrt((x - node.x)**2 + (y - node.y)**2)
                if node_dist < tol:
                    idx_near.append(idx)
                if len(idx_near) >= max_count: break
    else:
        _, pnts = kdtree.query((x, y), max_count, distance_upper_bound=tol)
        if tol <= 0.0:
            idx_near = [pnts[0]]
        else:
            # Query seems to do something I don't like whereby the list is filled
            # to the length of the max_count with the value of the length of nodes
            # Going to manually remove any values which equal the number of nodes
            # (this means the final node cannot be used)
            idx_near = pnts
            idx_near = np.delete(idx_near, np.where(idx_near == len(nodes)))
            idx_near = np.unique(idx_near)
    return idx_near


class Wall(object):
    __slots__ = ['f', 'dfdx', 'x_min', 'x_max']

    def __init__(self, f, x_min, x_max, dfdx=None, dx=1.0e-6):
        """
        Accept a user-supplied function f(x) to define the wall.
        x_min, x_max: range of x that is of interest (for plotting).
        Optionally, accept a function to define the slope of the wall.

        If that function is not supplied,
        we construct a finite-difference approximation.
        """
        if not callable(f):
            raise RuntimeError("Expected a callable function for f.")
        else:
            self.f = f
        if not dfdx:
            # Construct a finite-difference approximation.
            self.dfdx = lambda x, f=f, dx=dx: (f(x+dx)-f(x))/dx
        elif callable(dfdx):
            self.dfdx = dfdx
        else:
            raise RuntimeError("Expected a callable function for optional dfdx.")
        self.x_min = x_min
        self.x_max = x_max

    def __call__(self, x):
        return self.f(x)
