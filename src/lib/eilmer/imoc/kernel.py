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

debugLevel = 0

g = 1.4  # Ratio of specific heats for gas
p0 = 1.0 # Nondimensional total pressure for flow
T0 = 1.0 # Nondimensional total temperature for flow
axisymmetric = False

nodes = [] # Storage for the Node objects.
mesh_nodes = [] # Nodes that have been added to the characteristics mesh.

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
