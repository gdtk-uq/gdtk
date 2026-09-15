# kernel.py
"""
Kernel data and functions for the Isentropic Method-Of-Characteristics.
This is a regrowth of the old IMOC code that was implemented in C+Tcl/Tk.

Author(s):
  Peter J. (1) and Fabian Zander (2)
  1. Centre for Hypersonics, School of Mechanical Engineering, UofQ.
  2. Institute for Advanced Engineering and Space Sciences, USQ

Version:
  2019-12-28: Let's start coding in Python and see how it develops...
"""

import sys
import numpy as np
from scipy import spatial

debugLevel = 0

g = 1.4  # Ratio of specific heats for gas
p0 = 1.0 # Nondimensional total pressure for flow
T0 = 1.0 # Nondimensional total temperature for flow
axisymmetric = False

# Storage for all the Node objects.
nodes = []
# The indices that are used to link nodes into a characteristic mesh
# and into streamlines are the locations of the nodes in this list.
# Once a node is constructed and added to this list, it is important
# not to move it, else our simple linking will not work.
#
# A node may, however, be removed with delete_node().  Its slot is retained
# in the nodes list, so that every surviving node keeps its index, and the
# index is recorded in the deleted set below.  This mirrors the array of
# node pointers in the original moc C code, where deleting a node left a
# NULL entry behind that a later CreateNode could fill again.

# Indices of the slots in nodes that no longer hold a live node.
deleted = set()

# Indices of the nodes that have been added to the characteristics mesh.
char_mesh = []
# Indices of the starting nodes on the streamlines.
streamlines = []
# Storage for the user-defined functions that specify the walls.
walls = []

class Node(object):
    __slots__ = ['indx', 'x', 'y',
                 'theta', 'nu', 'mach',
                 'cplus_up', 'cplus_down',
                 'cminus_up', 'cminus_down',
                 'czero_up', 'czero_down']

    def __init__(self, indx=None, x=0.0, y=0.0,
                 theta=0.0, nu=0.0, mach=0.0,
                 cplus_up=None, cplus_down=None,
                 cminus_up=None, cminus_down=None,
                 czero_up=None, czero_down=None):
        """
        Initialize a Node.
        The default state will be empty and unconnected.

        We store the new Node in the nodes list at indx.
        If indx is not supplied, the new Node is appended to nodes.

        Note that we depend upon the index within the nodes storage list
        being the same as the value of the indx attribute.
        """
        if indx is None: indx = len(nodes)
        self.indx = indx
        self.x = x; self.y = y
        self.theta = theta; self.nu = nu; self.mach = mach
        # Indices for linking nodes.
        # A value of None indicates no link.
        self.cplus_up = cplus_up; self.cplus_down = cplus_down
        self.cminus_up = cminus_up; self.cminus_down = cminus_down
        self.czero_up = czero_up; self.czero_down = czero_down
        if indx < len(nodes):
            nodes[indx] = self
        elif indx == len(nodes):
            nodes.append(self)
        else:
            # Appending here would store the Node at a position that does not
            # match its indx attribute and quietly corrupt every link into it.
            raise RuntimeError("Node indx=%d is beyond the end of the nodes list"
                               " (len=%d); use create_node() to reach it."
                               % (indx, len(nodes)))
        return

    def __repr__(self):
        """
        A string representation of a Node.
        """
        strng = "Node(indx=%d, x=%g, y=%g" % (self.indx, self.x, self.y)
        strng += ", theta=%g, nu=%g, mach=%g" % (self.theta, self.nu, self.mach)
        strng += ", cplus_up=%s, cplus_down=%s" % (self.cplus_up, self.cplus_down)
        strng += ", cminus_up=%s, cminus_down=%s" % (self.cminus_up, self.cminus_down)
        strng += ", czero_up=%s, czero_down=%s)" % (self.czero_up, self.czero_down)
        return strng


def _as_index(i):
    """
    Accept either a node index or a Node and return the index.
    """
    if isinstance(i, Node): return i.indx
    if isinstance(i, (int, np.integer)): return int(i)
    raise RuntimeError("Not an index nor a Node: %s" % repr(i))


def valid_node(i):
    """
    Returns True if i refers to a slot that holds a live node.

    Port of ValidNode from the original moc C code.
    """
    if isinstance(i, Node): i = i.indx
    if not isinstance(i, (int, np.integer)): return False
    return 0 <= i < len(nodes) and i not in deleted


def number_of_nodes():
    """
    Returns the number of live nodes.

    Port of GetNumberOfNodes from the original moc C code.  Note that this is
    the count of live nodes and not the length of the nodes list; the two
    differ once any node has been deleted.
    """
    return len(nodes) - len(deleted)


def get_next_node_id(id_start=-1):
    """
    Search for the next live node at an index greater than id_start.
    Send -1 to find the first one.  Returns -1 if there is no such node.

    Port of GetNextNodeId from the original moc C code.  Useful for walking a
    nodes list that has holes in it.
    """
    if id_start < -1: id_start = -1
    for i in range(id_start+1, len(nodes)):
        if i not in deleted: return i
    return -1


def create_node(indx=-1):
    """
    Make a new, empty node and return its index.

    indx <  0 : use the lowest-numbered free slot, appending if there is none.
    indx >= 0 : use that particular slot, deleting whatever was there.

    Port of CreateNode from the original moc C code.  Reusing the lowest free
    slot is not just tidiness; scripts written against the original kernel can
    depend on the particular indices that come back, so we reproduce its
    choice exactly.
    """
    if indx is None: indx = -1
    indx = int(indx)
    if indx < 0:
        if deleted:
            # Take the first available space, as the C kernel does.
            i = min(deleted)
            deleted.discard(i)
            return Node(indx=i).indx
        return Node().indx
    # A particular slot has been requested.
    while len(nodes) <= indx:
        # Pad with dead slots so that nodes[indx] becomes addressable.
        # The C kernel had a preallocated array of node pointers, so any index
        # below MAX_NODES could be filled directly.
        deleted.add(Node().indx)
    if indx not in deleted:
        # If something is already there, destroy it.
        delete_node(indx)
    deleted.discard(indx)
    return Node(indx=indx).indx


def delete_node(i):
    """
    Remove nodes[i] from the mesh and leave its slot free for reuse.

    Port of DeleteNode from the original moc C code.  As in the C code, the
    characteristic and streamline connections of the remaining nodes are
    retained: the node's up-neighbour and down-neighbour are joined to each
    other, so that a line running through the deleted node stays walkable.
    """
    i = _as_index(i)
    if not valid_node(i): return False
    n = nodes[i]
    czero_down = n.czero_down
    for up, down in (('cplus_up', 'cplus_down'),
                     ('cminus_up', 'cminus_down'),
                     ('czero_up', 'czero_down')):
        u = getattr(n, up); d = getattr(n, down)
        # Guard each neighbour, in case a link into a slot that has already
        # been deleted has been left behind somewhere.
        if valid_node(u) and getattr(nodes[u], down) == i:
            setattr(nodes[u], down, d)
        if valid_node(d) and getattr(nodes[d], up) == i:
            setattr(nodes[d], up, u)
        setattr(n, up, None); setattr(n, down, None)
    while i in char_mesh: char_mesh.remove(i)
    if i in streamlines:
        # Keep the streamline registered, starting from the next node along.
        streamlines.remove(i)
        if czero_down is not None and valid_node(czero_down):
            streamlines.append(czero_down)
    deleted.add(i)
    return True


def delete_all_nodes():
    """
    Discard every node, starting afresh.

    The walls are left alone, unlike the DeleteAll procedure of the Tcl layer.
    """
    nodes.clear(); deleted.clear()
    char_mesh.clear(); streamlines.clear()
    return


def create_kd_tree():
    """
    Build a KDTree over the node positions, for fast searching.

    The tree holds a point for every slot in the nodes list, including the
    deleted ones, so that the indices it reports are node indices.
    find_nodes_near() filters the deleted slots out of the search results.
    """
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
            idx_near  - a list of node indices for nearby nodes
    """
    idx_near = []
    if kdtree is None:
        # Just do a slow search.
        if tol <= 0.0:
            # Find the nearest node.
            idx_near.append(-1) # This single value should be overwritten.
            dist_near = sys.float_info.max
            for idx, node in enumerate(nodes):
                if idx in deleted: continue
                dist = np.sqrt((x - node.x)**2 + (y - node.y)**2)
                if dist < dist_near:
                    dist_near = dist
                    idx_near[0] = idx
            if idx_near[0] < 0:
                raise RuntimeError("No nearest node found.")
        else:
            # Collect an array of the closest nodes as tuples of distance and index.
            close_nodes = []
            for idx, node in enumerate(nodes):
                if idx in deleted: continue
                dist = np.sqrt((x - node.x)**2 + (y - node.y)**2)
                if dist < tol: close_nodes.append((dist,idx))
            close_nodes.sort(key=lambda t: t[0]) # Sort on distance.
            idx_near = [t[1] for t in close_nodes] # Keep just the indices.
            if len(idx_near) > max_count: idx_near = idx_near[0:max_count]
    else:
        # Use the kdtree to do a fast search.
        # The tree carries a point for every slot in the nodes list, deleted
        # slots included, so ask for enough neighbours that the deleted ones
        # cannot crowd the live ones out, then drop them from the result.
        n_query = max_count
        if deleted: n_query = min(max_count + len(deleted), kdtree.n)
        if tol <= 0.0:
            # We want the nearest node, so do not bound the search distance.
            # Passing distance_upper_bound=0.0 finds nothing at all and leaves
            # the caller holding kdtree.n, which is not a valid node index.
            _, pnts = kdtree.query((x, y), n_query)
        else:
            _, pnts = kdtree.query((x, y), n_query, distance_upper_bound=tol)
        pnts = np.atleast_1d(pnts)
        if tol <= 0.0:
            idx_near = [i for i in pnts
                        if i != kdtree.n and i not in deleted][0:1]
        else:
            # Query seems to do something I don't like whereby the list is filled
            # to the length of the max_count with the value of the length of nodes
            # Going to manually remove any values which equal the number of nodes
            # (this means the final node cannot be used)
            idx_near = pnts
            idx_near = np.delete(idx_near, np.where(idx_near == kdtree.n))
            idx_near = np.unique(idx_near)
            if deleted:
                idx_near = idx_near[~np.isin(idx_near, list(deleted))]
                if len(idx_near) > max_count: idx_near = idx_near[0:max_count]
    return idx_near

def register_node_in_mesh(i):
    """
    Add nodes[i] to the characteristic mesh.
    """
    if isinstance(i, int):
        char_mesh.append(i)
    elif isinstance(i, Node):
        char_mesh.append(i.indx)
    else:
        raise RuntimeError("Not an index nor a Node")
    return len(char_mesh)

def register_streamline_start(i):
    """
    Register nodes[i] as the starting node for a streamline.
    """
    if isinstance(i, int):
        streamlines.append(i)
    elif isinstance(i, Node):
        streamlines.append(i.indx)
    else:
        raise RuntimeError("Not an index nor a Node")
    return len(streamlines)

def get_streamline_nodes(i):
    """
    Returns the indices of the streamline nodes,
    for the streamline going through nodes[i].
    """
    indices = [i,]
    j = nodes[i].czero_down
    while j is not None:
        indices.append(j)
        j = nodes[j].czero_down
    j = nodes[i].czero_up
    while j is not None:
        indices.insert(0, j)
        j = nodes[j].czero_up
    return indices


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
