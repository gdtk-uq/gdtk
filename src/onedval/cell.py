# Author: Rowan J. Gollan
# Place: UQ, Brisbane, Queensland, Australia
# Date: 25-Jun-2012

from cfpylib.geom.minimal_geometry import *

def tri_centroid(p0, p1, p2):
    """
    Compute centroid of triangle given 3 points
    
           p1
           +
          / \
         /   \
        /     \
       +-------+
       p0      p2
    """
    return (1.0/3.0)*(p0 + p1 + p2)

def tri_area(p0, p1, p2):
    """
    Compute the area of a triangle.
    """
    vector_area = 0.5*cross(p1-p0, p2-p1)
    return abs(vector_area)

def tri_normal(p0, p1, p2):
    """
    Compute the normal direction of a triangle in a plane.
    """
    vector_area = 0.5*cross(p1-p0, p2-p1)
    return unit(vector_area)

class Cell(object):
    """
    Python base class to abstract the behaviour of
    all cells.
    """
    def __init__(self, data):
        self._data = data.copy()
        return

    def variables(self):
        return self._data.keys()

    def get(self, k, default='no value'):
        return self._data.get(k, default)

    def area(self):
        return self._area

    def normal(self):
        return self._normal

    def centroid(self):
        return self._centroid

class QuadCell(Cell):
    """
    A quadrilateral cell in a 2D plane.
    """
    count = 0
    def __init__(self, data, pts):
        Cell.__init__(self, data)
        self._centroid = quad_centroid(pts[0], pts[1], pts[2], pts[3])
        self._normal = unit(quad_normal(pts[0], pts[1], pts[2], pts[3]))
        self._area = quad_area(pts[0], pts[1], pts[2], pts[3])
        QuadCell.count = QuadCell.count + 1
        return

class TriCell(Cell):
    """
    A triangular cell in a 2D plane.
    """
    count = 0
    def __init__(self, data, pts):
        Cell.__init__(self, data)
        self._centroid = tri_centroid(pts[0], pts[1], pts[2])
        self._normal = tri_normal(pts[0], pts[1], pts[2])
        self._area = tri_area(pts[0], pts[1], pts[2])
        TriCell.count = TriCell.count + 1
        return

class AxiCell(Cell):
    """
    An axisymmetric cell defined by inner and outer radii.
    """
    count = 0
    def __init__(self, data, pts):
        Cell.__init__(self, data)
        self._centroid = 0.5*(pts[0] + pts[1])
        self._normal = Vector(1.0, 0.0, 0.0)
        self._area = self._centroid.y*(pts[1].y - pts[0].y)
        AxiCell.count = AxiCell.count + 1
        return

def unique(lst):
    seen = set()
    unique = []
    for item in lst:
        if not item in seen:
            unique.append(item)
            seen.add(item)
    return unique

def create_cell(idx, data, cell_cnrs, var_map, scale):
    c_list = unique(cell_cnrs[idx])
    pts = []
    xlabel = var_map['x']
    ylabel = var_map['y']
    zlabel = var_map['z']
    for cnr in c_list:
        i = cnr-1
        pts.append(scale*Vector(data[xlabel][i], data[ylabel][i], data[zlabel][i]))

    if len(pts) == 3:
        pC = tri_centroid(pts[0], pts[1], pts[2])
        area = tri_area(pts[0], pts[1], pts[2])
        a0 = tri_area(pts[1], pts[2], pC)
        a1 = tri_area(pts[0], pts[2], pC)
        a2 = tri_area(pts[0], pts[1], pC)
        d = {}
        for v in data.keys():
            d0 = data[v][c_list[0]-1]
            d1 = data[v][c_list[1]-1]
            d2 = data[v][c_list[2]-1]
            d[v] = (a0*d0 + a1*d1 + a2*d2)/area
        return TriCell(d, pts)
    elif len(pts) == 4:
        pC = quad_centroid(pts[0], pts[1], pts[2], pts[3])
        area = quad_area(pts[0], pts[1], pts[2], pts[3])
        p01 = 0.5*(pts[0]+pts[1])
        p12 = 0.5*(pts[1]+pts[2])
        p23 = 0.5*(pts[2]+pts[3])
        p30 = 0.5*(pts[3]+pts[0])
        a0 = quad_area(pC, p12, pts[2], p23)
        a1 = quad_area(pC, p23, pts[3], p30)
        a2 = quad_area(pC, p30, pts[0], p01)
        a3 = quad_area(pC, p01, pts[1], p12)
        d = {}
        for v in data.keys():
            d0 = data[v][c_list[0]-1]
            d1 = data[v][c_list[1]-1]
            d2 = data[v][c_list[2]-1]
            d3 = data[v][c_list[3]-1]
            d[v] = (a0*d0 + a1*d1 + a2*d2 + a3*d3)/area
        return QuadCell(d, pts)
    else:
        print "Unknown cell type with num points= ", len(pts)
        print "Bailing out!"
        import sys
        sys.exit(1)
    return

def create_cells_from_slice(fname, var_map, scale):
    f = open(fname, 'r')
    # Read Title line and discard
    f.readline()
    # Gather variables
    var_list = []
    while 1:
        line = f.readline()
        if line.startswith('ZONE'):
            break
        if line.startswith('DATASETAUXDATA'):
            continue
        if line.startswith('VARIABLES'):
            tks = line.split()
            var_list.append(tks[2].strip('"'))
        else:
            v = line.strip().strip('"')
            var_list.append(v)
    # Read Strand and Soltime info and discard
    f.readline()
    # Read Nodes, elements -- pick these up
    line = f.readline()
    tks = line.split()
    nnodes = int(tks[0].split("=")[1][:-1])
    nelems = int(tks[1].split("=")[1][:-1])
    # Read datapacking and discard
    f.readline()
    # Read datatype and discard
    f.readline()
    # Read all data and place in tks
    tks = f.read().split()
    f.close()
    # Now pick up some data
    data = {}
    pos = 0
    for v in var_list:
        data[v] = []
        for i in range(nnodes):
            data[v].append(float(tks[pos]))
            pos = pos + 1
    cell_cnrs = []
    for i in range(nelems):
        cell_cnrs.append([int(tks[pos]), int(tks[pos+1]), int(tks[pos+2]), int(tks[pos+3])])
        pos = pos + 4
    
    cells = []
    for i in range(nelems):
        try:
            cells.append(create_cell(i, data, cell_cnrs, var_map, scale))
        except ZeroDivisionError:
            print "An attempt to create a cell with essentially zero area was made."
            
    return cells

def create_cells_from_line(fname, var_map, scale):
    f = open(fname, 'r')
    # Read Title line and discard
    f.readline()
    # Gather variables
    var_list = []
    while 1:
        line = f.readline()
        if line.startswith('ZONE'):
            break
        if line.startswith('VARIABLES'):
            tks = line.split()
            var_list.append(tks[2].strip('"'))
        else:
            v = line.strip().strip('"')
            var_list.append(v)
    # Read Strand and Soltime info and discard
    f.readline()
    # Read number of points in i
    tks = f.readline().split()
    npoints = int(tks[0].split("=")[1][:-1])
    # Read data packing and discard
    f.readline()
    # Read datatype and discard
    f.readline()
    # Read all data and place it in tks
    tks = f.read().split()
    f.close()
    # Now pick up data
    data = {}
    pos = 0
    for v in var_list:
        data[v] = []
    for i in range(npoints):
        for v in var_list:
            data[v].append(float(tks[pos]))
            pos = pos + 1
    # Create cells; first and last are special
    xlabel = var_map['x']
    ylabel = var_map['y']
    zlabel = var_map['z']
    cells = []
    # First cell:
    pt0 = Vector(data[xlabel][0], data[ylabel][0], data[zlabel][0])
    pt1 = Vector(data[xlabel][1], data[ylabel][1], data[zlabel][1])
    pts = [ pt0, 0.5*(pt0+pt1) ]
    d = {}
    for v in data.keys():
        d[v] = data[v][0]
    cells.append(AxiCell(d, pts))
    # Interior cells:
    for i in range(1,npoints-1):
        ptA = Vector(data[xlabel][i-1], data[ylabel][i-1], data[zlabel][i-1])
        ptB = Vector(data[xlabel][i], data[ylabel][i], data[zlabel][i])
        ptC = Vector(data[xlabel][i+1], data[ylabel][i+1], data[zlabel][i+1])
        pts = [ 0.5*(ptA+ptB), 0.5*(ptB+ptC) ]
        d = {}
        for v in data.keys():
            d[v] = data[v][i]
        cells.append(AxiCell(d, pts))
    # Last cell
    ptB = Vector(data[xlabel][npoints-2], data[ylabel][npoints-2], data[zlabel][npoints-2])
    ptC = Vector(data[xlabel][npoints-1], data[ylabel][npoints-1], data[zlabel][npoints-1])
    pts = [ 0.5*(ptB+ptC), ptC ]
    for v in data.keys():
        d[v] = data[v][npoints-1]
    cells.append(AxiCell(d, pts))
    
    return cells

def area(cells):
    A = 0.0
    for c in cells:
        A = A + c.area()
    return A


    




    
    
