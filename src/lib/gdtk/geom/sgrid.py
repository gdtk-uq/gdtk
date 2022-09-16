# sgrid.py
"""
StructuredGrid class for Eilmer calculations.

Similar to the Dlang equivalent class.

PJ, 2020-07-05 Initial code
    2022-09-16 Add Volume grid
"""
import math
from abc import ABC, abstractmethod
from copy import copy
import gzip
from gdtk.geom.vector3 import Vector3
from gdtk.geom.path import Path
from gdtk.geom.surface import ParametricSurface
from gdtk.geom.volume import ParametricVolume
from gdtk.geom.cluster import *

class StructuredGrid():
    """
    Structured grid.
    """
    _slots_ = ['dimensions', 'niv', 'njv', 'nkv', 'vertices', 'label']

    def __init__(self, **kwargs):
        """
        Initialize by either reading a file or by discretizing a ParametricSurface.
        """
        # print("kwargs=", kwargs)
        if "psurf" in kwargs.keys():
            psurf = kwargs['psurf']
            niv = kwargs.get('niv', 1)
            njv = kwargs.get('njv', 1)
            cf_list = kwargs.get('cf_list', [None, None, None, None])
            cf_list = [cf if isinstance(cf, ClusterFunction) else LinearFunction()
                       for cf in cf_list]
            self.make_from_psurface(psurf, niv, njv, cf_list)
        elif "pvolume" in kwargs.keys():
            pvolume = kwargs['pvolume']
            niv = kwargs.get('niv', 1)
            njv = kwargs.get('njv', 1)
            nkv = kwargs.get('nkv', 1)
            cf_list = kwargs.get('cf_list', [None, None, None])
            cf_list = [cf if isinstance(cf, ClusterFunction) else LinearFunction()
                       for cf in cf_list]
            self.make_from_pvolume(pvolume, niv, njv, nkv, cf_list)
        elif "gzfile" in kwargs.keys():
            self.read_from_gzip_file(kwargs.get('gzfile'))
        else:
            raise Exception("Do not know how to make grid.")
        self.label = "unknown"
        return

    def __repr__(self):
        str = "StructuredGrid("
        str += f"dimensions={self.dimensions}, niv={self.niv}, njv={self.njv}, nkv={self.nkv}"
        # [FIX-ME] limit how much is written for large number of vertices.
        str += f", vertices={self.vertices}"
        str += ")"
        return str

    def make_from_psurface(self, psurf, niv, njv, cf_list):
        if not isinstance(psurf, ParametricSurface):
            raise Exception("Need to supply a ParametricSurface to construct the grid.")
        if niv < 2:
            raise Exception(f"niv is too small: {niv}")
        if njv < 2:
            raise Exception(f"njv is too small: {njv}")
        rNorth = cf_list[0].distribute_parameter_values(niv)
        sEast = cf_list[1].distribute_parameter_values(njv)
        rSouth = cf_list[2].distribute_parameter_values(niv)
        sWest = cf_list[3].distribute_parameter_values(njv)
        self.niv = niv
        self.njv = njv
        self.nkv = 1
        self.dimensions = 2
        dr = 1.0/(niv-1)
        ds = 1.0/(njv-1)
        self.vertices = []
        for i in range(0,niv):
            r = dr*i
            self.vertices.append([])
            for j in range(0,njv):
                s = ds*j
                sdash = (1.0-r) * sWest[j] + r * sEast[j];
                rdash = (1.0-s) * rSouth[i] + s * rNorth[i];
                self.vertices[i].append(psurf(rdash, sdash))
        return

    def make_from_pvolume(self, pvolume, niv, njv, nkv, cf_list):
        if not isinstance(psurf, ParametricVolume):
            raise Exception("Need to supply a ParametricVolume to construct the grid.")
        if niv < 2:
            raise Exception(f"niv is too small: {niv}")
        if njv < 2:
            raise Exception(f"njv is too small: {njv}")
        if nkv < 2:
            raise Exception(f"nkv is too small: {nkv}")
        self.niv = niv
        self.njv = njv
        self.nkv = nkv
        self.dimensions = 3
        # Single cluster function for each index direction.
        # This is different to the cluster-function per edge for Eilmer.
        rs = cf_list[0].distribute_parameter_values(niv)
        ss = cf_list[1].distribute_parameter_values(njv)
        ts = cf_list[2].distribute_parameter_values(nkv)
        self.vertices = []
        for i in range(0,niv):
            r = rs[i]
            self.vertices.append([])
            for j in range(0,njv):
                s = ss[j]
                self.vertices[i].append([])
                for k in range(0,nkv):
                    t = ts[k]
                    self.vertices[i][j].append(pvolume(r, s, t))
        return

    def read_from_gzip_file(self, file_name):
        with gzip.open(file_name, "rt") as f:
            line = f.readline(); items = line.split()
            assert items[1] == "1.0", "incorrect structured_grid version"
            line = f.readline(); items = line.split()
            self.label = items[1]
            line = f.readline(); items = line.split()
            self.dimensions = int(items[1])
            line = f.readline(); items = line.split()
            self.niv = int(items[1])
            line = f.readline(); items = line.split()
            self.njv = int(items[1])
            line = f.readline(); items = line.split()
            self.nkv = int(items[1])
            self.vertices = []
            if self.dimensions == 1:
                # A single list.
                for i in range(self.niv):
                    line = f.readline(); items = line.split()
                    x = float(items[0]); y = float(items[1])
                    z = float(items[2]) if len(items) > 2 else 0.0
                    self.vertices.append(Vector3(x, y, z))
            elif self.dimensions == 2:
                # A list of lists.
                for i in range(self.niv): self.vertices.append([])
                for j in range(self.njv):
                    for i in range(self.niv):
                        line = f.readline(); items = line.split()
                        x = float(items[0]); y = float(items[1])
                        z = float(items[2]) if len(items) > 2 else 0.0
                        self.vertices[i].append(Vector3(x, y, z))
            elif self.dimensions == 3:
                # A list of lists of lists.
                for i in range(self.niv):
                    self.vertices.append([])
                    for j in range(self.njv):
                        self.vertices[i].append([])
                # Get the actual data.
                for k in range(self.nkv):
                    for j in range(self.njv):
                        for i in range(self.niv):
                            line = f.readline(); items = line.split()
                            x = float(items[0]); y = float(items[1])
                            z = float(items[2]) if len(items) > 2 else 0.0
                            self.vertices[i][j].append(Vector3(x, y, z))
            else:
                raise RuntimeError("Invalid dimensions.")
        return

    def write_to_gzip_file(self, file_name):
        with gzip.open(file_name, "wt") as f:
            f.write("structured_grid 1.0\n")
            f.write(f"label: {self.label}\n")
            f.write(f"dimensions: {self.dimensions}\n")
            f.write(f"niv: {self.niv}\n")
            f.write(f"njv: {self.njv}\n")
            f.write(f"nkv: {self.nkv}\n")
            for k in range(self.nkv):
                for j in range(self.njv):
                    for i in range(self.niv):
                        vtx = self.vertices[i][j][k] if self.nkv > 1 else self.vertices[i][j]
                        f.write("%.18e %.18e %.18e\n" % (vtx.x, vtx.y, vtx.z))
        return

    def write_to_vtk_file(self, file_name):
        with open(file_name, "wt") as f:
            f.write("# vtk DataFile Version 2.0\n")
            f.write(self.label+'\n')
            f.write("ASCII\n")
            f.write("\n")
            f.write("DATASET STRUCTURED_GRID\n")
            f.write("DIMENSIONS %d %d %d\n" % (self.niv, self.njv, self.nkv))
            f.write("POINTS %d float\n" % (self.niv*self.njv*self.nkv))
            for k in range(self.nkv):
                for j in range(self.njv):
                    for i in range(self.niv):
                        vtx = self.vertices[i][j][k] if self.nkv > 1 else self.vertices[i][j]
                        f.write("%.18e %.18e %.18e\n" % (vtx.x, vtx.y, vtx.z))
        return
