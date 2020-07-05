# sgrid.py
"""
StructuredGrid class for Eilmer calculations.

Similar to the Dlang equivalent class.

PJ, 2020-07-05
"""
import math
from abc import ABC, abstractmethod
from copy import copy
from eilmer.geom.vector3 import Vector3
from eilmer.geom.path import Path
from eilmer.geom.surface import ParametricSurface


class StructuredGrid():
    """
    Structured grid.
    """
    _slots_ = ['niv', 'njv', 'nkv', 'vertices', 'label']

    def __init__(self, **kwargs):
        """
        Initialize by either reading a file or by discretizing a ParametricSurface.
        """
        print("kwargs=", kwargs)
        if "psurf" in kwargs.keys():
            psurf = kwargs['psurf']
            niv = kwargs.get('niv', 1)
            njv = kwargs.get('njv', 1)
            cf_list = kwargs.get('cf_list', [])
            self.make_from_psurface(psurf, niv, njv, cf_list)
        elif "gzfile" in kwargs.keys():
            print("Read grid from gzip file")
        else:
            raise Exception("Do not know how to make grid.")
        self.label = "unknown"
        return

    def __repr__(self):
        str = "StructuredGrid("
        str += f"niv={self.niv}, njv={self.njv}, nkv={self.nkv}"
        # [FIX-ME] limit for large number of vertices.
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
        self.niv = niv
        self.njv = njv
        self.nkv = 1
        dr = 1.0/(niv-1)
        ds = 1.0/(njv-1)
        self.vertices = []
        for i in range(0,niv):
            r = dr*i
            self.vertices.append([])
            for j in range(0,njv):
                s = ds*j
                self.vertices[i].append(psurf(r,s))
        return

    def read_from_gzip_file(self, file_name):
        # [FIX-ME]
        return

    def write_to_gzip_file(self, file_name):
        # [FIX-ME]
        return

    def write_to_vtk_file(self, file_name):
        f = open(file_name, "w");
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
                    if self.nkv > 1:
                        vtx = self.vertices[i][j][k]
                    else:
                        vtx = self.vertices[i][j]
                    f.write("%.18e %.18e %.18e\n" % (vtx.x, vtx.y, vtx.z))
        f.close()
        return
