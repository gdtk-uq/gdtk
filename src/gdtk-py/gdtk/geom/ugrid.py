# ugrid.py
"""
UnstructuredGrid class for Eilmer. WIP.

Similar to the Dlang equivalent class.

NNG 2025
"""

import numpy as np
from abc import ABC, abstractmethod
from copy import copy
import gzip
import re
from gdtk.geom.vector3 import Vector3
from gdtk.geom.path import Path
from gdtk.geom.surface import ParametricSurface, CoonsPatch
from gdtk.geom.volume import ParametricVolume, TFIVolume
from gdtk.geom.sgrid import StructuredGrid
from gdtk.geom.cluster import *

face_name = {0:'west', 1:'east', 2:'south', 3:'north', 4:'bottom', 5:'top'}
face_types = {2:3}

class UnstructuredGrid(object):
    """
    Storage space for unstructured grid objects in Python. Note that for the moment,
    we are mostly interested in welding structured grids together.
    """
    _slots_ = ['label', 'dimensions', 'ncells', 'nvertices', 'nfaces', 'vertices',
               'celltypes', 'cell2vtx', 'boundaries', 'tags']

    def __init__(self, sg):
        """
        """
        self.label = sg.label
        self.dimensions = sg.dimensions

        if self.dimensions==2:
            self.nvertices = sg.niv * sg.njv;
            self.nfaces = (sg.niv)*(sg.njv-1) + (sg.niv-1)*(sg.njv);
            self.ncells = (sg.niv-1)*(sg.njv-1);

            self.vertices = copy(sg.vertices)

            # Eilmer's ugrid has i varying fastests, so let's transform the vertex order
            self.vertices.x = self.vertices.x.swapaxes(0, len(sg.vertices.x.shape)-1).flatten()
            self.vertices.y = self.vertices.y.swapaxes(0, len(sg.vertices.x.shape)-1).flatten()
            self.vertices.z = self.vertices.z.swapaxes(0, len(sg.vertices.z.shape)-1).flatten()

            vtx_id = np.arange(self.nvertices, dtype=np.int32).reshape((sg.njv, sg.niv))

            # FIXME: This might be a bug in the normal grid writer. It's 9? or 2?
            self.celltypes = np.zeros(self.ncells, dtype=np.int32) + 9
            self.cell2vtx = []
            for j in range(sg.njv-1):
                for i in range(sg.niv-1):
                    self.cell2vtx.append([vtx_id[j,i], vtx_id[j,i+1],
                                          vtx_id[j+1,i+1], vtx_id[j+1,i]])

            self.boundaries = []
            self.tags = []

            tag = "west";
            if sg.tags[0]!="": tag = sg.tags[0]
            bnd = []
            for j in range(sg.njv-1):
                bnd.append([vtx_id[j,0], vtx_id[j+1,0]])
            self.boundaries.append(bnd)
            self.tags.append(tag)

            tag = "east";
            if sg.tags[1]!="": tag = sg.tags[1]
            bnd = []
            for j in range(sg.njv-1):
                bnd.append([vtx_id[j,-1], vtx_id[j+1,-1]])
            self.boundaries.append(bnd)
            self.tags.append(tag)

            tag = "south";
            if sg.tags[2]!="": tag = sg.tags[2]
            bnd = []
            for i in range(sg.niv-1):
                bnd.append([vtx_id[0,i+1], vtx_id[0,i]])
            self.boundaries.append(bnd)
            self.tags.append(tag)

            tag = "north";
            if sg.tags[3]!="": tag = sg.tags[3]
            bnd = []
            for i in range(sg.niv-1):
                bnd.append([vtx_id[-1,i+1], vtx_id[-1,i]])
            self.boundaries.append(bnd)
            self.tags.append(tag)

        return

    def __repr__(self):
        str = "UnStructuredGrid("
        str += f"dimensions={self.dimensions}, nvertices={self.nvertices}, ncells={self.ncells}, nfaces={self.nfaces}"
        str += f", tags={self.tags}"
        str += ")"
        return str

    def write_to_su2_file(self, file_name):
        """
        Specific format for SU2-capable readers.
        """
        f = open(file_name, "wt")
        f.write("NDIME= {:d}\n".format(self.dimensions))
        f.write("\n")

        # TODO: Slow loop here because of ragged arrays. Might need to do better.
        f.write("NELEM= {:d}\n".format(self.ncells))
        for i in range(self.ncells):
            celltype = self.celltypes[i]
            vtxs = self.cell2vtx[i]
            fmt = ' '.join(['{:d}' for _ in range(len(vtxs)+2)])
            line = fmt.format(celltype, *vtxs, i)
            f.write(line)
            f.write("\n")

        f.write("\n")

        col = lambda arr : arr.reshape((arr.size,1))

        idx = np.arange(self.vertices.x.size)
        points = [col(self.vertices.x), col(self.vertices.y)]
        if self.dimensions>2: points.append(col(self.vertices.z))
        points.append(col(idx))
        points = np.concatenate(points, axis=1)

        npoins,ncols = points.shape
        fmt = ' '.join(['%16.16e' for _ in range(ncols-1)] + ['%d'])
        f.write("NPOIN= {:d}\n".format(self.nvertices))
        np.savetxt(f, points, fmt)
        f.write("\n")

        f.write("NMARK= {:d}\n".format(len(self.boundaries)))
        for tag,faces in zip(self.tags, self.boundaries):
            f.write("MARKER_TAG= {:s}\n".format(tag))
            f.write("MARKER_ELEMS= {:d}\n".format(len(faces)))
            for face in faces:
                nvtxs = len(face)
                facetype = face_types[nvtxs]
                fmt = ' '.join(['{:d}' for _ in range(len(face)+1)])
                line = fmt.format(facetype, *face)
                f.write(line)
                f.write("\n")
             
        f.close()
        return


if __name__=='__main__':

    L = 600.0e-3 # metres
    H = 0.20 * L

    #        d---------c
    #        |         |
    #        |         |
    #        |         |
    #        |         |
    #        a---------b
    #
    a = Vector3(x=0.0, y=0.0); b = Vector3(x=L, y=0.0);
    c = Vector3(x=L, y=H);     d = Vector3(x=0.0, y=H);

    patch = CoonsPatch(p00=a, p10=b, p11=c, p01=d)

    niv = 16+1; njv = 17+1;

    sgrd = StructuredGrid(psurf=patch, niv=niv, njv=njv,
                          tags = ['inflow', 'outflow', 'inflow', 'wall'])
    ugrd = UnstructuredGrid(sgrd)
    ugrd.write_to_su2_file('grid.su2')

    # FIXME: Unit tests back in here

