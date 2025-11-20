# sgrid.py
"""
StructuredGrid class for Eilmer, Chicken and Lorikeet calculations.

Similar to the Dlang equivalent class.

PJ, 2020-07-05 Initial code
    2022-09-16 Add Volume grid
NNG 2022-11-01 Arrayification (Canberra, ACT)
PJ  2022-11-03 Binary files
"""

import gzip
import warnings
from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np
from parse import parse

from gdtk.geom.cluster import ClusterFunction, LinearFunction
from gdtk.geom.surface import CoonsPatch, ParametricSurface
from gdtk.geom.vector3 import Vector3
from gdtk.geom.volume import ParametricVolume, TFIVolume

# Shorthand
CFList = Tuple[Optional[ClusterFunction], ...]


class _SGridConstructor(type):
    def __call__(cls: "type[StructuredGrid]", **kwargs) -> "StructuredGrid":
        # Previously, this was a custom __init__ to create & assign properties as needed.
        # Now, we delegate construction to several static-constructors.
        # This is a limit of Python code... We should deprecate the old constructor style ASAP

        deprecation_message = (
            "Variable '__init__' construction will be removed in the future. '"
            + cls.__name__
            + ".{}' should be used instead."
        )

        if "psurf" in kwargs.keys():
            warnings.warn(deprecation_message.format("from_psurface"), FutureWarning, stacklevel=2)
            psurf = kwargs.pop("psurf")
            instance = cls.from_psurface(psurf=psurf, **kwargs)
        elif "pvolume" in kwargs.keys():
            warnings.warn(deprecation_message.format("from_pvolume"), FutureWarning, stacklevel=2)
            pvolume = kwargs.pop("pvolume")
            instance = cls.from_pvolume(pvolume=pvolume, **kwargs)
        elif "gzfile" in kwargs.keys():
            warnings.warn(deprecation_message.format("read_from_gzip_file"), FutureWarning, stacklevel=2)
            file_name = kwargs.pop("gzfile")
            instance = cls.read_from_gzip_file(file_name=file_name, **kwargs)
        elif "binaryfile" in kwargs.keys():
            warnings.warn(deprecation_message.format("read_from_binary_file"), FutureWarning, stacklevel=2)
            instance = cls.read_from_binary_file(kwargs.get("binaryfile"))
        else:
            # A round-about way of getting back to the standard constructor
            instance = type.__call__(cls, **kwargs)

        return instance


# We use a metaclass to overload the __call__ on the *class*, to point to
# our custom constructors, rather than directly to __init__
@dataclass(slots=True)
class StructuredGrid(metaclass=_SGridConstructor):
    """
    A structured grid can be constructed on a parametric surface or volume or
    can be read from a file.

    Once in memory, the vertex coordinates are available as the numpy array, vertices.

    There is a service function to return a subgrid, and functions to write the
    grid to file.  Output formats include the native format and VTK format.
    """

    dimensions: int
    vertices: Vector3
    niv: int = 1
    njv: int = 1
    nkv: int = 1
    label: str = field(default="unknown", repr=False)
    tags: tuple[str, ...] = field(default_factory=tuple)

    @staticmethod
    def from_psurface(
        psurf: ParametricSurface,
        niv: int = 1,
        njv: int = 1,
        cf_list: CFList = (None, None, None, None),
        tags: tuple[str, str, str, str] = ("", "", "", ""),
    ):
        if not isinstance(psurf, ParametricSurface):
            raise TypeError("Need to supply a ParametricSurface to construct the grid.")
        if niv < 2:
            raise ValueError(f"niv is too small: {niv}")
        if njv < 2:
            raise ValueError(f"njv is too small: {njv}")
        if len(cf_list) < 4:
            raise ValueError("Parameter 'cf_list' must have 4 elements")

        # Filter out None values
        cf_list = tuple(cf or LinearFunction() for cf in cf_list)

        # Start with uniformly-distributed sample points.
        r = np.fromfunction(lambda i, j: i, (niv, njv), dtype=float) / (niv - 1)
        s = np.fromfunction(lambda i, j: j, (niv, njv), dtype=float) / (njv - 1)
        # Compute independent cluster function along each edge.
        # [FIX-ME] PJ 2023-04-07 Order of list items.
        # AM - Using NamedTuples would sort this out
        rNorth = cf_list[0](r)
        sEast = cf_list[1](s)
        rSouth = cf_list[2](r)
        sWest = cf_list[3](s)
        # Blend the clustered sample points from each edge of the rs unit square.
        sdash = (1.0 - r) * sWest + r * sEast
        rdash = (1.0 - s) * rSouth + s * rNorth
        # Compute the xyz spatial coordinates of the surface.
        return StructuredGrid(
            dimensions=2, niv=niv, njv=njv, nkv=1, vertices=psurf(rdash, sdash), tags=tags
        )

    @staticmethod
    def from_pvolume(
        pvolume: ParametricVolume,
        niv: int = 1,
        njv: int = 1,
        nkv: int = 1,
        cf_list: CFList = (None, None, None),
        # Ew ew ew, NamedTuples would be so much nicer
        tags: tuple[str, str, str, str, str, str] = ("", "", "", "", "", ""),
    ):
        # if not isinstance(pvolume, ParametricVolume):
        #    raise Exception("Need to supply a ParametricVolume to construct the grid.")
        if niv < 2:
            raise ValueError(f"niv is too small: {niv}")
        if njv < 2:
            raise ValueError(f"njv is too small: {njv}")
        if nkv < 2:
            raise ValueError(f"nkv is too small: {nkv}")
        if len(cf_list) < 3:
            raise ValueError("Parameter 'cf_list' must have 4 elements")

        # Filter out None values
        cf_list = tuple(cf or LinearFunction() for cf in cf_list)

        rs = np.fromfunction(lambda i, j, k: i, (niv, njv, nkv)) / (niv - 1)
        ss = np.fromfunction(lambda i, j, k: j, (niv, njv, nkv)) / (njv - 1)
        ts = np.fromfunction(lambda i, j, k: k, (niv, njv, nkv)) / (nkv - 1)
        # Single cluster function for each index direction.
        # This is different to the cluster-function per edge for Eilmer.
        rdash = cf_list[0](rs)
        sdash = cf_list[1](ss)
        tdash = cf_list[2](ts)
        return StructuredGrid(
            dimensions=3, niv=niv, njv=njv, nkv=nkv, vertices=pvolume(rdash, sdash, tdash), tags=tags
        )

    def subgrid(self, i0: int = 0, j0: int = 0, k0: int = 0, niv: int = 1, njv: int = 1, nkv: int = 1):
        """
        Returns a copy of a subgrid of vertices.
        Start at vertex i0,j0,k0 and extend for niv,njv,nkv vertices.
        """
        if self.dimensions == 3:
            newXs = self.vertices.x[i0 : i0 + niv, j0 : j0 + njv, k0 : k0 + nkv].copy()
            newYs = self.vertices.y[i0 : i0 + niv, j0 : j0 + njv, k0 : k0 + nkv].copy()
            newZs = self.vertices.z[i0 : i0 + niv, j0 : j0 + njv, k0 : k0 + nkv].copy()
            subgrid_vertices = Vector3(newXs, newYs, newZs)
        elif self.dimensions == 2:
            newXs = self.vertices.x[i0 : i0 + niv, j0 : j0 + njv].copy()
            newYs = self.vertices.y[i0 : i0 + niv, j0 : j0 + njv].copy()
            newZs = self.vertices.z[i0 : i0 + niv, j0 : j0 + njv].copy()
            subgrid_vertices = Vector3(newXs, newYs, newZs)
        elif self.dimensions == 1:
            newXs = self.vertices.x[i0 : i0 + niv].copy()
            newYs = self.vertices.y[i0 : i0 + niv].copy()
            newZs = self.vertices.z[i0 : i0 + niv].copy()
            subgrid_vertices = Vector3(newXs, newYs, newZs)
        # [TODO] PJ 2023-01-22, Think about the edge case where the number
        # of vertices in a dimension is reduced to 1.
        # The dimensions set below no longer be consistent with the array shape,
        # that is inherited from the original grid.
        if niv == 1 and njv == 1 and nkv == 1:
            dimensions = 0
        elif njv == 1 and nkv == 1:
            dimensions = 1
        elif nkv == 1:
            dimensions = 2
        else:
            dimensions = 3

        return StructuredGrid(dimensions=dimensions, niv=niv, njv=njv, nkv=nkv, vertices=subgrid_vertices)

    @staticmethod
    def from_gzip_file(file_name, tags: Optional[tuple[str, ...]] = None):
        """
        Native format for Eilmer, Chicken and Lorikeet.
        """
        with gzip.open(file_name, "rt") as file:
            lines = (line.strip() for line in file)

            (format_version,) = parse("structured_grid {}", next(lines))
            assert format_version in ("1.1", "1.0")
            (label,) = parse("label: {}", next(lines)) or ("",)
            (dimensions,) = parse("dimensions: {:d}", next(lines))
            (niv,) = parse("niv: {:d}", next(lines))
            (njv,) = parse("njv: {:d}", next(lines))
            (nkv,) = parse("nkv: {:d}", next(lines))

            # We now want to load the text associated with the vertex coordinates,
            # followed by the boundary tags, if they are present.

            x, y, z = np.loadtxt(file, max_rows=nkv * njv * niv, unpack=True)

            file_tags = []
            if format_version == "1.1":
                (ntags,) = parse("ntags: {:d}", next(lines))
                for i in range(ntags):
                    _, tag = parse("tag[{:d}]: {}", next(lines))
                    file_tags.append(tag)

        # The serialized data in the file has loops in the VTK index order,
        # with k as the outer loop, then j and then i as the innermost loop
        x = x.reshape((nkv, njv, niv)).transpose().squeeze()
        y = y.reshape((nkv, njv, niv)).transpose().squeeze()
        z = z.reshape((nkv, njv, niv)).transpose().squeeze()

        return StructuredGrid(
            dimensions=dimensions,
            vertices=Vector3(x=x, y=y, z=z),
            niv=niv,
            njv=njv,
            nkv=nkv,
            tags=(tags or tuple(file_tags)),
            label=label,
        )

    def read_from_binary_file(self, file_name):
        """
        Bare-bones binary reading for Chicken and Lorikeet.
        Not for Eilmer. No reading of boundary tags.
        """
        data = np.fromfile(file_name, dtype=float)
        data = data.reshape((data.shape[0] // 3, 3))
        self.dimensions = int(data[0, 0])
        self.niv = int(data[1, 0])
        self.njv = int(data[1, 1])
        self.nkv = int(data[1, 2])
        x = data[2:, 0]
        y = data[2:, 1]
        z = data[2:, 2]
        if self.dimensions == 1:
            pass
        elif self.dimensions == 2:
            x = x.reshape((self.njv, self.niv)).transpose()
            y = y.reshape((self.njv, self.niv)).transpose()
            z = z.reshape((self.njv, self.niv)).transpose()
        elif self.dimensions == 3:
            x = x.reshape((self.nkv, self.njv, self.niv)).transpose()
            y = y.reshape((self.nkv, self.njv, self.niv)).transpose()
            z = z.reshape((self.nkv, self.njv, self.niv)).transpose()
        else:
            raise RuntimeError("Invalid dimensions.")
        self.vertices = Vector3(x=x.copy(), y=y.copy(), z=z.copy())
        return

    def write_to_gzip_file(self, file_name, format_version="1.0"):
        """
        Native format for Eilmer, Chicken and Lorikeet.
        """
        f = gzip.open(file_name, "wt")
        f.write("structured_grid %s\n" % format_version)
        f.write(f"label: {self.label}\n")
        f.write(f"dimensions: {self.dimensions}\n")
        f.write(f"niv: {self.niv}\n")
        f.write(f"njv: {self.njv}\n")
        f.write(f"nkv: {self.nkv}\n")
        # The serialized data in the file has loops in the VTK index order,
        # with k as the outer loop, then j and then i as the innermost loop
        data = np.zeros((self.nkv * self.njv * self.niv, 3), dtype=float)
        data[:, 0] = self.vertices.x.transpose().flatten()
        data[:, 1] = self.vertices.y.transpose().flatten()
        data[:, 2] = self.vertices.z.transpose().flatten()
        np.savetxt(f, data)
        if format_version == "1.1":
            f.write(f"ntags: {len(self.tags)}\n")
            for i, tag in enumerate(self.tags):
                f.write(f"tag[{i}]: {tag}\n")
        f.close()
        return

    def write_to_binary_file(self, file_name):
        """
        Bare-bones binary writing for Chicken and Lorikeet.
        Not for Eilmer. No writing of boundary tags.
        """
        data = np.zeros((self.nkv * self.njv * self.niv + 2, 3), dtype=float)
        # Pack the metadata into the first two rows.
        data[0, :] = [float(self.dimensions), 0.0, 0.0]
        data[1, :] = [float(self.niv), float(self.njv), float(self.nkv)]
        # Pack the main data into the remaining rows.
        data[2:, 0] = self.vertices.x.transpose().flatten()
        data[2:, 1] = self.vertices.y.transpose().flatten()
        data[2:, 2] = self.vertices.z.transpose().flatten()
        data.tofile(file_name)
        return

    def write_to_vtk_file(self, file_name):
        """
        Generic format for VTK-capable readers.
        """
        f = open(file_name, "wt")
        f.write("# vtk DataFile Version 2.0\n")
        f.write(self.label + "\n")
        f.write("ASCII\n")
        f.write("\n")
        f.write("DATASET STRUCTURED_GRID\n")
        f.write("DIMENSIONS %d %d %d\n" % (self.niv, self.njv, self.nkv))
        f.write("POINTS %d float\n" % (self.niv * self.njv * self.nkv))
        # The serialized data in the file has k as the outer loop,
        # then j and then i as the innermost loop.
        data = np.zeros((self.nkv * self.njv * self.niv, 3), dtype=float)
        data[:, 0] = self.vertices.x.transpose().flatten()
        data[:, 1] = self.vertices.y.transpose().flatten()
        data[:, 2] = self.vertices.z.transpose().flatten()
        np.savetxt(f, data)
        f.close()
        return


class StructuredGrid_old:
    """
    Structured grid.
    """

    _slots_ = ["dimensions", "niv", "njv", "nkv", "vertices", "label"]

    def __init__(self, **kwargs):
        """
        Initialize by either reading a file or by discretizing a ParametricSurface.
        """
        # print("kwargs=", kwargs)
        if "psurf" in kwargs.keys():
            psurf = kwargs["psurf"]
            niv = kwargs.get("niv", 1)
            njv = kwargs.get("njv", 1)
            cf_list = kwargs.get("cf_list", [None, None, None, None])
            cf_list = [cf if isinstance(cf, ClusterFunction) else LinearFunction() for cf in cf_list]
            self.from_psurface(psurf, niv, njv, cf_list)
        elif "pvolume" in kwargs.keys():
            pvolume = kwargs["pvolume"]
            niv = kwargs.get("niv", 1)
            njv = kwargs.get("njv", 1)
            nkv = kwargs.get("nkv", 1)
            cf_list = kwargs.get("cf_list", [None, None, None])
            cf_list = [cf if isinstance(cf, ClusterFunction) else LinearFunction() for cf in cf_list]
            self.from_pvolume(pvolume, niv, njv, nkv, cf_list)
        elif "gzfile" in kwargs.keys():
            self.read_from_gzip_file(kwargs.get("gzfile"))
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

    def from_psurface(self, psurf, niv, njv, cf_list):
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
        dr = 1.0 / (niv - 1)
        ds = 1.0 / (njv - 1)
        self.vertices = []
        for i in range(0, niv):
            r = dr * i
            self.vertices.append([])
            for j in range(0, njv):
                s = ds * j
                sdash = (1.0 - r) * sWest[j] + r * sEast[j]
                rdash = (1.0 - s) * rSouth[i] + s * rNorth[i]
                self.vertices[i].append(psurf(rdash, sdash))
        return

    def from_pvolume(self, pvolume, niv, njv, nkv, cf_list):
        if not isinstance(pvolume, ParametricVolume):
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
        for i in range(0, niv):
            r = rs[i]
            self.vertices.append([])
            for j in range(0, njv):
                s = ss[j]
                self.vertices[i].append([])
                for k in range(0, nkv):
                    t = ts[k]
                    self.vertices[i][j].append(pvolume(r, s, t))
        return

    def read_from_gzip_file(self, file_name):
        with gzip.open(file_name, "rt") as f:
            line = f.readline()
            items = line.split()
            assert items[1] == "1.0", "incorrect structured_grid version"
            line = f.readline()
            items = line.split()
            self.label = items[1]
            line = f.readline()
            items = line.split()
            self.dimensions = int(items[1])
            line = f.readline()
            items = line.split()
            self.niv = int(items[1])
            line = f.readline()
            items = line.split()
            self.njv = int(items[1])
            line = f.readline()
            items = line.split()
            self.nkv = int(items[1])
            self.vertices = []
            if self.dimensions == 1:
                # A single list.
                for i in range(self.niv):
                    line = f.readline()
                    items = line.split()
                    x = float(items[0])
                    y = float(items[1])
                    z = float(items[2]) if len(items) > 2 else 0.0
                    self.vertices.append(Vector3(x, y, z))
            elif self.dimensions == 2:
                # A list of lists.
                for i in range(self.niv):
                    self.vertices.append([])
                for j in range(self.njv):
                    for i in range(self.niv):
                        line = f.readline()
                        items = line.split()
                        x = float(items[0])
                        y = float(items[1])
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
                            line = f.readline()
                            items = line.split()
                            x = float(items[0])
                            y = float(items[1])
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
            f.write(self.label + "\n")
            f.write("ASCII\n")
            f.write("\n")
            f.write("DATASET STRUCTURED_GRID\n")
            f.write("DIMENSIONS %d %d %d\n" % (self.niv, self.njv, self.nkv))
            f.write("POINTS %d float\n" % (self.niv * self.njv * self.nkv))
            for k in range(self.nkv):
                for j in range(self.njv):
                    for i in range(self.niv):
                        vtx = self.vertices[i][j][k] if self.nkv > 1 else self.vertices[i][j]
                        f.write("%.18e %.18e %.18e\n" % (vtx.x, vtx.y, vtx.z))
        return


if __name__ == "__main__":
    p0 = Vector3(x=0.0, y=0.0)
    p1 = Vector3(x=1.0, y=0.0)
    p2 = Vector3(x=1.0, y=1.0)
    p3 = Vector3(x=0.0, y=1.0)
    patch = CoonsPatch(p00=p0, p10=p1, p11=p2, p01=p3)
    grid = StructuredGrid.from_psurface(psurf=patch, niv=11, njv=11)
    grid2 = StructuredGrid_old(psurf=patch, niv=11, njv=11)
    grid.write_to_gzip_file("test.gz")
    grid2.write_to_gzip_file("test2.gz")

    gridb = StructuredGrid.from_gzip_file(file_name="test.gz")
    assert grid.vertices.x.shape == gridb.vertices.x.shape
    assert np.all(np.isclose(grid.vertices.x, gridb.vertices.x))
    assert np.all(np.isclose(grid.vertices.y, gridb.vertices.y))
    assert np.all(np.isclose(grid.vertices.z, gridb.vertices.z))

    p0 = Vector3(x=0.0, y=0.0, z=0.0)
    p1 = Vector3(x=1.0, y=0.0, z=0.0)
    p2 = Vector3(x=1.0, y=1.0, z=0.0)
    p3 = Vector3(x=0.0, y=1.0, z=0.0)
    p4 = Vector3(x=0.0, y=0.0, z=1.0)
    p5 = Vector3(x=1.0, y=0.0, z=1.0)
    p6 = Vector3(x=1.0, y=1.0, z=1.0)
    p7 = Vector3(x=0.0, y=1.0, z=1.0)
    volume = TFIVolume(p000=p0, p100=p1, p110=p2, p010=p3, p001=p4, p101=p5, p111=p6, p011=p7)
    grid3 = StructuredGrid.from_pvolume(volume, niv=11, njv=11, nkv=11)
    grid4 = StructuredGrid_old(pvolume=volume, niv=11, njv=11, nkv=11)
    grid3.write_to_gzip_file("test3.gz")
    grid4.write_to_gzip_file("test4.gz")

    grid3b = StructuredGrid.from_gzip_file("test3.gz")
    assert grid3.vertices.x.shape == grid3b.vertices.x.shape
    assert np.all(np.isclose(grid3.vertices.x, grid3b.vertices.x))
    assert np.all(np.isclose(grid3.vertices.y, grid3b.vertices.y))
    assert np.all(np.isclose(grid3.vertices.z, grid3b.vertices.z))
