# surface.py
"""
ParametricSurface classes for geometric modelling.

These are made to work like the Dlang equivalent classes.

PJ, 2020-07-05
"""
import numpy as np
from abc import ABC, abstractmethod
from copy import deepcopy
from gdtk.geom.vector3 import Vector3, approxEqualVectors
from gdtk.geom.path import Line


class ParametricSurface(ABC):
    """
    Base class for the family of surfaces.
    """
    @abstractmethod
    def __repr__(self):
        pass

    @abstractmethod
    def __call__(self, r, s):
        pass


class CoonsPatch(ParametricSurface):
    """
    Surface using transfinite interpolation of the edges.
    """
    _slots_ = ['north', 'east', 'south', 'west',
               'p00', 'p10', 'p11', 'p01', 'defined_by_corners',
               'offset']

    def __init__(self, north=None, east=None, south=None, west=None,
                 p00=None, p10=None, p11=None, p01=None, offset=Vector3(0,0,0)):
        """
        Initialise from edges or corner points.
        """
        if all([north, east, south, west]):
            self.north = deepcopy(north)
            self.east = deepcopy(east)
            self.south = deepcopy(south)
            self.west = deepcopy(west)
            self.p00 = self.south(0.0)
            self.p10 = self.south(1.0)
            self.p01 = self.north(0.0)
            self.p11 = self.north(1.0)
            p00_alt = self.west(0.0)
            p10_alt = self.east(0.0)
            p01_alt = self.west(1.0)
            p11_alt = self.east(1.0)
            if not approxEqualVectors(self.p00, p00_alt):
                raise Exception("CoonsPatch open corner p00={self.p00}, {p00_alt}")
            if not approxEqualVectors(self.p10, p10_alt):
                raise Exception("CoonsPatch open corner p10={self.p10}, {p10_alt}")
            if not approxEqualVectors(self.p01, p01_alt):
                raise Exception("CoonsPatch open corner p01={self.p01}, {p01_alt}")
            if not approxEqualVectors(self.p11, p11_alt):
                raise Exception("CoonsPatch open corner p11={self.p11}, {p11_alt}")
            self.defined_by_corners = False
        elif all([p00, p10, p11, p01]):
            self.north = Line(p01, p11)
            self.east = Line(p10, p11)
            self.south = Line(p00, p10)
            self.west = Line(p00, p01)
            self.p00 = deepcopy(p00)
            self.p10 = deepcopy(p10)
            self.p11 = deepcopy(p11)
            self.p01 = deepcopy(p01)
            self.defined_by_corners = True
        else:
            raise Exception("CoonsPatch should be defined by four edges or four corners.")
        self.offset = offset
        return

    def __repr__(self):
        str = "CoonsPatch("
        if self.defined_by_corners:
            str += f"p00={self.p00}, p10={self.p10}, p11={self.p10}, p01={self.p10}"
        else:
            str += f"north={self.north}, east={self.east}, south={self.south}, west={self.west}"
        str += f", offset={self.offset})"
        return str

    def __call__(self, r, s):
        """
        Transfinite interpolation to an interior point, p.
        """
        south_r = self.south(r)
        north_r = self.north(r)
        west_s = self.west(s)
        east_s = self.east(s)
        p = south_r*(1.0-s) + north_r*s + west_s*(1.0-r) + east_s*r - \
            (self.p00*(1.0-r)*(1.0-s) + self.p01*(1.0-r)*s +
             self.p10*r*(1.0-s) + self.p11*r*s) + self.offset
        return p

    def __add__(self, offset):
        """
        Returns a copy of the original, displaced by a Vector3 object.
        """
        if not isinstance(offset, Vector3):
            raise Exception(f"Cannot add a {type(offset)} to a CoonsPatch.")
        new_patch = deepcopy(self)
        new_patch.offset += offset
        return new_patch

    def __sub__(self, offset):
        """
        Returns a copy of the original, displaced by a Vector3 object.
        """
        if not isinstance(offset, Vector3):
            raise Exception(f"Cannot subtract a {type(offset)} from a CoonsPatch.")
        new_patch = deepcopy(self)
        new_patch.offset -= offset
        return new_patch


if __name__=='__main__':
    p0 = Vector3(x=0.0, y=0.0)
    p1 = Vector3(x=1.0, y=0.0)
    p2 = Vector3(x=1.0, y=1.0)
    p3 = Vector3(x=0.0, y=1.0)
    patch = CoonsPatch(p00=p0, p10=p1, p11=p2, p01=p3)
    x = patch(0.5, 0.5)
    xx= patch(np.array([0.5, 0.5]), np.array([0.5, 0.5]))
    assert(np.isclose(x.x, xx.x[0]))
    assert(np.isclose(x.y, xx.y[0]))
    #
    n = Line(p3, p2)
    s = Line(p0, p1)
    e = Line(p1, p2)
    w = Line(p0, p3)
    patch = CoonsPatch(north=n, south=s, east=e, west=w)
    x = patch(0.5, 0.5)
    xx= patch(np.array([0.5, 0.5]), np.array([0.5, 0.5]))
    print(xx)
    assert(np.isclose(x.x, xx.x[0]))
    assert(np.isclose(x.y, xx.y[0]))
    #
    patch2 = patch + Vector3(1.0, 0.0, 0.0)
    x2 = patch2(0.5, 0.5)
    print(x2)
    assert(np.isclose(x2.x, x.x+1.0))
    assert(np.isclose(x2.y, x.y))
    assert(np.isclose(x2.z, x.z))
    #
    patch3 = patch - Vector3(0.0, 1.0, 0.0)
    x3 = patch3(0.5, 0.5)
    print(x3)
    assert(np.isclose(x3.x, x.x))
    assert(np.isclose(x3.y, x.y-1.0))
    assert(np.isclose(x3.z, x.z))
