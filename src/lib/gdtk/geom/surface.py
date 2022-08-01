# surface.py
"""
ParametricSurface classes for geometric modelling.

These are made to work like the Dlang equivalent classes.

PJ, 2020-07-05
"""
import math
from abc import ABC, abstractmethod
from copy import copy
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
               'p00', 'p10', 'p11', 'p10', 'defined_by_corners']

    def __init__(self, north=None, east=None, south=None, west=None,
                 p00=None, p10=None, p11=None, p01=None):
        """
        Initialise from edges or corner points.
        """
        if all([north, east, south, west]):
            self.north = copy(north)
            self.east = copy(east)
            self.south = copy(south)
            self.west = copy(west)
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
            self.p00 = copy(p00)
            self.p10 = copy(p10)
            self.p11 = copy(p11)
            self.p01 = copy(p01)
            self.defined_by_corners = True
        else:
            raise Exception("CoonsPatch should be defined by four edges or four corners.")
        return

    def __repr__(self):
        str = "CoonsPatch("
        if self.defined_by_corners:
            str += f"p00={self.p00}, p10={self.p10}, p11={self.p10}, p01={self.p10}"
        else:
            str += f"north={self.north}, east={self.east}, south={self.south}, west={self.west}"
        str += ")"
        return str

    def __call__(self, r, s):
        """
        Transfinite interpolation to an interior point, p.
        """
        south_r = self.south(r)
        north_r = self.north(r)
        west_s = self.west(s)
        east_s = self.east(s)
        p = (1.0-s)*south_r + s*north_r + (1.0-r)*west_s + r*east_s - \
            ((1.0-r)*(1.0-s)*self.p00 + (1.0-r)*s*self.p01 +
             r*(1.0-s)*self.p10 + r*s*self.p11)
        return p
