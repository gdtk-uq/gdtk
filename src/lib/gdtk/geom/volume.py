# volume.py
"""
ParametricVolume classes for geometric modelling.

These are made to work like the Dlang equivalent classes.

Notes: NNG made minor changes to this module (and surface.py)
to allow for arrayification as part of the chicken project.
Principly, this is because of a strange quirk of python's
operator overloading, where the __mult__ method of the
leftmost object in a a*b operation is called first.

For the geom library, we always want GeomObject*array
rather than array*Geomobject, to make sure that GeomObject's
__mult__ is called, which correctly applies the array
to the individual components of the object.

PJ, 2022-09-16
NNG, 2022-11-01
"""
import numpy as np
from abc import ABC, abstractmethod
from copy import copy
from gdtk.geom.vector3 import Vector3, approxEqualVectors
from gdtk.geom.path import Path, Line
from gdtk.geom.surface import ParametricSurface, CoonsPatch


class ParametricVolume(ABC):
    """
    Base class for the family of volumes.
    """
    @abstractmethod
    def __repr__(self):
        pass

    @abstractmethod
    def __call__(self, r, s):
        pass


class TFIVolume(ParametricVolume):
    """
    Volume using transfinite interpolation of the surfaces.

    The volume may be constructed by supplying 6 ParametricSurfaces or
    by supplying 8 Vector3 objects for the corner locations.
    """
    _slots_ = ['iminus', 'iplus', 'jminus', 'jplus', 'kminus', 'kplus',
               'p000', 'p100', 'p110', 'p010',
               'p001', 'p101', 'p111', 'p011',
               'defined_by_corners']

    def __init__(self, iminus=None, iplus=None, jminus=None, jplus=None, kminus=None, kplus=None,
                 p000=None, p100=None, p110=None, p010=None,
                 p001=None, p101=None, p111=None, p011=None):
        """
        Initialise from surfaces or corner points.
        """
        if all([iminus, iplus, jminus, jplus, kminus, kplus]):
            # Copy the bounding surfaces that have been provided.
            self.iminus = copy(iminus)
            self.iplus = copy(iplus)
            self.jminus = copy(jminus)
            self.jplus = copy(jplus)
            self.kminus = copy(kminus)
            self.kplus = copy(kplus)
            #
            # Compute corners from bounding surfaces.
            self.p000 = self.kminus(0.0, 0.0)
            self.p100 = self.kminus(1.0, 0.0)
            self.p110 = self.kminus(1.0, 1.0)
            self.p010 = self.kminus(0.0, 1.0)
            self.p001 = self.kplus(0.0, 0.0)
            self.p101 = self.kplus(1.0, 0.0)
            self.p111 = self.kplus(1.0, 1.0)
            self.p011 = self.kplus(0.0, 1.0)
            #
            p000_alt1 = self.iminus(0.0, 0.0)
            p000_alt2 = self.jminus(0.0, 0.0)
            if not (approxEqualVectors(self.p000, p000_alt1) and approxEqualVectors(self.p000, p000_alt2)):
                raise Exception("TFIVolume open corner p000={self.p000}, {p000_alt1}, {p000_alt2}")
            p100_alt1 = self.iplus(0.0, 0.0)
            p100_alt2 = self.jminus(1.0, 0.0)
            if not (approxEqualVectors(self.p100, p100_alt1) and approxEqualVectors(self.p100, p100_alt2)):
                raise Exception("TFIVolume open corner p100={self.p100}, {p100_alt1}, {p100_alt2}")
            p110_alt1 = self.iplus(1.0, 0.0)
            p110_alt2 = self.jplus(1.0, 0.0)
            if not (approxEqualVectors(self.p110, p110_alt1) and approxEqualVectors(self.p110, p110_alt2)):
                raise Exception("TFIVolume open corner p110={self.p110}, {p110_alt1}, {p110_alt2}")
            p010_alt1 = self.iminus(1.0, 0.0)
            p010_alt2 = self.jplus(0.0, 0.0)
            if not (approxEqualVectors(self.p010, p010_alt1) and approxEqualVectors(self.p010, p010_alt2)):
                raise Exception("TFIVolume open corner p010={self.p010}, {p010_alt1}, {p010_alt2}")
            p001_alt1 = self.iminus(0.0, 1.0)
            p001_alt2 = self.jminus(0.0, 1.0)
            if not (approxEqualVectors(self.p001, p001_alt1) and approxEqualVectors(self.p001, p001_alt2)):
                raise Exception("TFIVolume open corner p001={self.p001}, {p001_alt1}, {p001_alt2}")
            p101_alt1 = self.iplus(0.0, 1.0)
            p101_alt2 = self.jminus(1.0, 1.0)
            if not (approxEqualVectors(self.p101, p101_alt1) and approxEqualVectors(self.p101, p101_alt2)):
                raise Exception("TFIVolume open corner p101={self.p101}, {p101_alt1}, {p101_alt2}")
            p111_alt1 = self.iplus(1.0, 1.0)
            p111_alt2 = self.jplus(1.0, 1.0)
            if not (approxEqualVectors(self.p111, p111_alt1) and approxEqualVectors(self.p111, p111_alt2)):
                raise Exception("TFIVolume open corner p111={self.p111}, {p111_alt1}, {p111_alt2}")
            p011_alt1 = self.iminus(1.0, 1.0)
            p011_alt2 = self.jplus(0.0, 1.0)
            if not (approxEqualVectors(self.p011, p011_alt1) and approxEqualVectors(self.p011, p011_alt2)):
                raise Exception("TFIVolume open corner p011={self.p011}, {p011_alt1}, {p011_alt2}")
            self.defined_by_corners = False
        elif all([p000, p100, p110, p010, p001, p101, p111, p011]):
            self.iminus = CoonsPatch(p00=p000, p10=p010, p11=p011, p01=p001)
            self.iplus = CoonsPatch(p00=p100, p10=p110, p11=p111, p01=p101)
            self.jminus = CoonsPatch(p00=p000, p10=p100, p11=p101, p01=p001)
            self.jplus = CoonsPatch(p00=p010, p10=p110, p11=p111, p01=p011)
            self.kminus = CoonsPatch(p00=p000, p10=p100, p11=p110, p01=p010)
            self.kplus = CoonsPatch(p00=p001, p10=p101, p11=p111, p01=p011)
            self.p000 = copy(p000)
            self.p100 = copy(p100)
            self.p110 = copy(p110)
            self.p010 = copy(p010)
            self.p001 = copy(p001)
            self.p101 = copy(p101)
            self.p111 = copy(p111)
            self.p011 = copy(p011)
            self.defined_by_corners = True
        else:
            raise Exception("TFIVolume should be defined by six faces or eight corners.")
        return

    def __repr__(self):
        str = "TFIVolume("
        if self.defined_by_corners:
            str += f"p000={self.p000}, p100={self.p100}, p110={self.p100}, p010={self.p100}, "
            str += f"p001={self.p001}, p101={self.p101}, p111={self.p101}, p011={self.p101}"
        else:
            str += f"iminus={self.iminus}, iplus={self.iplus}, "
            str += f"jminus={self.jminus}, jplus={self.jplus}, "
            str += f"kminus={self.kminus}, kplus={self.kplus}, "
        str += ")"
        return str

    def __call__(self, r, s, t):
        """
        Transfinite interpolation to an interior point, p.
        """
        iminus_st = self.iminus(s, t)
        iplus_st = self.iplus(s, t)
        jminus_rt = self.jminus(r, t)
        jplus_rt = self.jplus(r, t)
        kminus_rs = self.kminus(r, s)
        kplus_rs = self.kplus(r, s)
        omr = 1.0-r; oms = 1.0-s; omt = 1.0-t;
        BigC = self.p000*(omr*oms*omt) + self.p001*(omr*oms*t) + \
            self.p010*(omr*s*omt) + self.p011*(omr*s*t) + \
            self.p100*(r*oms*omt) + self.p101*(r*oms*t) + \
            self.p110*(r*s*omt) + self.p111*(r*s*t)
        p_rst = 0.5*(iminus_st*omr + iplus_st*r + \
                     jminus_rt*oms + jplus_rt*s + \
                     kminus_rs*omt + kplus_rs*t) - 0.5*BigC
        return p_rst


class SweptSurfaceVolume(ParametricVolume):
    """
    Volume constructed from a kminus face and an edge from p0 (p000) to p4 (p001).

    The resulting volume has p000 located at edge04(0.0) and p001 at edge(1.0).
    """

    def __init__(self, face0123, edge04):
        if not isinstance(face0123, ParametricSurface):
            raise RuntimeError("SweptSurfaceVolume: face0123 is not a ParametricSurface.")
        if not isinstance(edge04, Path):
            raise RuntimeError("SweptSurfaceVolume: edge04 is not a Path object.")
        self.face0123 = copy(face0123)
        self.edge04 = copy(edge04)
        return

    def __repr__(self):
        return f"SweptSurfaceVolume(face0123={self.face0123}, edge04={self.edge04})"

    def __call__(self, r, s, t):
        """
        Evaluate point(s) in the volume.
        """
        return self.edge04(t) + self.face0123(r, s) - self.face0123(0.0, 0.0)


if __name__=='__main__':
    p0 = Vector3(x=0.0, y=0.0, z=0.0)
    p1 = Vector3(x=1.0, y=0.0, z=0.0)
    p2 = Vector3(x=1.0, y=1.0, z=0.0)
    p3 = Vector3(x=0.0, y=1.0, z=0.0)
    p4 = Vector3(x=0.0, y=0.0, z=1.0)
    p5 = Vector3(x=1.0, y=0.0, z=1.0)
    p6 = Vector3(x=1.0, y=1.0, z=1.0)
    p7 = Vector3(x=0.0, y=1.0, z=1.0)
    volume = TFIVolume(p000=p0, p100=p1, p110=p2, p010=p3,
                       p001=p4, p101=p5, p111=p6, p011=p7)
    #
    x = volume(0.5, 0.5, 0.5)
    xx= volume(np.array([0.5, 0.5]), np.array([0.5, 0.5]), np.array([0.5, 0.5]))
    #
    assert(np.isclose(x.x, xx.x[0]))
    assert(np.isclose(x.y, xx.y[0]))
    assert(np.isclose(x.z, xx.z[0]))
    #
    face0123 = CoonsPatch(p00=p0, p10=p1, p11=p2, p01=p3)
    edge04 = Line(p0, p4)
    volume2 = SweptSurfaceVolume(face0123=face0123, edge04=edge04)
    # print("volume2=", volume2)
    #
    x2 = volume2(0.5, 0.5, 0.5)
    xx2 = volume2(np.array([0.5, 0.5]), np.array([0.5, 0.5]), np.array([0.5, 0.5]))
    # print("x2=", x2, " xx2=", xx2)
    assert(np.isclose(x2.x, 0.5))
    assert(np.isclose(x2.y, 0.5))
    assert(np.isclose(x2.z, 0.5))
    assert(np.isclose(xx2.x[0], 0.5))
    assert(np.isclose(xx2.y[0], 0.5))
    assert(np.isclose(xx2.z[0], 0.5))
