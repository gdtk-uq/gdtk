# path.py
"""
Path classes for geometric modelling.

These are made to work like the Dlang equivalent classes.

PJ, 2020-02-05
"""
import math
from eilmer.geom.vector3 import Vector3, cross

class Path(object):
    pass


class Line(Path):
    """
    Straight line between two points.
    """
    __slots__ = ['p0', 'p1']

    def __init__(self, p0, p1):
        self.p0 = p0
        self.p1 = p1
        return

    def __str__(self):
        text = "Line(p0=%s, p1=%s)" % (self.p0, self.p1)
        return text
    
    def __call__(self, t):
        return self.p0*(1-t) + self.p1*t

    def length(self):
        return abs(self.p1 - self.p0)
    
    # end class Line


class Arc(Path):
    __slots__ = ['a', 'b', 'c']

    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c
        return

    def __str__(self):
        text = "Arc(a=%s, b=%s, c=%s)" % (self.a, self.b, self.c)
        return text

    def __call__(self, t):
        p, L = self.evaluate_position_and_length(t)
        return p

    def length(self):
        p, L = self.evaluate_position_and_length(1.0)
        return L

    def evaluate_position_and_length(self, t):
        L = 0.0;
        ca = self.a - self.c; ca_mag = abs(ca)
        cb = self.b - self.c; cb_mag = abs(cb)
        if abs(ca_mag - cb_mag) > 1.0e-5:
            raise Exception("Arc: radii do not match ca=%s cb=%s" % (ca, cb))
        # First vector in plane.
        tangent1 = Vector3(ca); tangent1.normalize() 
        # Compute unit normal to plane of all three points.
        n = cross(ca, cb)
        if abs(n) > 0.0:
            n.normalize()
        else:
            raise Exception("Arc: cannot find plane of three points.")
        # Third (orthogonal) vector is in the original plane.
        tangent2 = cross(n, tangent1) 
        # Now transform to local coordinates so that we can do 
        # the calculation of the point along the arc in 
        # the local xy-plane, with ca along the x-axis.
        cb_local = Vector3(cb)
        cb_local.transform_to_local_frame(tangent1, tangent2, n)
        if abs(cb_local.z) > 1.0e-6:
            raise Exception("Arc: problem with transformation cb_local=%s" % cb_local)
        # Angle of the final point on the arc is in the range -pi < th <= +pi.
        theta = math.atan2(cb_local.y, cb_local.x)
        # The length of the circular arc.
        L = theta * cb_mag
        # Move the second point around the arc in the local xy-plane.
        theta *= t
        loc = Vector3(math.cos(theta)*cb_mag, math.sin(theta)*cb_mag, 0.0)
        # Transform back to global xyz coordinates
        # and remember to add the centre coordinates.
        loc.transform_to_global_frame(tangent1, tangent2, n, self.c);
        return loc, L
    
    # end class Arc


class Polyline(Path):
    pass


class ArcLengthParameterizedPath(Path):
    pass
