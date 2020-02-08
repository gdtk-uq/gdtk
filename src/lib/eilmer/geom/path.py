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
        self.p0 = Vector3(p0)
        self.p1 = Vector3(p1)
        return

    def __str__(self):
        return "Line(p0={}, p1={})".format(self.p0, self.p1)
    
    def __call__(self, t):
        return self.p0*(1-t) + self.p1*t

    def length(self):
        return abs(self.p1 - self.p0)
    
    # end class Line


class Arc(Path):
    __slots__ = ['a', 'b', 'c']

    def __init__(self, a, b, c):
        self.a = Vector3(a)
        self.b = Vector3(b)
        self.c = Vector3(c)
        return

    def __repr__(self):
        return "Arc(a={}, b={}, c={})".format(self.a, self.b, self.c)

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
    """
    Collection of Path segments.
    """
    __slots__ = ['segments', 't_values', 'isclosed']

    def __init__(self, segments, closed=False, tolerance=1.0e-10):
        self.segments = []
        for seg in segments: self.segments.append(seg)
        self.isclosed = closed
        if self.isclosed:
            p0 = self.segments[0](0.0)
            p1 = self.segments[-1](1.0)
            if (abs(p1-p0) > tolerance):
                self.segments.append(Line(p0,p1))
        self.reset_breakpoints()
        return

    def reset_breakpoints(self):
        self.t_values = []
        t_total = 0.0
        for seg in self.segments:
            t_total += seg.length()
            self.t_values.append(t_total)
        for i in range(len(self.t_values)): self.t_values[i] /= t_total
        return
    
    def __repr__(self):
        text = "Polyline(segments=["
        n = len(self.segments)
        for i in range(n):
            text += '{}'.format(self.segments[i])
            text += ', ' if i < n-1 else ']'
        return text
    
    def __call__(self, t):
        n = len(self.segments)
        if n == 1: return self.segments[0](t)
        i = 0
        while (i < n) and (t < self.t_values[i]): i += 1
        i = min(i, n-1)
        if i == 0:
            t_local = t/self.t_values[0]
        else:
            t_local = (t-self.t_values[i-1])/(self.t_values[i]-self.t_values[i-1])
        return self.segments[i](t_local)

    def length(self):
        L = 0.0
        for seg in self.segments: L += seg.length()
        return L
    
    # end class Polyline


class ArcLengthParameterizedPath(Path):
    pass
