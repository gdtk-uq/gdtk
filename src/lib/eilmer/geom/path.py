# path.py
"""
Path classes for geometric modelling.

These are made to work like the Dlang equivalent classes.

PJ, 2020-02-05
"""
import math
from abc import ABC, abstractmethod
from eilmer.geom.vector3 import Vector3, cross

class Path(ABC):
    """
    Base class for the family of paths.
    """
    @abstractmethod
    def __repr__(self):
        pass

    @abstractmethod
    def __call__(self, t):
        pass

    def length(self, n=20):
        """
        Crude evaluation of path length by sampling and summing.

        Subclasses may fall back to using this method.
        """
        L = 0.0
        p0 = self.__call__(0.0)
        dt = 1.0/n
        for i in range(n):
            p1 = self.__call__((i+1)*dt)
            L += abs(p1-p0)
            p0 = p1
        return L


class Line(Path):
    """
    Straight line between two points.
    """
    __slots__ = ['p0', 'p1']

    def __init__(self, p0, p1):
        self.p0 = Vector3(p0)
        self.p1 = Vector3(p1)
        return

    def __repr__(self):
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


class Bezier(Path):
    """
    Bezier curve defined on a list of points.
    """
    __slots__ = ['B']

    def __init__(self, B):
        try:
            self.B = [Vector3(p) for p in B]
        except Exception as err:
            raise ValueError(f"Was expecting to get a list of points for B, but got {B}")
        return

    def __repr__(self):
        return f"Bezier(B={self.B})"

    def __call__(self, t):
        if len(self.B) == 1: return self.B[0]
        n_order = len(self.B) - 1
        # Apply de Casteljau's algorithm.
        Q = self.B.copy() # work array will be overwritten
        for k in range(n_order):
            for i in range(n_order-k):
                Q[i] = (1.0-t)*Q[i] + t*Q[i+1]
        return Q[0]

    # end class Bezier


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
        while (i < n) and (t > self.t_values[i]): i += 1
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
    """
    A Path reparameterized such that equal increments in t correspond
    to approximately equal increments in arc length.
    """
    __slots__ = ['underlying_path', 'arc_lengths', 't_values', '_n']

    def __init__(self, underlying_path, n=1000):
        if isinstance(underlying_path, Path):
            self.underlying_path = underlying_path
            if n < 1: raise RuntimeError("Should have at least one arc-length sample.")
            self._n = n
            self.set_arc_lengths()
        else:
            raise NotImplementedError("underlying_path should be a type of Path")
        return

    def set_arc_lengths(self):
        """
        Compute the arc lengths for a number of sample points along the Path
        (in equally-spaced values of t) so that these can later be used to do
        a reverse interpolation on the evaluation parameter.
        """
        dt = 1.0/self._n
        L = 0.0
        self.arc_lengths = [0.0,]
        self.t_values = [0.0,]
        p0 = self.underlying_path(0.0)
        for i in range(1, self._n+1):
            p1 = self.underlying_path(dt*i)
            L += abs(p1-p0)
            self.arc_lengths.append(L)
            self.t_values.append(dt*i)
            p0 = p1
        return

    def underlying_t(self, t):
        """
        Search the pieces of arc length to find the piece containing the
        desired point and then interpolate the local value of t for that piece.
        """
        # The incoming parameter value, t, is proportional to arc_length fraction.
        if t < 0.0: return 0.0
        if t > 1.0: return 1.0
        L_target = t * self.arc_lengths[-1]
        # Starting from the right-hand end,
        # let's try to find a point to the left of L_target.
        # If the value is out of range, this should just result in
        # us extrapolating one of the end segments -- that's OK.
        i = self._n - 1
        while (L_target < self.arc_lengths[i]) and (i > 0): i -= 1
        frac = (L_target - self.arc_lengths[i]) / \
               (self.arc_lengths[i+1]-self.arc_lengths[i])
        return (1.0-frac)*self.t_values[i] + frac*self.t_values[i+1]

    def __repr__(self):
        return "ArcLengthParameterizedPath(underlying_path={}, n={})".format(
            self.underlying_path, self._n)

    def __call__(self, t):
        return self.underlying_path(self.underlying_t(t))

    def length(self):
        return self.underlying_path.length()

    # end class ArcLengthParameterizedPath

