# path.py
"""
Path classes for geometric modelling.

These are made to work like the Dlang equivalent classes.

PJ, 2020-02-05
Arrayification by NNG, 2022-10-07 (Pokolbin, NSW)
FnPath by PJ, 2022-11-14
"""
import numpy as np
from abc import ABC, abstractmethod
from gdtk.geom.vector3 import Vector3, cross


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


class FnPath(Path):
    """
    A path defined by a user-supplied Python function.
    """
    __slots__ = ['fn']

    def __init__(self, fn):
        """
        Construct from a user-supplied function fn(t) that yields Vector3 values.

        Note that, as of Nick's arrayification of the gridding functions in Nov 2022,
        this user-supplied function needs to be able to accept a numpy.ndarray of t values.
        The shape of the array may have one, two or three dimensions, depending on which
        grid-generation function is calling it.
        """
        if not callable(fn): raise RuntimeError("FnPath: fn is not callable.")
        if not isinstance(fn(0.0), Vector3):
            raise RuntimeError("FnPath: fn did not return a Vector3 instance for t=0.0.")
        t_ends = np.array([0.0, 1.0])
        p_ends = fn(t_ends)
        if not isinstance(p_ends, Vector3):
            raise RuntimeError("FnPath: fn did not return a Vector3 instance for array of t values.")
        self.fn = fn
        return

    def __repr__(self):
        return "FnPath(fn={})".format(self.fn)

    def __call__(self, t):
        return self.fn(t)

    # end class FnPath


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
        if np.any(np.absolute(cb_local.z)) > 1.0e-6:
            raise Exception("Arc: problem with transformation cb_local=%s" % cb_local)
        # Angle of the final point on the arc is in the range -pi < th <= +pi.
        theta = np.arctan2(cb_local.y, cb_local.x)
        # The length of the circular arc.
        L = theta * cb_mag
        # Move the second point around the arc in the local xy-plane.
        theta *= t
        loc = Vector3(np.cos(theta)*cb_mag, np.sin(theta)*cb_mag, 0.0*theta)
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
                Q[i] = Q[i]*(1.0-t) + Q[i+1]*t
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
        self.t_values = [0.0]
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

        f = self.segments[0](t) # Hmmmmmmmm FIXME
        f*= 0.0

        for i, tl, tu in zip(range(n),self.t_values[:-1], self.t_values[1:]):
            t_local = (t-tl)/(tu-tl)
            isokay = np.logical_and(t_local>=-1e-9, t_local<1.000000001)
            f += self.segments[i](t_local)*isokay
        return f

    def length(self):
        L = 0.0
        for seg in self.segments: L += seg.length()
        return L

    # end class Polyline

class Polyline_old(Path):
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

class Spline(Polyline):
    """
    Construct a spline of Bezier segments from a sequence of points.
    """
    def __init__(self, points, closed=False, tolerance=1.0e-10):
        self.points = [Vector3(p) for p in points]
        if closed and (abs(self.points[0]-self.points[-1]) > tolerance):
            self.points.append(Vector3(points[0]))
        self.closed = closed

        # Given m+1 interpolation points p, determine the m-segment
        # Bezier polyline that interpolates these points as a spline.
        # This is done by first determining the array of weight points
        # which define the spline and then evaluating the cubic
        # Bezier segments.
        # Reference:
        #   G. Engelin & F. Uhlig (1996)
        #   Numerical Algorithms with C
        #   Springer, Berlin
        #   Section 12.3.1
        #
        # For a natural spline, the first and last weight points
        # are also the first and last interpolation points.
        # And, for the initial guess at the remaining weight points,
        # just use the supplied data points.
        # This amounts to copying the whole p collection.
        m = len(self.points)-1
        d = [Vector3(p) for p in self.points]
        #
        # Apply Gauss-Seidel iteration until
        # the internal weight points converge.
        for j in range(50):
            max_diff = 0.0
            for i in range(1, m):
                old_p = Vector3(d[i])
                d[i] = 0.25*(6.0*self.points[i] - d[i-1] - d[i+1])
                max_diff = max(max_diff, abs(d[i]-old_p))
            if max_diff < tolerance: break
        #
        # Final stage; calculate the cubic Bezier segments.
        segments = [Bezier([self.points[i],
                            (2.0*d[i]+d[i+1])/3.0,
                            (d[i]+2.0*d[i+1])/3.0,
                            self.points[i+1]])
                    for i in range(m)]
        super().__init__(segments, closed, tolerance)
        return

    def __repr__(self):
        return f"Spline(points={self.points}, closed={self.closed})"


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
        self.t_values = np.array(self.t_values)
        self.arc_lengths = np.array(self.arc_lengths)
        return

    def underlying_t(self, t):
        """
        Search the pieces of arc length to find the piece containing the
        desired point and then interpolate the local value of t for that piece.
        """
        # The incoming parameter value, t, is proportional to arc_length fraction.
        L_target = t * self.arc_lengths[-1]

        # Do a single variable linear interpolation to approximate an ordinary t value
        ut = np.interp(L_target, self.arc_lengths, self.t_values, left=0.0, right=1.0)
        return ut

    def __repr__(self):
        return "ArcLengthParameterizedPath(underlying_path={}, n={})".format(
            self.underlying_path, self._n)

    def __call__(self, t):
        return self.underlying_path(self.underlying_t(t))

    def length(self):
        return self.underlying_path.length()

    # end class ArcLengthParameterizedPath

class ArcLengthParameterizedPath_old(Path):
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

if __name__=='__main__':
    # Test code by NNG
    p0 = Vector3(0.0, 0.0)
    p1 = Vector3(1.0, 1.0)
    line = Line(p0, p1)
    x = line(0.5)
    xx = line(np.array([0.5, 0.5]))
    assert(np.isclose(x.x, xx.x[0]))

    # FnPath PJ 2022-11-14
    def myFun(t):
        end0 = Vector3(0.0, 0.0, 0.0)
        end1 = Vector3(1.0, 2.0, 3.0)
        return end0*(1.0-t) + end1*t
    funpth = FnPath(myFun)
    x = funpth(0.5)
    xx = funpth(np.array([0.5, 0.5]))
    assert(np.isclose(x.x, 0.5))
    assert(np.isclose(xx.x[0], 0.5))

    a = Vector3(1.0, 0.0)
    b = Vector3(0.0, 1.0)
    c = Vector3(0.0, 0.0)
    arc = Arc(a,b,c)
    x = arc(0.5)
    xx = arc(np.array([0.5, 0.5]))
    assert(np.isclose(x.x, xx.x[0]))

    a = Vector3(0.0, 0.0)
    b = Vector3(0.25, 0.25)
    c = Vector3(1.0, 1.0)
    bez = Bezier([a,b,c])
    x = bez(0.5)
    xx = bez(np.array([0.5, 0.5]))
    assert(np.isclose(x.x, xx.x[0]))

    a = Vector3(0.0, 0.0)
    b = Vector3(1.0, 1.0)
    c = Vector3(4.0, 4.0)
    l0 = Line(a, b)
    l1 = Line(b, c)
    polyline = Polyline([l0, l1])
    polyline2 = Polyline_old([l0, l1])
    assert(np.isclose(polyline(0.0).x,      polyline2(0.0).x))
    assert(np.isclose(polyline(0.25/2.0).x, polyline2(0.25/2.0).x))
    assert(np.isclose(polyline(0.75).x,     polyline2(0.75).x))
    assert(np.isclose(polyline(1.0).x,      polyline2(1.0).x))

    xx = polyline(np.array([0.0, 0.25/2.0, 0.75, 1.0]))
    assert(np.isclose(polyline(0.0).x,      xx.x[0]))
    assert(np.isclose(polyline(0.25/2.0).x, xx.x[1]))
    assert(np.isclose(polyline(0.75).x,     xx.x[2]))
    assert(np.isclose(polyline(1.0).x,      xx.x[3]))

    a = [0.0, 0.0]
    b = [0.25, 0.25]
    c = [1.0, 1.0]
    spline = Spline([a,b,c])
    x = spline(0.5)
    xx = spline(np.array([0.5, 0.5]))
    assert(np.isclose(x.x, xx.x[0]))

    a = Vector3(0.0, 0.0)
    b = Vector3(1.0, 1.0)
    c = Vector3(4.0, 4.0)
    l0 = Line(a, b)
    l1 = Line(b, c)
    polyline = Polyline([l0, l1])
    ppath = ArcLengthParameterizedPath(polyline)
    ppath2 = ArcLengthParameterizedPath_old(polyline)
    x = ppath(0.5)
    x2 = ppath2(0.5)
    assert(np.isclose(x.x, x2.x))
    assert(np.isclose(x.y, x2.y))

    x = ppath(0.0)
    x2 = ppath2(0.0)
    assert(np.isclose(x.x, x2.x))
    assert(np.isclose(x.y, x2.y))

    x = ppath(1.0)
    x2 = ppath2(1.0)
    assert(np.isclose(x.x, x2.x))
    assert(np.isclose(x.y, x2.y))

    x = ppath(0.5)
    xx = ppath(np.array([0.5, 0.5]))
    assert(np.isclose(x.x, xx.x[0]))
    assert(np.isclose(x.y, xx.y[0]))
