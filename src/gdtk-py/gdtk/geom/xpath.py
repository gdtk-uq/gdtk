# xpath.py
"""
Class to collect functions to act as segments for a path (x, y(x)).

Peter J.
2022-02-01 : initially built for the Puffin flow calculator.
2022-02-03 : moveto and lineto methods
2022-02-10 : bezier2to and bezier3to methods
2022-02-11 : XBezier
"""

from copy import copy


class XPath(object):
    """
    Segmented function y(x).
    """
    __slots__ = ['xs', 'ys', 'fs']

    def __init__(self, xs=None, fs=None):
        """
        We may construct an XPath as a set of x values and a set of callable functions.
        Input:
        xs : sequence of n+1 x-values that define the extents of the segments.
        fs : sequence of n functions defining the y-coordinate y(x) for each segment.
        """
        if xs and fs:
            # We have been given information to start buiding the path.
            nseg = len(fs)
            assert len(xs) >= nseg+1, "Too few x values."
            self.xs = [float(x) for x in xs]
            assert all([callable(f) for f in fs]), "Not all f items are callable."
            self.fs = copy(fs)
            ys = [fs[0](self.xs[0]),]
            for i in range(nseg): ys.append(self.fs[i](self.xs[i+1]))
        else:
            # We start with no segments.
            nseg = 0
            xs = []
            fs = []
        return

    def __repr__(self):
        return "XPath(xs={}, fs={})".format(self.xs, self.fs)

    def __call__(self, x):
        """
        Returns y(x).
        Note that outside the range of xs values, we extrapolate
        with the function for the nearest available segment.
        """
        nseg = len(self.fs)
        if nseg == 1: return self.fs[0](x)
        i = 0
        while (i < nseg) and (x > self.xs[i]): i += 1
        i = max(0, min(i-1, nseg-1))
        return self.fs[i](x)

    def moveto(self, x, y):
        """
        Starts an XPath at the given coordinates.
        Returns the object so that we may chain methods.
        """
        self.xs = [x,]
        self.ys = [y,]
        self.fs = []
        return self

    def lineto(self, x, y):
        """
        Appends a linear segment, given the coordinates of end point of segment.
        Returns the object so that we may chain methods.
        """
        assert len(self.xs) > 0, "No starting point to which we can append."
        x0 = self.xs[-1]; y0 = self.ys[-1]
        x1 = x; y1 = y
        def fn(x):
            frac = (x-x0)/(x1-x0)
            return y0*(1.0-frac) + y1*frac
        self.xs.append(x1)
        self.ys.append(y1)
        self.fs.append(fn)
        return self

    def bezier2to(self, x1, y1, x2, y2):
        """
        Appends a quadratic Bezier segment, given the Bezier control points.
        Returns the object so that we may chain methods.
        """
        assert len(self.xs) > 0, "No starting point to which we can append."
        x0 = self.xs[-1]; y0 = self.ys[-1]
        def xbez(t): return x0*(1-t)**2 + 2*x1*(1-t)*t + x2*t**2
        def dxdt(t): return 2*((x1-x0)*(1-t) + (x2-x1)*t)
        def ybez(t): return y0*(1-t)**2 + 2*y1*(1-t)*t + y2*t**2
        def fn(x):
            # Initial guess for t assumes linear distribution.
            t = (x-x0)/(x2-x0)
            def g(t): return xbez(t)-x
            dt = -g(t)/dxdt(t) # Newton-Raphson increment
            count = 0
            while abs(dt) > 1.0e-11 and count < 20:
                t += dt
                dt = -g(t)/dxdt(t)
                count += 1
            # At this point we should have t that corresponds to x.
            return ybez(t)
        self.xs.append(x2)
        self.ys.append(y2)
        self.fs.append(fn)
        return self

    def bezier3to(self, x1, y1, x2, y2, x3, y3):
        """
        Appends a cubic Bezier segment, given the Bezier control points.
        Returns the object so that we may chain methods.
        """
        assert len(self.xs) > 0, "No starting point to which we can append."
        x0 = self.xs[-1]; y0 = self.ys[-1]
        def xbez(t): return x0*(1-t)**3 + 3*x1*(1-t)**2*t + 3*x2*(1-t)*t**2 + x3*t**3
        def dxdt(t): return 3*((x1-x0)*(1-t)**2 + 2*(x2-x1)*(1-t)*t + (x3-x2)*t**2)
        def ybez(t): return y0*(1-t)**3 + 3*y1*(1-t)**2*t + 3*y2*(1-t)*t**2 + y3*t**3
        def fn(x):
            # Initial guess for t assumes linear distribution.
            t = (x-x0)/(x3-x0)
            def g(t): return xbez(t)-x
            dt = -g(t)/dxdt(t) # Newton-Raphson increment
            count = 0
            while abs(dt) > 1.0e-11 and count < 20:
                t += dt
                dt = -g(t)/dxdt(t)
                count += 1
            # At this point we should have t that corresponds to x.
            return ybez(t)
        self.xs.append(x3)
        self.ys.append(y3)
        self.fs.append(fn)
        return self

    # end class XPath


class XBezier(object):
    """
    Computes y(x) based on parametric Bezier polynomials x=x(t) y=y(t).
    """
    __slots__ = ['xbs', 'ybs', 'dxbs']

    def __init__(self, xbs=None, ybs=None):
        """
        Construct from a sequence of x-coordinates and a sequence of y-coordinates
        for the Bezier control points.

        We assume that the x coordinates are in order of increasing value
        and close to being equally distributed.
        """
        assert len(xbs) > 1, "Too few defining points."
        assert len(xbs) == len(ybs), "Unequal sequences of coordinates."
        self.xbs = [float(x) for x in xbs]
        self.ybs = [float(y) for y in ybs]
        # Retain differences for later calculation of the derivative dxdt(x).
        self.dxbs = [self.xbs[i+1]-self.xbs[i] for i in range(len(self.xbs)-1)]
        return

    def __repr__(self):
        return "XBezier(xbs={}, ybs={})".format(self.xbs, self.ybs)

    def bez(self, b_array, t):
        """
        Generic Bezier evaluation given a specific array of control points.
        """
        B = [b for b in b_array]
        while len(B) > 1:
            B = [B[i]*(1-t) + B[i+1]*t for i in range(len(B)-1)]
        return B[0]

    def xbez(self, t): return self.bez(self.xbs, t)

    def ybez(self, t): return self.bez(self.ybs, t)

    def dxdt(self, t): return len(self.dxbs)*self.bez(self.dxbs, t)

    def __call__(self, x):
        """
        Returns y(x).

        Note that, for x outside the range of xbs values,
        we return the nearest end-point value.
        """
        if x <= self.xbs[0]: return self.ybs[0]
        if x >= self.xbs[-1]: return self.ybs[-1]
        # Initial guess for t assumes linear distribution.
        t = (x-self.xbs[0])/(self.xbs[-1]-self.xbs[0])
        def g(t): return self.xbez(t)-x
        dt = -g(t)/self.dxdt(t) # Newton-Raphson increment
        count = 0
        while abs(dt) > 1.0e-11 and count < 20:
            t += dt
            dt = -g(t)/self.dxdt(t)
            count += 1
        # At this point we should have t that corresponds to x.
        return self.ybez(t)

    # end class XBezier


if __name__ == '__main__':
    print("Try out XPath class.")
    def f0(x): return -1.0
    def f1(x): return 2.0
    xp = XPath([0.0, 1.0, 2.0], [f0, f1])
    xtest = [-0.1, 0.0, 0.1, 0.2, 1.3, 2.1]
    ytest = [xp(x) for x in xtest]
    print("xtest=", xtest)
    print("ytest=", ytest)
    #
    print("Exercise moveto, lineto.")
    xp.moveto(0.0, 3.0).lineto(1.0, 3.5).lineto(2.0, 4.0)
    ytest = [xp(x) for x in xtest]
    print("xtest=", xtest)
    print("ytest=", ytest)
    #
    print("bezier2to: Equally-spaced x locations.")
    xp2 = XPath().moveto(0.0,0.0).bezier2to(0.5,1.0, 1.0,0.0)
    xtest = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]
    ytest = [xp2(x) for x in xtest]
    print("xtest=", xtest)
    print("ytest=", ytest)
    #
    print("bezier3to: Unequally-spaced x locations")
    xp3 = XPath().moveto(0.0,0.0).bezier3to(1.0/4,1.0, 3.0/4,-1.0, 1.0,0.0)
    xtest = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]
    ytest = [xp3(x) for x in xtest]
    print("xtest=", xtest)
    print("ytest=", ytest)
    #
    print("XBezier: Unequally-spaced x locations.")
    xp4 = XBezier([0.0,1.0/4,3.0/4,1.0],[0.0,1.0,-1.0,0.0])
    xtest = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]
    ytest = [xp4(x) for x in xtest]
    print("xtest=", xtest)
    print("ytest=", ytest)
