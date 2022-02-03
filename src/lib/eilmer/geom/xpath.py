# xpath.py
"""
Class to collect functions to act as segments for a path (x, y(x)).

Peter J.
2022-02-01 : initially built for the Puffin flow calculator.
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
        Note that outside the range is xs values, we extrapolate
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

    # end class XPath

if __name__ == '__main__':
    print("Try out XPath")
    def f0(x): return -1.0
    def f1(x): return 2.0
    xp = XPath([0.0, 1.0, 2.0], [f0, f1])
    xtest = [-0.1, 0.0, 0.1, 0.2, 1.3, 2.1]
    ytest = [xp(x) for x in xtest]
    print("xtest=", xtest)
    print("ytest=", ytest)
    #
    xp.moveto(0.0, 3.0).lineto(1.0, 3.5).lineto(2.0, 4.0)
    ytest = [xp(x) for x in xtest]
    print("xtest=", xtest)
    print("ytest=", ytest)
