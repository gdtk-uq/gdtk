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
    __slots__ = ['xs', 'fs']

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
        else:
            # We start with no segments.
            nseg = 0
            xs = []
            fs = []
        return

    def __repr__(self):
        return "XPath(xs={}, fs={})".format(self.xs, self.fs)

    def __call__(self, x):
        nseg = len(self.fs)
        if nseg == 1: return self.fs[0](x)
        i = 0
        while (i < nseg) and (x > self.xs[i]): i += 1
        i = max(0, min(i-1, nseg-1))
        return self.fs[i](x)

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
