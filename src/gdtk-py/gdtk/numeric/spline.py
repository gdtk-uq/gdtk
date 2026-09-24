# spline.py
"""
Simple implementation of interpolation with a natural cubic spline.

Peter J. 20-Oct-2011, to go with the mech3750 lecture.
         16-Aug-2014, Python3 port
         01-Feb-2022, Small edits when including in the Eilmer library.
"""

from numpy import array, zeros, linalg, linspace

class CubicSpline(object):
    """
    Interpolatory cubic spline, storing coefficients for the polynomial pieces.
    """
    def __init__(self, xi, yi):
        """
        Sets up the interpolatory cubic spline through the xi, yi points.

        xi : sequence of x-coordinates
        yi : sequence of y-coordinates

        There will be n+1 knots, n segments and n-1 interior knots.
        """
        n = len(xi) - 1   # number of segments
        assert len(yi) == n + 1
        self.n = n
        self.x = array(xi)
        self.y = array(yi)
        self.a = zeros(n, float)
        self.b = zeros(n, float)
        self.c = zeros(n, float)
        # for d, just use the y array
        h = self.x[1:n+1] - self.x[0:n]
        # Solve for the second-derivative at the interior knots.
        # The sigma[i], i=1...n-1 (inclusive) are the unknowns since
        # the natural end conditions are sigma[0]=sigma[n]=0.0
        sigma = zeros(n+1, float)
        m = n-1 # number of unknowns (and constraint equations)
        A = zeros((m,m), float)
        rhs = zeros(m, float)
        for k in range(m):
            i = k + 1 # i-index as in the lecture notes
            if k > 0: A[k,k-1] = h[i-1]
            A[k,k] = 2.0*(h[i-1]+h[i])
            if k < m-1: A[k,k+1] = h[i]
            rhs[k] = 6.0*((yi[i+1]-yi[i])/h[i] - (yi[i]-yi[i-1])/h[i-1])
        z = linalg.solve(A, rhs)
        # Put sigma values in place in the full array.
        for i in range(m): sigma[i+1] = z[i]
        # Evaluate and store coefficients for the polynomial segments.
        for i in range(n):
            self.a[i] = (sigma[i+1]-sigma[i])/(6.0*h[i])
            self.b[i] = sigma[i]/2.0
            self.c[i] = (yi[i+1]-yi[i])/h[i] - (2.0*sigma[i] + sigma[i+1])*h[i]/6.0
        return

    def __call__(self, x):
        """
        Evaluates the spline at point x by first searching for the
        relevant segment and then evaluating the local polynomial piece.
        """
        i = 0
        if x < self.x[0]:
            i = 0 # The point is off range to the left.
        else:
            i = self.n-1 # Start search from the right-hand end.
            while x < self.x[i] and i > 0: i -= 1
        # Now that we found the relevant segment, evaluate locally.
        dx = x - self.x[i]
        y = ((self.a[i]*dx + self.b[i])*dx + self.c[i])*dx + self.y[i]
        return y

if __name__ == "__main__":
    print("Begin Cubic_spline demo.")
    def runge(x): return 1.0/(1.0 + 25* x * x)
    N = 10
    x_sample = linspace(-1.0, 1.0, N)
    y_sample = runge(x_sample)
    #
    x_plot = linspace(-1.0, 1.0, 200)
    y_plot = runge(x_plot)
    #
    s = CubicSpline(x_sample, y_sample)
    y_spline = [s(x) for x in x_plot]
    #
    from pylab import plot, xlabel, ylabel, title, show
    plot(x_plot, y_plot, '-',
         x_plot, y_spline, '--',
         x_sample, y_sample, 'o')
    xlabel('x'); ylabel('y')
    title("Cubic spline approximation to Runge's function with N=%d" % N)
    show()
    #
    print("Done.")
