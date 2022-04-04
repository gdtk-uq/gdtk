// spline.d
// Simple implementation of interpolation with a natural cubic spline.
//
// Peter J. 20-Oct-2011, to go with the mech3750 lecture.
//          16-Aug-2014, Python3 port
//          01-Feb-2022, Small edits when including in the Eilmer library.
//          04-Apr-2022, Dlang variant.
//

module nm.spline;
import nm.bbla: Matrix, zeros, gaussJordanElimination;
import std.math: fabs;
import std.stdio: writeln, writefln;


class CubicSpline {
    // Interpolatory cubic spline, storing coefficients for the polynomial pieces.

public:
    double[] xs, ys; // Interpolation points of the spline.

    this(const double[] xs, const double[] ys)
    // Sets up the interpolatory cubic spline through the xi, yi points.
    //   xs : sequence of x-coordinates
    //   ys : sequence of y-coordinates
    //
    // There will be n+1 knots, n segments and n-1 interior knots.
    {
        n = xs.length - 1;  // number of segments
        assert(ys.length == n+1, "Unequal xs, ys lengths.");
        this.xs = xs.dup();
        this.ys = ys.dup();
        a.length = n;
        b.length = n;
        c.length = n;
        // for d, just use the y array
        double[] h; h.length = n;
        foreach (i; 0 .. n) { h[i] = xs[i+1] - xs[i]; }
        // Solve for the second-derivative at the interior knots.
        // The sigma[i], i=1...n-1 (inclusive) are the unknowns since
        // the natural end conditions are sigma[0]=sigma[n]=0.0
        size_t m = n-1; // number of unknowns (and constraint equations)
        auto A = zeros!double(m, m+1); // Augmented matrix
        foreach (k; 0 .. m) {
            size_t i = k + 1; // i-index as in the lecture notes
            if (k > 0) { A[k,k-1] = h[i-1]; }
            A[k,k] = 2.0*(h[i-1]+h[i]);
            if (k < m-1) { A[k,k+1] = h[i]; }
            A[k,m] = 6.0*((ys[i+1]-ys[i])/h[i] - (ys[i]-ys[i-1])/h[i-1]); // RHS elements
        }
        gaussJordanElimination!double(A);
        // Put sigma values in place in the full array.
        double[] sigma; foreach (i; 0 .. n+1) { sigma ~= 0.0; }
        foreach (i; 0 .. m) { sigma[i+1] = A[i,m]; }
        // Evaluate and store coefficients for the polynomial segments.
        foreach (i; 0 .. n) {
            a[i] = (sigma[i+1]-sigma[i])/(6.0*h[i]);
            b[i] = sigma[i]/2.0;
            c[i] = (ys[i+1]-ys[i])/h[i] - (2.0*sigma[i] + sigma[i+1])*h[i]/6.0;
        }
    } // end constructor

    this(ref const(CubicSpline) other)
    {
        this(other.xs, other.ys);
    }

    double opCall(double x) const
    // Evaluates the spline at point x by first searching for the
    // relevant segment and then evaluating the local polynomial piece.
    {
        int i = 0;
        if (x < xs[0]) {
            i = 0; // The point is off range to the left.
        } else {
            i = cast(int)n - 1; // Start search from the right-hand end.
            while ((x < xs[i]) && (i > 0)) { i -= 1; }
        }
        // Now that we found the relevant segment, evaluate locally.
        double dx = x - xs[i];
        return ((a[i]*dx + b[i])*dx + c[i])*dx + ys[i];
    }

    double xmin() const { return xs[0]; }

    double xmax() const { return xs[$-1]; }

private:
    size_t n;
    double[] a, b, c;
} // end class CubicSpline


version(spline_test) {
    import util.msg_service;
    import std.conv;
    import nm.number;
    int main() {
        double runge(double x) { return 1.0/(1.0 + 25* x * x); }
        int N = 20;
        double x0 = -1.0;
        double x1 = 1.0;
        double dx = (x1-x0)/(N-1);
        double[] x_sample, y_sample;
        foreach (i; 0 .. N) {
            double xx = x0 + dx*i;
            x_sample ~= xx;
            y_sample ~= runge(xx);
        }
        auto s = new CubicSpline(x_sample, y_sample);
        N = 100;
        dx = (x1-x0)/(N-1);
        foreach (i; 0 .. N) {
            double xx = x0 + dx*i;
            double y_runge = runge(xx);
            double y_spline = s(xx);
            double dy = y_spline - y_runge;
            // writefln("%g %g %g %g", xx, y_runge, y_spline, dy);
            assert(fabs(dy) < 0.02, failedUnitTest());
        }
        return 0;
    } // end main()
}
