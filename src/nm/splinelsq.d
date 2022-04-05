// splinelsq.d
// Brute-force implementation of a least-squares fitted spline.
//
// Peter J.
// 2022-04-04: First code.
//

module nm.splinelsq;
import std.math: fabs, fmax;
import std.stdio: writeln, writefln;
import std.algorithm.searching: minElement, maxElement;
import std.conv: to;
import nm.spline;
import nm.nelmin;

class CubicSplineLsq {
    // A CubicSpline that has been fitted to the original data in a least-squares error sense.

public:
    double[] xd, yd; // Copy of the original data.
    double[] wd; // Weights for the data points.
    double[] xs; // x-coordinates of the spline knots
    int nseg; // Number of segments in the underlying spline.

    this(const double[] xd, const double[] yd, const double[] wd, const double[] xs, int nseg)
    // Construct a spline of nseg segments to approximate the xd, yd points.
    // wd is the array of weights for the data points,
    // for use when computing the sum-square error.
    // xs array may be used to set the x-locations of the spline knots.
    // If it is supplied, it's length has to be nseg+1.
    {
        // Check for reasonable arrays coming in.
        this.xd = xd.dup();
        size_t nd = xd.length;
        assert(nd > nseg+1, "Too few data points.");
        this.yd = yd.dup();
        assert(this.yd.length == nd, "yd array not same length as xd.");
        this.wd = wd.dup();
        if (this.wd.length == 0) {
            this.wd.length = nd;
            foreach (i; 0 .. nd) { this.wd[i] = 1.0; }
        } else {
            assert(this.wd.length == nd, "Inconsistent length for wd.");
        }
        if (xs.length > 0) {
            assert(xs.length == nseg+1, "Inconsistent length for xs.");
        }
        this.xs = xs.dup();
        this.nseg = nseg;
        assert(this.nseg > 1, "Too few segments.");
        //
        // Set up an initial guess for the spline as a straight line.
        double[] ys;
        if (this.xs.length == 0) {
            foreach (i; 0 .. nseg+1) {
                double frac = to!double(i)/nseg;
                this.xs ~= xd[0]*(1.0-frac) + xd[$-1]*frac;
                ys ~= yd[0]*(1.0-frac) + yd[$-1]*frac;
            }
        } else {
            foreach (i; 0 .. nseg+1) {
                double frac = (this.xs[i]-xd[0])/(xd[$-1]-xd[0]);
                ys ~= yd[0]*(1.0-frac) + yd[$-1]*frac;
            }
        }
        // Use the range of the y data to estimate a suitable perturbation.
        double dy = maxElement(yd) - minElement(yd); dy *= 0.1;
        double[] delys; foreach (i; 0 .. nseg+1) { delys ~= dy; }
        // Set up the residual function.
        double sse(double[] ps) {
            // The parameter vector is the vector of y-values for the spline knots.
            assert(ps.length == this.xs.length, "Unexpected length for parameter vector.");
            spl = new CubicSpline(this.xs, ps);
            double sumerr = 0.0;
            foreach (j; 0 .. xd.length) {
                double e = yd[j] - spl(xd[j]);
                sumerr += this.wd[j]*e*e;
            }
            return sumerr;
        }
        // Optimize the spline by adjusting the ys values.
        double f;
        int nfe, nres;
        double tol = 1.0e-6;
        int P = 2;
        int nfe_max = 1200;
        bool conv_flag = minimize!(sse,double)(ys, f, nfe, nres, delys, tol, P, nfe_max);
        // writeln("ys= ", ys, " f= ", f);
        // writeln("convergence-flag= ", conv_flag);
        // writeln("number-of-fn-evaluations= ", nfe);
        // writeln("number-of-restarts= ", nres);
        spl = new CubicSpline(this.xs, ys);
    }

    this(const(CubicSplineLsq) other)
    {
        spl = new CubicSpline(other.spl);
        nseg = other.nseg;
        xs = other.xs.dup();
        xd = other.xd.dup();
        yd = other.yd.dup();
        wd = other.wd.dup();
    }

    double opCall(double x) const
    // Evaluates the underlying spline at point x.
    {
        return spl(x);
    }

    double xmin() const
    {
        return spl.xmin();
    }

    double xmax() const
    {
        return spl.xmax();
    }

private:
    CubicSpline spl; // The underlying spline model.
} // end class CubicSplineLsq


version(splinelsq_test) {
    import util.msg_service;
    import std.conv;
    import nm.number;
    int main() {
        double runge(double x) { return 1.0/(1.0 + 25* x * x); }
        int N = 200;
        double x0 = -1.0;
        double x1 = 1.0;
        double dx = (x1-x0)/(N-1);
        double[] x_sample, y_sample, w_sample;
        foreach (i; 0 .. N) {
            double xx = x0 + dx*i;
            x_sample ~= xx;
            y_sample ~= runge(xx);
            w_sample ~= 1.0;
        }
        w_sample[0] = 100.0; w_sample[$-1]=100.0;
        double[] xs;
        auto s = new CubicSplineLsq(x_sample, y_sample, w_sample, xs, 10);
        N = 100;
        dx = (x1-x0)/(N-1);
        double max_dy = 0.0;
        foreach (i; 0 .. N) {
            double xx = x0 + dx*i;
            double y_runge = runge(xx);
            double y_spline = s(xx);
            double dy = y_spline - y_runge;
            max_dy = fmax(max_dy, fabs(dy));
            // writefln("%g %g %g %g", xx, y_runge, y_spline, dy);
            assert(fabs(dy) < 0.02, failedUnitTest());
        }
        // writeln("max_dy=", max_dy);
        return 0;
    } // end main()
}
