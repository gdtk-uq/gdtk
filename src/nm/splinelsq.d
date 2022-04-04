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
    this(double[] xd, double[] yd, double[] wd, double[] xs, int nseg)
    // Construct a spline of nseg segments to approximate the xd, yd points.
    // wd is the array of weights for the data points,
    // for use when computing the sum-square error.
    // xs array may be used to set the x-locations of the spline knots.
    // If it is supplied, it's length has to be nseg+1.
    {
        // Check for reasonable arrays coming in.
        size_t nd = xd.length;
        assert(nd > nseg+1, "Too few data points.");
        assert(yd.length == nd, "yd array not same length as xd.");
        if (wd.length == 0) {
            wd.length = nd;
            foreach (i; 0 .. nd) { wd[i] = 1.0; }
        }
        assert(wd.length == nd, "Inconsistent length for wd.");
        if (xs.length > 0) {
            assert(xs.length == nseg+1, "Inconsistent length for xs.");
        }
        assert(nseg > 1, "Too few segments.");
        //
        // Set up an initial guess for the spline as a straight line.
        double[] ys;
        if (xs.length == 0) {
            foreach (i; 0 .. nseg+1) {
                double frac = to!double(i)/nseg;
                xs ~= xd[0]*(1.0-frac) + xd[$-1]*frac;
                ys ~= yd[0]*(1.0-frac) + yd[$-1]*frac;
            }
        } else {
            foreach (i; 0 .. nseg+1) {
                double frac = (xs[i]-xd[0])/(xd[$-1]-xd[0]);
                ys ~= yd[0]*(1.0-frac) + yd[$-1]*frac;
            }
        }
        // Use the range of the y data to estimate a suitable perturbation.
        double dy = maxElement(yd) - minElement(yd); dy *= 0.1;
        double[] delys; foreach (i; 0 .. nseg+1) { delys ~= dy; }
        // Set up the residual function.
        double sse(double[] ps) {
            // The parameter vector is the vector of y-values for the spline knots.
            assert(ps.length == xs.length, "Unexpected length for parameter vector.");
            spl = new CubicSpline(xs, ps);
            double sumerr = 0.0;
            foreach (j; 0 .. xd.length) {
                double e = yd[j] - spl(xd[j]);
                sumerr += wd[j]*e*e;
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
        spl = new CubicSpline(xs, ys);
    }

    double opCall(double x)
    // Evaluates the underlying spline at point x.
    {
        return spl(x);
    }

private:
    CubicSpline spl; // The underlying spline model.
    int nseg; // Number of segments in the underlying spline.
    double[] xd, yd; // Copy of the original data.
    double wd; // Weights for the data points.
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
