// xsplinelsq.d
// A Path based on the cubic-spline y(x) that has been fitted to data points.
//
// Peter J.
// 2022-04-04 : first code, making use of nm.spline and nm.splinelsq.
//
module geom.gpath.xsplinelsq;

import std.conv;
import std.math;
import std.stdio: File;
import std.string;

import ntypes.complex;
import nm.number;
import nm.spline;
import nm.splinelsq;

import geom.elements;
import geom.gpath.path;

class XSplineLsq : Path {
    // A Path based the cubic-spline y(x).

public:
    CubicSplineLsq spl;

    this(double[] xd, double[] yd, double[] wd, double[] xs, int nseg)
    // Constructs a spline fitted in a least-squares error sense
    // specified (x,y) data points.
    {
        spl = new CubicSplineLsq(xd, yd, wd, xs, nseg);
    }

    this(string fileName, double[] xs, int nseg)
    // Contructs a spline from a file containing x(,y(,w)) coordinates.
    {
        // This function takes a filename and processes it assuming that each
        // line contains (x,y,w) tuples (space-delimited).
        // If any y-values are missing on a given line, they are assumed to be 0.0.
        // If any weights are missing, they are assumed to be 1.0.
        double[] xd, yd, wd;
        auto f = File(fileName, "r");
        foreach (line; f.byLine) {
            auto tokens = line.strip().split();
            if (tokens.length == 0) continue; // ignore blank lines
            if (tokens[0] == "#") continue; // ignore comment lines
            double x = to!double(tokens[0]);
            double y = 0.0; if (tokens.length > 1) { y = to!double(tokens[1]); }
            double w = 1.0; if (tokens.length > 2) { w = to!double(tokens[2]); }
            xd ~= x; yd ~= y; wd ~= w;
        }
        assert (nseg > 2, "Too few segments specified.");
        if (xs.length > 0) {
            // Already have xs locations.
            assert(xs.length == nseg+1, "Incorrect number of knots specified.");
        } else {
            // Distribute xs locations uniformly over data range.
            xs.length = nseg+1;
            foreach (i; 0 .. nseg+1) {
                double t = to!double(i)/nseg;
                xs[i] = xd[0]*(1.0-t) + xd[$-1]*t;
            }
        }
        spl = new CubicSplineLsq(xd, yd, wd, xs, nseg);
    } // end constructor

    this(const(XSplineLsq) other)
    {
        spl = new CubicSplineLsq(other.spl);
    }

    override XSplineLsq dup() const
    {
        return new XSplineLsq(this);
    }

    override Vector3 opCall(double t) const
    {
        double x = spl.xmin() * (1.0-t) + spl.xmax() * t;
        double y = spl(x);
        double z = 0.0;
        return Vector3(x, y, z);
    } // end opCall()

    override string toString() const
    {
        return "XSplineLsq(xd=" ~ to!string(spl.xd) ~
            ", yd=" ~ to!string(spl.yd) ~
            ", xs=" ~ to!string(spl.xs) ~
            ", nseg=" ~ to!string(spl.nseg) ~
            ")";
    }

    override string classString() const
    {
        return "XSplineLsq";
    }

} // end class XSplineLsq

version(xsplinelsq_test) {
    import util.msg_service;
    import std.stdio: writeln, writefln;
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
        auto s = new XSplineLsq(x_sample, y_sample, w_sample, xs, 10);
        N = 100;
        double max_dy = 0.0;
        double dt = 1.0/(N-1);
        foreach (i; 0 .. N) {
            double t = dt * i;
            Vector3 p = s(t);
            double y_runge = runge(p.x);
            double dy = p.y - y_runge;
            max_dy = fmax(max_dy, fabs(dy));
            // writefln("%g %g %g %g", p.x, y_runge, p.y, dy);
            assert(fabs(dy) < 0.02, failedUnitTest());
        }
        // writeln("max_dy=", max_dy);
        return 0;
    }
} // end xsplinelsq_test

