// xspline.d
// A Path based on the cubic-spline y(x).
//
// Peter J.
// 2022-04-04 : first code, making use of nm.spline and nm.splinelsq.
//
module geom.gpath.xspline;

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

class XSpline : Path {
    // A Path based the cubic-spline y(x).

public:
    CubicSpline spl;

    this(double[] xs, double[] ys)
    // Constructs a spline through specified (x,y) points.
    {
        spl = new CubicSpline(xs, ys);
    }

    this(string fileName)
    // Contructs a spline from a file containing x(,y) coordinates.
    {
        // This function takes a filename and processes it assuming that each
        // line contains (x,y) tuples (space-delimited).  If any y-values are
        // missing on a given line, they are assumed to be 0.0.
        double[] xs, ys;
        auto f = File(fileName, "r");
        foreach (line; f.byLine) {
            auto tokens = line.strip().split();
            if (tokens.length == 0) continue; // ignore blank lines
            if (tokens[0] == "#") continue; // ignore comment lines
            double x = to!double(tokens[0]);
            double y = 0.0; if (tokens.length > 1) { y = to!double(tokens[1]); }
            xs ~= x; ys ~= y;
        }
        spl = new CubicSpline(xs, ys);
    }

    this(const(XSpline) other)
    {
        spl = new CubicSpline(other.spl);
    }

    override XSpline dup() const
    {
        return new XSpline(this);
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
        return "XSpline(xs=" ~ to!string(spl.xs) ~ ", ys=" ~ to!string(spl.ys) ~ ")";
    }

    override string classString() const
    {
        return "XSpline";
    }

} // end class XSpline

version(xspline_test) {
    import util.msg_service;
    import std.stdio: writefln;
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
        auto s = new XSpline(x_sample, y_sample);
        N = 100;
        double dt = 1.0/(N-1);
        foreach (i; 0 .. N) {
            double t = dt * i;
            double xx = s(0).x * (1.0-t) + s(1).x * t;
            double y_runge = runge(xx);
            Vector3 p = s(t);
            // writefln("%g %g %g", xx, y_runge, p.y);
            assert(approxEqualVectors(p, Vector3(xx, y_runge), 0.02), failedUnitTest());
        }
        return 0;
    }
} // end xspline_test

