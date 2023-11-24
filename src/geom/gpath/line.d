// line.d
// Peter J. 2017-11-29: Split out of the original path module.

module geom.gpath.line;

import std.math;
import std.conv: to;
import ntypes.complex;
import nm.number;

import geom.elements;
import geom.gpath.path;


class Line : Path {
public:
    Vector3 p0; // end-point at t=0
    Vector3 p1; // end-point at t=1
    this(in Vector3 p0, in Vector3 p1)
    {
        this.p0 = p0; this.p1 = p1;
    }
    this(ref const(Line) other)
    {
        p0 = other.p0; p1 = other.p1;
    }
    override Line dup() const
    {
        return new Line(p0, p1);
    }
    override Vector3 opCall(double t) const
    {
        return (1.0-t)*p0 + t*p1;
    }
    override Vector3 dpdt(double t) const
    {
        return p1 - p0;
    }
    override Vector3 d2pdt2(double t) const
    {
        return Vector3(0.0,0.0,0.0);
    }
    override string toString() const
    {
        return "Line(p0=" ~ to!string(p0) ~ ", p1=" ~ to!string(p1) ~ ")";
    }
    override string classString() const
    {
        return "Line";
    }
    override number partial_length(double ta, double tb) const
    {
        Vector3 dp = p1 - p0;
        return fabs(tb - ta) * geom.abs(dp);
    }
    override Vector3 point_from_length(number length, out double t) const
    {
        Vector3 dp = p1 - p0;
        t = (length/geom.abs(dp)).re;
        return this.opCall(t);
    }
} // end class Line


version(line_test) {
    import util.msg_service;
    int main() {
        // Some simple evaluations.
        auto a = Vector3([1.0, 2.2, 3.0]);
        auto b = Vector3(1.0);
        auto ab = new Line(a, b);
        auto c = ab(0.5);
        assert(approxEqualVectors(c, Vector3(1.0, 1.1, 1.5)), failedUnitTest());
        auto ab2 = ab.dup();
        auto d = ab2(0.5);
        assert(approxEqualVectors(c, d), failedUnitTest());
        //
        // Set up an intersection that should be found on Line.
        auto pth = new Line(Vector3(0.0,1.0), Vector3(1.0,1.0));
        auto ps = Vector3(0.5,0.5);
        auto dir = Vector3(0.0,1.0);
        double t;
        auto found = pth.intersect2D(ps, dir, t);
        assert(found, failedUnitTest());
        // intersect2D parametric location on Line
        assert(isClose(t,0.5), failedUnitTest());
        //
        version(complex_numbers) {
            // Try out the complex derivative evaluation.
            auto p0 = Vector3(0.0, 0.0);
            double alpha = 1.0;
            auto p1 = Vector3(alpha, alpha);
            auto line0 = new Line(p0, p1);
            double h = 1.0e-20;
            number ih = complex(0,h);
            auto p1_dash = Vector3(alpha+ih, alpha+ih);
            auto line1 = new Line(p0, p1_dash);
            // What we want to compute is the sensitivity
            // of the midpoint with respect to alpha.
            double dpmid_da_x = line1(0.5).x.im / h;
            double dpmid_da_y = line1(0.5).y.im / h;
            double dpmid_da_z = line1(0.5).z.im / h;
            // import std.stdio;
            // writeln("dpmid_da x:", dpmid_da_x, " y:", dpmid_da_y, " z:", dpmid_da_z);
            assert(isClose(dpmid_da_x,0.5), failedUnitTest());
            assert(isClose(dpmid_da_y,0.5), failedUnitTest());
            assert(isClose(dpmid_da_z,0.0), failedUnitTest());
        }
        //
        return 0;
    }
} // end line_test
