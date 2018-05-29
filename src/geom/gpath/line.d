// line.d
// Peter J. 2017-11-29: Split out of the original path module.

module geom.gpath.line;

import std.math;
import std.conv: to;
import nm.complex;
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
        return to!number(1.0-t)*p0 + to!number(t)*p1;
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
    override double partial_length(double ta, double tb) const
    {
        Vector3 dp = p1 - p0;
        return fabs(tb - ta) * geom.abs(dp).re;
    }
    override Vector3 point_from_length(double length, out double t) const
    {
        Vector3 dp = p1 - p0;
        t = length/(geom.abs(dp).re);
        return this.opCall(t);
    }
} // end class Line


version(line_test) {
    import util.msg_service;
    int main() {
        auto a = Vector3([1.0, 2.2, 3.0]);
        auto b = Vector3(1.0);
        auto ab = new Line(a, b);
        auto c = ab(0.5);
        assert(approxEqualVectors(c, Vector3(1.0, 1.1, 1.5)), failedUnitTest());
        auto ab2 = ab.dup();
        auto d = ab2(0.5);
        assert(approxEqualVectors(c, d), failedUnitTest());
        auto pth = new Line(Vector3(0.0,1.0), Vector3(1.0,1.0));
        auto ps = Vector3(0.5,0.5);
        auto dir = Vector3(0.0,1.0);
        double t;
        auto found = pth.intersect2D(ps, dir, t);
        assert(found, failedUnitTest()); // "intersect2D not found on Line"
        // intersect2D parametric location on Line
        assert(approxEqual(t,0.5), failedUnitTest());
        return 0;
    }
} // end line_test
