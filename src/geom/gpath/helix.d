// helix.d
// Peter J. 2017-11-29: Split out of the original path module.

module geom.gpath.helix;

import std.conv;
import std.math;
import ntypes.complex;
import nm.number;

import geom.elements;
import geom.gpath.path;


class Helix : Path {
public:
    Vector3 a0; // beginning point on local z-axis (t = 0)
    Vector3 a1; // end point on local z-axis (t = 1)
    Vector3 xdsh; // local x-axis, unit vector
    Vector3 ydsh; // local y-axis, unit vector
    Vector3 zdsh; // local z-axis, unit vector
    double r0, r1; // starting and ending radii
    double theta01; // angle (in radians) from starting point to ending point,
    // assuming the right-hand screw convention.
    //
    // Helix constructed from fundamental parameters.
    this(in Vector3 a0, in Vector3 a1, in Vector3 xlocal,
         double r0, double r1, double dtheta)
    {
        this.a0 = a0; this.a1 = a1;
        this.r0 = r0; this.r1 = r1;
        this.theta01 = dtheta;
        // set up local unit vectors at p0
        xdsh = unit(xlocal);
        zdsh = a1 - a0; zdsh.normalize(); // along the axis of the helix
        ydsh = cross(zdsh, xdsh); // normal to both
    }
    // Helix constructed from point_start to point_end about an axis
    // from axis0 to axis1.
    // We will compute the fundamantal parameters from these points.
    this(in Vector3 point_start, in Vector3 point_end,
         in Vector3 axis0, in Vector3 axis1)
    {
        // Local vectors relative to axis0.
        Vector3 a = axis1 - axis0;
        Vector3 b = point_start - axis0;
        Vector3 c = point_end - axis0;
        zdsh = unit(a);
        Vector3 a0b = b - dot(b,a)*zdsh;
        xdsh = unit(a0b);
        ydsh = cross(zdsh,xdsh);
        a0 = axis0 + dot(b,zdsh)*zdsh;
        a1 = axis0 + dot(c,zdsh)*zdsh;
        r0 = dot(b,xdsh).re;
        Vector3 a1c = c - dot(c,zdsh)*zdsh;
        r1 = geom.abs(a1c).re;
        Vector3 origin = Vector3(0.0, 0.0, 0.0);
        a1c.transform_to_local_frame(xdsh, ydsh, zdsh, origin);
        theta01 = atan2(a1c.y.re, a1c.x.re);
    }
    this(ref const(Helix) other)
    {
        a0 = other.a0; a1 = other.a1;
        r0 = other.r0; r1 = other.r1;
        theta01 = other.theta01;
        xdsh = other.xdsh; ydsh = other.ydsh; zdsh = other.zdsh;
    }
    override Helix dup() const
    {
        return new Helix(a0, a1, xdsh, r0, r1, theta01);
    }
    override Vector3 opCall(double t) const
    {
        double r = r0*(1.0-t) + r1*t;
        double theta = theta01 * t;
        Vector3 p = r*cos(theta)*xdsh + r*sin(theta)*ydsh + a0*(1.0-t) + a1*t;
        return p;
    }
    override string toString() const
    {
        return "Helix(a0=" ~ to!string(a0) ~
            ", a1=" ~ to!string(a1) ~
            ", xdsh=" ~ to!string(xdsh) ~
            ", r0=" ~ to!string(r0) ~
            ", r1=" ~ to!string(r1) ~
            ", theta01=" ~ to!string(theta01) ~
            ")";
    }
    override string classString() const
    {
        return "Helix";
    }
} // end class Helix

version(helix_test) {
    import util.msg_service;
    int main() {
        auto axis0 = Vector3([0.0, 0.0, 0.0]);
        auto axis1 = Vector3([1.0, 0.0, 0.0]);
        auto pstart = Vector3([0.0, 1.0, 0.0]);
        auto pend = Vector3([1.0, 0.0, 1.0]);
        auto h1 = new Helix(pstart, pend, axis0, axis1);
        auto p = h1(0.5);
        assert(approxEqualVectors(p, Vector3(0.5, 0.7071068, 0.7071068)), failedUnitTest());
        auto a0 = Vector3([0.0, 0.0, 0.0]);
        auto a1 = Vector3([1.0, 0.0, 0.0]); // axis is in global-frame x-direction
        auto xlocal = Vector3([0.0, 1.0, 0.0]); // pointing at start point
        auto r0 = 1.0;
        auto r1 = 1.0;
        auto dtheta = PI/2;
        auto h2 = new Helix(a0, a1, xlocal, r0, r1, dtheta);
        auto p2 = h2(0.5);
        assert(approxEqualVectors(p2, Vector3(0.5, 0.7071068, 0.7071068)),
               failedUnitTest()); // "Helix from fundamental parameters"
        return 0;
    }
} // end helix_test
