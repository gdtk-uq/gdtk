// subrangedsurface.d

module geom.surface.subrangedsurface;

import std.conv;
import geom.elements;
import geom.surface.parametricsurface;


class SubRangedSurface : ParametricSurface {
public:
    ParametricSurface underlying_surface;
    double r0; // to subrange r, when evaluating a point on the surface
    double r1;
    double s0;
    double s1;

    this(const ParametricSurface psurf,
         double r0=0.0, double r1=1.0, double s0=0.0, double s1=1.0)
    {
        underlying_surface = psurf.dup();
        this.r0 = r0;
        this.r1 = r1;
        this.s0 = s0;
        this.s1 = s1;
    }
    override Vector3 opCall(double r, double s) const
    {
        r = r0 + (r1-r0)*r; // subrange the parameter
        s = s0 + (s1-s0)*s;
        return underlying_surface(r,s);
    }
    override ParametricSurface dup() const
    {
        return new SubRangedSurface(underlying_surface, r0, r1, s0, s1);
    }
    override string toString() const
    {
        return "SubRangedSurface(underlying_surface=" ~
            to!string(underlying_surface) ~
            ", r0=" ~ to!string(r0) ~ ", r1=" ~ to!string(r1) ~
            ", s0=" ~ to!string(s0) ~ ", s1=" ~ to!string(s1) ~
            ")";
    }
} // end class SubRangedSurface

// [TODO] version(subrangedsurface_test)
