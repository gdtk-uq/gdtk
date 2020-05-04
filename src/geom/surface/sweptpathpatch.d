// sweptpathpatch.d

module geom.surface.sweptpathpatch;

import std.conv;

import geom.elements;
import geom.gpath;
import geom.surface.parametricsurface;

class SweptPathPatch : ParametricSurface {
public:
    Path cA; // Path to sweep (think of it as a west edge)
    Path cB; // Path along which cA is swept (think of it as a south edge)
    Vector3 p00, p10, p11, p01;

    this(const Path cA, const Path cB)
    {
        this.cA = cA.dup();
        this.cB = cB.dup();
        p00 = cB(0.0);
        p10 = cB(1.0);
        p01 = cB(0.0) + cA(1.0) - cA(0.0);
        p11 = cB(1.0) + cA(1.0) - cA(0.0);
    }

    this(ref const(SweptPathPatch) other)
    {
        cA = other.cA.dup();
        cB = other.cB.dup();
        p00 = other.p00;
        p10 = other.p10;
        p01 = other.p01;
        p11 = other.p11;
    }

    override SweptPathPatch dup() const
    {
        return new SweptPathPatch(this.cA, this.cB);
    }

    override Vector3 opCall(double r, double s) const
    {
        Vector3 p = cB(r) + cA(s) - cA(0.0);
        return p;
    }

    override string toString() const
    {
        return "SweptPathPatch(cA=" ~ to!string(cA) ~ ", cB=" ~ to!string(cB) ~ ")";
    }
} // end class SweptPathPatch

// [TODO] version(sweptpathpatch)
