// coonspatch.d

module geom.surface.coonspatch;

import std.conv;
import ntypes.complex;
import nm.number;

import geom.elements;
import geom.gpath;
import geom.surface.parametricsurface;


class CoonsPatch : ParametricSurface {
public:
    Path north, east, south, west; // bounding paths
    Vector3 p00, p10, p11, p01;    // corners

    this(in Vector3 p00, in Vector3 p10, in Vector3 p11, in Vector3 p01)
    {
        north = new Line(p01, p11);
        east = new Line(p10, p11);
        south = new Line(p00, p10);
        west = new Line(p00, p01);
        this(south, north, west, east);
    }

    this(in Path south, in Path north, in Path west, in Path east)
    // The particular order for the boundary surfaces goes way back
    // to the original grid generation paper, so it doesn't match the
    // default order of NESW in the rest of the flow code.
    {
        this.north = north.dup();
        this.east = east.dup();
        this.south = south.dup();
        this.west = west.dup();
        p00 = south(0.0);
        p10 = south(1.0);
        p01 = north(0.0);
        p11 = north(1.0);
        // Check alternate evaluation of corners for consistency.
        Vector3 p00_alt = west(0.0);
        Vector3 p10_alt = east(0.0);
        Vector3 p01_alt = west(1.0);
        Vector3 p11_alt = east(1.0);
        if (!approxEqualVectors(p00, p00_alt)) {
            throw new Error(text("CoonsPatch open corner p00= ", p00,
                                 " p00_alt= ", p00_alt));
        }
        if (!approxEqualVectors(p10, p10_alt)) {
            throw new Error(text("CoonsPatch open corner p10= ", p10,
                                 " p10_alt= ", p10_alt));
        }
        if (!approxEqualVectors(p11, p11_alt)) {
            throw new Error(text("CoonsPatch open corner p11= ", p11,
                                 " p11_alt= ", p11_alt));
        }
        if (!approxEqualVectors(p01, p01_alt)) {
            throw new Error(text("CoonsPatch open corner p01= ", p01,
                                 " p01_alt= ", p01_alt));
        }
    }

    this(ref const(CoonsPatch) other)
    {
        this.north = other.north.dup();
        this.east = other.east.dup();
        this.south = other.south.dup();
        this.west = other.west.dup();
        p00 = other.p00;
        p10 = other.p10;
        p01 = other.p01;
        p11 = other.p11;
    }

    override CoonsPatch dup() const
    {
        return new CoonsPatch(this.south, this.north, this.west, this.east);
    }

    override Vector3 opCall(double r, double s) const
    {
        Vector3 south_r = south(r);
        Vector3 north_r = north(r);
        Vector3 west_s = west(s);
        Vector3 east_s = east(s);
        Vector3 p = (1.0-s)*south_r + s*north_r + (1.0-r)*west_s + r*east_s -
            ((1.0-r)*(1.0-s)*p00 + (1.0-r)*s*p01 + r*(1.0-s)*p10 + r*s*p11);
        return p;
    }

    override string toString() const
    {
        return "CoonsPatch(south=" ~ to!string(south) ~
            ", north=" ~ to!string(north) ~
            ", west=" ~ to!string(west) ~
            ", east=" ~ to!string(east) ~
            ")";
    }
} // end class CoonsPatch


version(coonspatch_test) {
    import util.msg_service;
    int main() {
        auto p00 = Vector3([0.0, 0.1, 3.0]);
        auto p10 = Vector3(1.0, 0.1, 3.0);
        auto p11 = Vector3(1.0, 1.1, 3.0);
        auto p01 = Vector3(0.0, 1.1, 3.0);
        auto my_patch = new CoonsPatch(p00, p10, p11, p01);
        auto c = my_patch(0.5, 0.5);
        assert(approxEqualVectors(c, Vector3(0.5, 0.6, 3.0)), failedUnitTest());
        c = my_patch(0.1, 0.1);
        assert(approxEqualVectors(c, Vector3(0.1, 0.2, 3.0)), failedUnitTest());
        return 0;
    }
}
