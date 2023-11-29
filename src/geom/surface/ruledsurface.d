// ruledsurface.d

module geom.surface.ruledsurface;

import std.conv;
import ntypes.complex;
import nm.number;

import geom.elements;
import geom.gpath;
import geom.surface.parametricsurface;

class RuledSurface : ParametricSurface {
public:
    Path edge0;
    Path edge1;
    string ruled_direction;
    bool pure2D;
    Vector3 p00, p10, p11, p01;

    this(const Path edge0, const Path edge1, string ruled_direction="s", bool pure2D=false)
    {
        this.edge0 = edge0.dup();
        this.edge1 = edge1.dup();
        this.ruled_direction = ruled_direction;
        this.pure2D = pure2D;
        switch (ruled_direction) {
        case "r":
        case "i":
            // edge0 is west boundary
            p00 = edge0(0.0); p01 = edge0(1.0);
            // edge1 is east boundary
            p10 = edge1(0.0); p11 = edge1(1.0);
            break;
        case "s":
        case "j":
            // edge0 is south boundary
            p00 = edge0(0.0); p10 = edge0(1.0);
            // edge1 in north boundary
            p01 = edge1(0.0); p11 = edge1(1.0);
            break;
        default:
            throw new Exception("Invalid string for ruled_direction: " ~ ruled_direction);
        }
    }

    this(ref const(RuledSurface) other)
    {
        edge0 = other.edge0.dup();
        edge1 = other.edge1.dup();
        ruled_direction = other.ruled_direction;
        pure2D = other.pure2D;
        p00 = other.p00;
        p10 = other.p10;
        p01 = other.p01;
        p11 = other.p11;
    }

    override RuledSurface dup() const
    {
        return new RuledSurface(this.edge0, this.edge1, this.ruled_direction, this.pure2D);
    }

    override Vector3 opCall(double r, double s) const
    {
        Vector3 p;
        switch (ruled_direction) {
        case "r":
        case "i":
            p = (1.0-r)*edge0(s) + r*edge1(s);
            break;
        case "s":
        case "j":
            p = (1.0-s)*edge0(r) + s*edge1(r);
            break;
        default:
            throw new Exception("Invalid string for ruled_direction: " ~ ruled_direction);
        }
        if (pure2D) { p.z = 0.0; }
        return p;
    }

    override string toString() const
    {
        return "RuledSurface(edge0=" ~ to!string(edge0) ~ ", edge1=" ~ to!string(edge1) ~
            ", ruled_direction=" ~ ruled_direction ~ ", pure2D=" ~ to!string(pure2D) ~ ")";
    }
} // end class RuledSurface

// [TODO] version(ruledsurface_test)
