// channelpatch.d

module geom.surface.channelpatch;

import std.conv;
import ntypes.complex;
import nm.number;

import geom.elements;
import geom.gpath;
import geom.surface.parametricsurface;

class ChannelPatch : ParametricSurface {
public:
    Path cSouth;
    Path cNorth;
    bool ruled;
    bool pure2D;
    Vector3 p00, p10, p11, p01;

    this(const Path cSouth, const Path cNorth, bool ruled=false, bool pure2D=false)
    {
        this.cSouth = cSouth.dup();
        this.cNorth = cNorth.dup();
        this.ruled = ruled;
        this.pure2D = pure2D;
        p00 = cSouth(0.0);
        p10 = cSouth(1.0);
        p01 = cNorth(0.0);
        p11 = cNorth(1.0);
    }

    this(ref const(ChannelPatch) other)
    {
        cSouth = other.cSouth.dup();
        cNorth = other.cNorth.dup();
        ruled = other.ruled;
        pure2D = other.pure2D;
        p00 = other.p00;
        p10 = other.p10;
        p01 = other.p01;
        p11 = other.p11;
    }

    override ChannelPatch dup() const
    {
        return new ChannelPatch(this.cSouth, this.cNorth, this.ruled, this.pure2D);
    }

    override Vector3 opCall(double r, double s) const
    {
        auto bridge_path = make_bridging_path(r);
        Vector3 p = bridge_path(s);
        if (pure2D) { p.z = 0.0; }
        return p;
    }

    override string toString() const
    {
        return "ChannelPatch(cSouth=" ~ to!string(cSouth) ~ ", cNorth=" ~ to!string(cNorth) ~
            ", ruled=" ~ to!string(ruled) ~ ", pure2D=" ~ to!string(pure2D) ~ ")";
    }

    Path make_bridging_path(double r) const
    // For any value 0.0 <= r <= 1.0, this function returns the Path
    // that bridges the patch from a point A on the south curve to a
    // corresponding point B on the north curve.
    //
    // When using this surface patch in a simulation with adjacent
    // patches, it may be useful to call this function with values
    // of r=0.0 and r=1.0 to get the west and east edges, respectively,
    // of this ChannelPatch.
    // These edges can be used to construct conforming patches.
    //
    // Note that the Bezier, with all 4 points colinear and equally spaced,
    // will be equivalent to a straight line between its end points.
    // You can make your life easy by carefully setting the end points
    // and slopes of your defining south and north edge paths.
    {
        Vector3 pA = cSouth(r);
        Vector3 pB = cNorth(r);
        if (pure2D) { pA.z = 0.0; pB.z = 0.0; }
        if (ruled) {
            // Bridge with a straight line for a ruled surface.
            return new Line(pA, pB);
        } else {
            // Bridge with a Bezier3 path that is normal to both defining curves.
            Vector3 pBminuspA = pB - pA;
            double L = abs(pBminuspA).re;
            Vector3 tangentA = cSouth.dpdt(r);
            Vector3 tangentB = cNorth.dpdt(r);
            // Out-of-plane vectors
            Vector3 oopvA = cross(tangentA, pBminuspA);
            Vector3 oopvB = cross(tangentB, pBminuspA);
            // Inward-facing normal vectors on the surface.
            Vector3 nA = cross(oopvA, tangentA); nA.normalize();
            Vector3 nB = cross(tangentB, oopvB); nB.normalize();
            // Intermediate control points for the cubic Bezier.
            Vector3 p1 = pA + L/3.0*nA;
            Vector3 p2 = pB + L/3.0*nB;
            if (pure2D) { p1.z = 0.0; p2.z = 0.0; }
            return new Bezier([pA, p1, p2, pB]);
        }
    } // end make_bridging_path()
} // end class ChannelPatch

// [TODO] version(channelpatch_test)
