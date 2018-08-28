// channelpatch.d

module geom.surface.channelpatch;

import std.conv;
import nm.complex;
import nm.number;

import geom.elements;
import geom.gpath;
import geom.surface.parametricsurface;

class ChannelPatch : ParametricSurface {
public:
    Path cA; // south edge
    Path cB; // north edge
    bool ruled;
    bool pure2D;
    Vector3 p00, p10, p11, p01;
    
    this(const Path cA, const Path cB, bool ruled=false, bool pure2D=false)
    {
        this.cA = cA.dup();
        this.cB = cB.dup();
        this.ruled = ruled;
        this.pure2D = pure2D;
        p00 = cA(0.0);
        p10 = cA(1.0);
        p01 = cB(0.0);
        p11 = cB(1.0);
    }
    
    this(ref const(ChannelPatch) other)
    {
        cA = other.cA.dup();
        cB = other.cB.dup();
        ruled = other.ruled;
        pure2D = other.pure2D;
        p00 = other.p00;
        p10 = other.p10;
        p01 = other.p01;
        p11 = other.p11;
    }

    override ChannelPatch dup() const
    {
        return new ChannelPatch(this.cA, this.cB, this.ruled, this.pure2D);
    }

    override Vector3 opCall(double r, double s) const 
    {
        auto bridge_path = make_bridging_path(r);
        Vector3 p = bridge_path(s);
        if (pure2D) { p.refz = 0.0; }
        return p;
    }

    override string toString() const
    {
        return "ChannelPatch(cA=" ~ to!string(cA) ~ ", cB=" ~ to!string(cB) ~
            ", ruled=" ~ to!string(ruled) ~ ", pure2D=" ~ to!string(pure2D) ~ ")";
    }
    
    Path make_bridging_path(double r) const
    {
        Vector3 cAr = cA(r); 
        Vector3 cBr = cB(r);
        if (pure2D) { cAr.refz = 0.0; cBr.refz = 0.0; }
        if (ruled) {
            // Bridge with a straight line for a ruled surface.
            return new Line(cAr, cBr);
        } else {
            // Bridge with a Bezier3 path that is normal to both defining curves.
            Vector3 pBminuspA = cBr - cAr;
            double L = abs(pBminuspA).re;
            Vector3 dcArdt = cA.dpdt(r);
            Vector3 dcBrdt = cB.dpdt(r);
            // Out-of-plane vectors
            Vector3 oopvA = cross(dcArdt, pBminuspA);
            Vector3 oopvB = cross(dcBrdt, pBminuspA);
            // Inward-facing normal vectors on the surface.
            Vector3 nA = cross(oopvA, dcArdt); nA.normalize();
            Vector3 nB = cross(dcBrdt, oopvB); nB.normalize();
            // Intermediate control points for the cubic Bezier.
            Vector3 p1 = cAr + L/3.0*nA;
            Vector3 p2 = cBr + L/3.0*nB;
            if (pure2D) { p1.refz = 0.0; p2.refz = 0.0; }
            return new Bezier([cAr, p1, p2, cBr]);
        }
    } // end make_bridging_path()
} // end class ChannelPatch

// [TODO] version(channelpatch_test)
