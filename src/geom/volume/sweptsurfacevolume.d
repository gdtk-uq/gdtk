// sweptsurfacevolume.d

module geom.volume.sweptsurfacevolume;

import std.conv;
import geom.elements;
import geom.gpath;
import geom.surface;
import geom.volume.parametricvolume;


class SweptSurfaceVolume : ParametricVolume {
public:
    ParametricSurface face0123; // The bottom surface.
    Path edge04; // The path along which points from the bottom surface will be swept.
    // Note that the line edge04(0.0) location anchors the p0 corner of the volume.
    // Effectively, the bottom face will be moved to that location.
    Vector3[8] p; // Corner points for the defined volume.

    this(const ParametricSurface face0123, const Path edge04)
    {
        this.face0123 = face0123.dup();
        this.edge04 = edge04.dup();
        p[0] = edge04(0.0);
        p[1] = edge04(0.0) + face0123(1.0, 0.0) - face0123(0.0, 0.0);
        p[2] = edge04(0.0) + face0123(1.0, 1.0) - face0123(0.0, 0.0);
        p[3] = edge04(0.0) + face0123(0.0, 1.0) - face0123(0.0, 0.0);
        p[4] = edge04(1.0);
        p[5] = edge04(1.0) + face0123(1.0, 0.0) - face0123(0.0, 0.0);
        p[6] = edge04(1.0) + face0123(1.0, 1.0) - face0123(0.0, 0.0);
        p[7] = edge04(1.0) + face0123(0.0, 1.0) - face0123(0.0, 0.0);
    }

    this(ref const(SweptSurfaceVolume) other)
    {
        face0123 = other.face0123.dup();
        edge04 = other.edge04.dup();
        foreach(i; 0 .. 8) this.p[i] = other.p[i].dup();
    }

    override SweptSurfaceVolume dup() const
    {
        return new SweptSurfaceVolume(this.face0123, this.edge04);
    }

    override Vector3 opCall(double r, double s, double t) const
    // Locate a point within the volume by blended linear interpolation.
    // Input:
    //     r: interpolation parameter i-direction west-->east, 0.0<=r<=1.0
    //     s: interpolation parameter j-direction south-->north, 0.0<=s<=1.0
    //     t: interpolation parameter k-direction bottom-->top, 0.0<=t<=1.0
    // Returns:
    //     a Vector3 value for the point.
    {
        Vector3 p_rst = edge04(t) + face0123(r, s) - face0123(0.0, 0.0);
        return p_rst;
    } // end opCall

    override string toString() const
    {
        string repr = "SweptSurfaceVolume(face0123=" ~ to!string(face0123);
        repr ~= ", edge04=" ~ to!string(edge04) ~ ")";
        return repr;
    } // end toString
} // end SweptSurfaceVolume

version(sweptsurfacevolume_test) {
    import util.msg_service;
    int main() {
        Vector3[8] p;
        p[0] = Vector3(0.0, 0.1, 0.0);
        p[1] = Vector3(1.0, 0.1, 0.0);
        p[2] = Vector3(1.0, 1.1, 0.0);
        p[3] = Vector3(0.0, 1.1, 0.0);

        p[4] = Vector3(0.0, 0.1, 3.0);
        p[5] = Vector3(1.0, 0.1, 3.0);
        p[6] = Vector3(1.0, 1.1, 3.0);
        p[7] = Vector3(0.0, 1.1, 3.0);

        ParametricSurface my_face = new CoonsPatch(p[0], p[1], p[2], p[3]);
        auto my_box = new SweptSurfaceVolume(my_face, new Line(p[0],p[4]));
        auto d = my_box(0.1, 0.1, 0.5);
        assert(approxEqualVectors(d, Vector3(0.1, 0.2, 1.5)), failedUnitTest());
        return 0;
    }
} // sweptsurfacevolume_test
