// twosurfacevolume.d

module geom.volume.twosurfacevolume;

import std.conv;
import geom.elements;
import geom.gpath;
import geom.surface;
import geom.volume.parametricvolume;


class TwoSurfaceVolume : ParametricVolume {
public:
    ParametricSurface face0123; // The bottom surface.
    ParametricSurface face4567; // The top surface.
    Vector3[8] p; // Corner points for the defined volume.

    this(const ParametricSurface face0123, const ParametricSurface face4567)
    {
        this.face0123 = face0123.dup();
        this.face4567 = face4567.dup();
        p[0] = face0123(0.0, 0.0);
        p[1] = face0123(1.0, 0.0);
        p[2] = face0123(1.0, 1.0);
        p[3] = face0123(0.0, 1.0);
        p[4] = face4567(0.0, 0.0);
        p[5] = face4567(1.0, 0.0);
        p[6] = face4567(1.0, 1.0);
        p[7] = face4567(0.0, 1.0);
    }

    this(ref const(TwoSurfaceVolume) other)
    {
        face0123 = other.face0123.dup();
        face4567 = other.face4567.dup();
        foreach(i; 0 .. 8) this.p[i] = other.p[i].dup();
    }

    override TwoSurfaceVolume dup() const
    {
        return new TwoSurfaceVolume(this.face0123, this.face4567);
    }

    override Vector3 opCall(double r, double s, double t) const
    // Locate a point within the volume by linear interpolation between
    // respective points on bottom and top surface.
    // Input:
    //     r: interpolation parameter i-direction west-->east, 0.0<=r<=1.0
    //     s: interpolation parameter j-direction south-->north, 0.0<=s<=1.0
    //     t: interpolation parameter k-direction bottom-->top, 0.0<=t<=1.0
    // Returns:
    //     a Vector3 value for the point.
    {
        Vector3 p_rst = face0123(r, s).scale(1.0-t) + face4567(r, s).scale(t);
        return p_rst;
    } // end opCall

    override string toString() const
    {
        string repr = "TwoSurfaceVolume(face0123=" ~ to!string(face0123);
        repr ~= ", face4567=" ~ to!string(face4567) ~ ")";
        return repr;
    } // end toString
} // end TwoSurfaceVolume

version(twosurfacevolume_test) {
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

        ParametricSurface my_face_bottom = new CoonsPatch(p[0], p[1], p[2], p[3]);
        ParametricSurface my_face_top = new CoonsPatch(p[4], p[5], p[6], p[7]);
        auto my_box = new TwoSurfaceVolume(my_face_bottom, my_face_top);
        auto d = my_box(0.1, 0.1, 0.5);
        assert(approxEqualVectors(d, Vector3(0.1, 0.2, 1.5)), failedUnitTest());
        return 0;
    }
} // twosurfacevolume_test
