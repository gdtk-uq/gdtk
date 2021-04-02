// twosurfacevolume.d
//
// Make a "ruled" volume between two surfaces, analogous to a ruled surface.
// A point within the volume is linearly interpolated between
// corresponding points on the faces.
//
// Authors: Ingo Jahn and Peter J.
//
// 2019-05-31 Ingo's original with face at t=0 and t=1.
// 2021-04-02 PJ generalize to have the pair of faces in any of the 3 directions.
//
// Note that the r,s coordinates local to the faces are aligned as per the
// faces bound to the TFIVolume.  The volume's r,s,t parameters are aligned
// with its i,j,k grid index directions, respectively.  Thes index directions
// are shown in the following ASCII figures.
//
// 3-------2  7-------6
// |^      |  |^      |
// |j  B   |  |j  T   |
// |   i-> |  |   i-> |
// 0-------1  4-------5
//
// 4-------5  7-------6
// |^      |  |^      |
// |k  S   |  |k  N   |
// |   i-> |  |   i-> |
// 0-------1  3-------2
//
// 4-------7  5-------6
// |^      |  |^      |
// |k  W   |  |k  E   |
// |   j-> |  |   j-> |
// 0-------3  1-------2
//
module geom.volume.twosurfacevolume;

import std.conv;
import geom.elements;
import geom.gpath;
import geom.surface;
import geom.volume.parametricvolume;


class TwoSurfaceVolume : ParametricVolume {
public:
    ParametricSurface face0; // The bottom surface.
    ParametricSurface face1; // The top surface.
    string ruled_direction; // one of "i", "j", "k", "r", "s", "t"
    Vector3[8] p; // Corner points for the defined volume.

    this(const ParametricSurface face0, const ParametricSurface face1, string ruled_direction="t")
    {
        this.face0 = face0.dup();
        this.face1 = face1.dup();
        this.ruled_direction = ruled_direction;
        switch (ruled_direction) {
        case "r":
        case "i":
            // face0 is west face
            p[0] = face0(0.0,0.0); p[3] = face0(1.0,0.0); p[7] = face0(1.0,1.0); p[4] = face0(0.0,1.0);
            // face1 is east face
            p[1] = face1(0.0,0.0); p[2] = face1(1.0,0.0); p[6] = face1(1.0,1.0); p[5] = face1(0.0,1.0);
            break;
        case "s":
        case "j":
            // face0 is south face
            p[0] = face0(0.0,0.0); p[1] = face0(1.0,0.0); p[5] = face0(1.0,1.0); p[4] = face0(0.0,1.0);
            // face1 in north face
            p[3] = face1(0.0,0.0); p[2] = face1(1.0,0.0); p[6] = face1(1.0,1.0); p[7] = face1(0.0,1.0);
            break;
        case "t":
        case "k":
            // face0 is bottom face
            p[0] = face0(0.0,0.0); p[1] = face0(1.0,0.0); p[2] = face0(1.0,1.0); p[3] = face0(0.0,1.0);
            // face1 is top face
            p[4] = face1(0.0,0.0); p[5] = face1(1.0,0.0); p[6] = face1(1.0,1.0); p[7] = face1(0.0,1.0);
            break;
        default:
            throw new Exception("Invalid string for ruled_direction: " ~ ruled_direction);
        }
    }

    this(ref const(TwoSurfaceVolume) other)
    {
        face0 = other.face0.dup();
        face1 = other.face1.dup();
        ruled_direction = other.ruled_direction;
        foreach(i; 0 .. 8) { this.p[i] = other.p[i].dup(); }
    }

    override TwoSurfaceVolume dup() const
    {
        return new TwoSurfaceVolume(face0, face1, ruled_direction);
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
        Vector3 p_rst;
        switch (ruled_direction) {
        case "r":
        case "i":
            p_rst = face0(s, t).scale(1.0-1) + face1(s, t).scale(r);
            break;
        case "s":
        case "j":
            p_rst = face0(r, t).scale(1.0-s) + face1(r, t).scale(s);
            break;
        case "t":
        case "k":
            p_rst = face0(r, s).scale(1.0-t) + face1(r, s).scale(t);
            break;
        default:
            throw new Exception("Invalid string for ruled_direction: " ~ ruled_direction);
        }
        return p_rst;
    } // end opCall

    override string toString() const
    {
        string repr = "TwoSurfaceVolume(face0=" ~ to!string(face0);
        repr ~= ", face1=" ~ to!string(face1);
        repr ~= ", ruled_direction=" ~ ruled_direction;
        repr ~= ")";
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
        auto my_box = new TwoSurfaceVolume(my_face_bottom, my_face_top, "t");
        auto d = my_box(0.1, 0.1, 0.5);
        assert(approxEqualVectors(d, Vector3(0.1, 0.2, 1.5)), failedUnitTest());
        return 0;
    }
} // twosurfacevolume_test
