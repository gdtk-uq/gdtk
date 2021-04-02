// slabvolume.d

module geom.volume.slabvolume;

import std.conv;
import geom.elements;
import geom.gpath;
import geom.surface;
import geom.volume.parametricvolume;
import geom.volume.sweptsurfacevolume;


class SlabVolume : SweptSurfaceVolume {
public:
    Vector3 dz; // slab thickness

    this(const ParametricSurface face0123, const Vector3 dz)
    // The simplest use case is that face0123 surface is in the x,y-plane and
    // the thickness vector is in the z-direction, however, the volume can be
    // made with more general orientation and the thickness vector just needs
    // to be somewhat orthogonal, so that the volume does not collapse.
    // The face0123 surface will form the bottom surface of the volume.
    {
        this.dz = dz;
        Vector3 p0 = face0123(0.0,0.0);
        super(face0123, new Line(p0, p0+dz));
    }

    this(ref const(SlabVolume) other)
    {
        dz = other.dz;
        Vector3 p0 = face0123(0.0,0.0);
        super(other.face0123, new Line(p0, p0+dz));
    }

    override SlabVolume dup() const
    {
        return new SlabVolume(this.face0123, this.dz);
    }

    override string toString() const
    {
        string repr = "SlabVolume(face0123=" ~ to!string(face0123);
        repr ~= ", dz=" ~ to!string(dz) ~ ")";
        return repr;
    } // end toString
} // end SlabVolume

version(slabvolume_test) {
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
        auto my_box = new SlabVolume(my_face, p[4]-p[0]);
        auto d = my_box(0.1, 0.1, 0.5);
        assert(approxEqualVectors(d, Vector3(0.1, 0.2, 1.5)), failedUnitTest());
        return 0;
    }
} // end slabvolume_test
