// subrangedvolume.d

module geom.volume.subrangedvolume;

import std.conv;
import geom.elements;
import geom.gpath;
import geom.surface;
import geom.volume.parametricvolume;


class SubRangedVolume : ParametricVolume {
public:
    ParametricVolume underlying_volume;
    double r0; // to subrange r, when evaluating a point in the volume
    double r1;
    double s0;
    double s1;
    double t0;
    double t1;

    this(const ParametricVolume pvolume,
         double r0=0.0, double r1=1.0,
         double s0=0.0, double s1=1.0,
         double t0=0.0, double t1=1.0)
    {
        underlying_volume = pvolume.dup();
        this.r0 = r0;
        this.r1 = r1;
        this.s0 = s0;
        this.s1 = s1;
        this.t0 = t0;
        this.t1 = t1;
    }
    override Vector3 opCall(double r, double s, double t) const
    {
        r = r0 + (r1-r0)*r; // subrange the parameter
        s = s0 + (s1-s0)*s;
        t = t0 + (t1-t0)*t;
        return underlying_volume(r, s, t);
    }
    override ParametricVolume dup() const
    {
        return new SubRangedVolume(underlying_volume,
                                   r0, r1, s0, s1, t0, t1);
    }
    override string toString() const
    {
        return "SubRangedVolume(underlying_volume=" ~
            to!string(underlying_volume) ~
            ", r0=" ~ to!string(r0) ~ ", r1=" ~ to!string(r1) ~
            ", s0=" ~ to!string(s0) ~ ", s1=" ~ to!string(s1) ~
            ", t0=" ~ to!string(t0) ~ ", t1=" ~ to!string(t1) ~
            ")";
    }
} // end class SubRangedVolume

// [TODO] test
