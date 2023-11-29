// modifiedpath.d
// Peter J. 2017-11-29: Split out of the original path module.

module geom.gpath.modifiedpath;

import std.conv;
import std.math;

import ntypes.complex;
import nm.number;

import geom.elements;
import geom.gpath.path;


class ReParameterizedPath : Path {
public:
    Path underlying_path;
    this(const Path other)
    {
        underlying_path = other.dup();
    }
    override Vector3 opCall(double t) const
    {
        double tdsh = underlying_t(t);
        return underlying_path(tdsh);
    }
    override Vector3 dpdt(double t) const
    {
        double tdsh = underlying_t(t);
        return underlying_path.dpdt(tdsh) * d_underlying_t_dt(tdsh);
    }
    override Vector3 d2pdt2(double t) const
    {
        double tdsh = underlying_t(t);
        return underlying_path.d2pdt2(tdsh) * pow(d_underlying_t_dt(tdsh),2) +
            underlying_path.dpdt(tdsh) * d2_underlying_t_dt2(tdsh);
    }
    abstract override string toString() const
    {
        return "ReParameterizedPath(underlying_path=" ~ to!string(underlying_path) ~ ")";
    }
    override string classString() const
    {
        return "ReParameterizedPath";
    }

protected:
    abstract double underlying_t(double t) const;
    abstract double d_underlying_t_dt(double t) const;
    abstract double d2_underlying_t_dt2(double t) const;

} // end class ReParameterizedPath

class ArcLengthParameterizedPath : ReParameterizedPath {
public:
    this(const Path other)
    {
        super(other);
        set_arc_length_vector(100);
    }
    this(ref const(ArcLengthParameterizedPath) other)
    {
        this(other.underlying_path);
    }
    override ArcLengthParameterizedPath dup() const
    {
        return new ArcLengthParameterizedPath(this.underlying_path);
    }
    override string toString() const
    {
        return "ArcLengthParameterizedPath(underlying_path=" ~ to!string(underlying_path) ~ ")";
    }
    override string classString() const
    {
        return "ArcLengthParametrizedPath";
    }

protected:
    double[] arc_length_vector;
    void set_arc_length_vector(int N)
    {
        // Compute the arc_lengths for a number of sample points
        // so that these can later be used to do a reverse interpolation
        // on the evaluation parameter.
        arc_length_vector.length = 0;
        if ( N == 0 ) return;
        double dt = 1.0 / N;
        double L = 0.0;
        arc_length_vector ~= L;
        Vector3 p0 = underlying_path(0.0);
        Vector3 p1;
        foreach (i; 1 .. N+1) {
            p1 = underlying_path(dt * i);
            Vector3 dp = p1 - p0;
            L += geom.abs(dp).re;
            arc_length_vector ~= L;
            p0 = p1;
        }
    } // end set_arc_length_vector()
    override double underlying_t(double t) const
    {
        // The incoming parameter value, t, is proportional to arc_length fraction.
        if (t <= 0.0) { return 0.0; }
        if (t >= 1.0) { return 1.0; }
        // Do a reverse look-up from the arc_length fraction to the original t parameter
        // of the underlying Path.
        double L_target = t * arc_length_vector[$-1];
        // Starting from the right-hand end,
        // let's try to find a point to the left of L_target.
        // If the value is out of range, this should just result in
        // us extrapolating one of the end segments -- that's OK.
        int i = to!int(arc_length_vector.length) - 1;
        double dt = 1.0 / (arc_length_vector.length - 1);
        while ( L_target < arc_length_vector[i] && i > 0 ) i--;
        double frac = (L_target - arc_length_vector[i]) /
            (arc_length_vector[i+1] - arc_length_vector[i]);
        return (1.0 - frac) * dt*i + frac * dt*(i+1);
    }
    override double d_underlying_t_dt(double t) const
    {
        // input "t" is underlying_t
        // derivative of inverse fn
        return 1.0/d_arc_f_dt(t);
    }
    override double d2_underlying_t_dt2(double t) const
    {
        // input "t" is underlying_t
        // derivative of inverse fn
        return -d2_arc_f_dt2(t)/pow(d_arc_f_dt(t),3);
    }
private:
    double d_arc_f_dt(double t) const
    {
        // "t" is underlying_t
        Vector3 dpdt = underlying_path.dpdt(t);
        return geom.abs(dpdt).re / arc_length_vector[$-1];
    }
    double d2_arc_f_dt2(double t) const
    {
        // "t" is underlying_t
        //chain rule on d_arc_f_dt
        Vector3 dpdt = underlying_path.dpdt(t);
        dpdt.normalize(); // direction only
        Vector3 d2pdt2 = underlying_path.d2pdt2(t);
        return dot(dpdt,d2pdt2).re / arc_length_vector[$-1];
    }
} // end class ArcLengthParameterizedPath


class SubRangedPath : ReParameterizedPath {
public:
    double t0;
    double t1;
    this(const Path other, double newt0, double newt1)
    {
        super(other);
        t0 = newt0;
        t1 = newt1;
    }
    this(ref const(SubRangedPath) other)
    {
        this(other.underlying_path, other.t0, other.t1);
    }
    override SubRangedPath dup() const
    {
        return new SubRangedPath(this.underlying_path, this.t0, this.t1);
    }
    override string toString() const
    {
        return "SubRangedPath(underlying_path=" ~ to!string(underlying_path)
            ~ ", t0=" ~ to!string(t0) ~ " t1=" ~ to!string(t1) ~ ")";
    }
    override string classString() const
    {
        return "SubRangedPath";
    }
protected:
    override double underlying_t(double t) const
    {
        return t0 + (t1 - t0)*t;
    }
    override double d_underlying_t_dt(double t) const
    {
        return t1 - t0;
    }
    override double d2_underlying_t_dt2(double t) const
    {
        return 0.0;
    }
} // end class SubRangedPath


class ReversedPath : SubRangedPath {
    // Just a particular case of SubRangedPath
    this(const Path other)
    {
        super(other, 1.0, 0.0);
    }
    override ReversedPath dup() const
    {
        return new ReversedPath(this.underlying_path);
    }
} // end class ReversedPath


class TransformedPath : Path {
public:
    Path original_path;
    this(const Path other)
    {
        original_path = other.dup();
    }
    override Vector3 opCall(double t) const
    {
        Vector3 p = original_path(t);
        return apply_transform(p);
    }
    abstract override string toString() const
    {
        return "TransformedPath(original_path=" ~ to!string(original_path) ~ ")";
    }
    override string classString() const
    {
        return "TransformedPath";
    }

protected:
    abstract Vector3 apply_transform(ref Vector3 p) const;

} // end class TransformedPath


class TranslatedPath : TransformedPath {
    Vector3 shift;
    this(const Path other, const Vector3 shift)
    {
        super(other);
        this.shift = shift;
    }
    override TranslatedPath dup() const
    {
        return new TranslatedPath(this.original_path, this.shift);
    }
    override string toString() const
    {
        return "TranslatedPath(original_path=" ~ to!string(original_path)
            ~ ", shift=" ~ to!string(shift) ~ ")";
    }
    override string classString() const
    {
        return "TranslatedPath";
    }

protected:
    override Vector3 apply_transform(ref Vector3 p) const
    {
        return p+shift;
    }
} // end class TranslatedPath


class MirrorImagePath : TransformedPath {
    Vector3 point;
    Vector3 normal;
    this(const Path other, const Vector3 point, const Vector3 normal)
    {
        super(other);
        this.point = point;
        this.normal = normal;
    }
    override MirrorImagePath dup() const
    {
        return new MirrorImagePath(this.original_path, this.point, this.normal);
    }
    override string toString() const
    {
        return "MirrorImagePath(original_path=" ~ to!string(original_path)
            ~ ", point=" ~ to!string(point) ~ ", normal=" ~ to!string(normal) ~ ")";
    }
    override string classString() const
    {
        return "MirrorImagePath";
    }

protected:
    override Vector3 apply_transform(ref Vector3 p) const
    {
        return p.mirror_image(point, normal);
    }
} // end class MirrorImagePath


class RotatedAboutZAxisPath : TransformedPath {
    double dtheta; // in radians
    this(const Path other, double angle)
    {
        super(other);
        dtheta = angle;
    }
    override RotatedAboutZAxisPath dup() const
    {
        return new RotatedAboutZAxisPath(this.original_path, this.dtheta);
    }
    override string toString() const
    {
        return "RotatedAboutZAxisPath(original_path=" ~ to!string(original_path)
            ~ ", angle=" ~ to!string(dtheta) ~ ")";
    }
    override string classString() const
    {
        return "RotatedAboutZAxisPath";
    }

protected:
    override Vector3 apply_transform(ref Vector3 p) const
    {
        return p.rotate_about_zaxis(dtheta);
    }
} // end class RotatedAboutZAxisPath


version(modifiedpath_test) {
    import util.msg_service;
    int main() {
        import geom.gpath.line;
        import geom.gpath.arc;
        import geom.gpath.polyline;
        import geom.gpath.bezier;
        auto a = Vector3([0.0, 0.0, 0.0]);
        auto b = Vector3([1.0, 1.0, 1.0]);
        auto c = Vector3([4.0, 4.0, 4.0]);
        auto abc = new Bezier([a, b, c]);
        auto abc_dsh = new ArcLengthParameterizedPath(abc);
        auto f = abc_dsh(0.5);
        assert(approxEqualVectors(f, Vector3(2,2,2)), failedUnitTest());

        a = Vector3([2.0, 2.0, 0.0]);
        b = Vector3([1.0, 2.0, 1.0]);
        c = Vector3([1.0, 2.0, 0.0]);
        auto acb = new ArcLengthParameterizedPath(new Bezier([a, c, b]));
        auto L = acb.underlying_path.length();
        auto dA = Vector3(-1, 0, 1);
        auto dAdt = abs(dA)/L;
        Vector3 d2A = Vector3(-1, 0, 1);
        auto d2Adt2 = dot(unit(d2A),Vector3(2,0,2))/L;
        // check to finite-difference in Path
        assert(approxEqualVectors(acb.dpdt(0.5), acb.Path.dpdt(0.5)),
               failedUnitTest());
        assert(approxEqualVectors(acb.d2pdt2(0.5), acb.Path.d2pdt2(0.5)),
               failedUnitTest());
        // the following checks are kind of redundant since they just follow
        // the same math as the function definitions
        assert(approxEqualVectors(acb.dpdt(0.5), Vector3(-1, 0, 1)/dAdt),
               failedUnitTest());
        assert(approxEqualVectors(acb.d2pdt2(0.5),
                                  Vector3(2,0,2)/pow(dAdt,2)-
                                  Vector3(-1, 0, 1)*d2Adt2/pow(dAdt,3)),
               failedUnitTest());

        a = Vector3([2.0, 2.0, 0.0]);
        b = Vector3([1.0, 2.0, 1.0]);
        c = Vector3([1.0, 2.0, 0.0]);
        abc = new Arc(a, b, c);
        auto polyline = new Polyline([abc, new Line(b, c)]);
        auto rev_poly = new ReversedPath(polyline);
        assert(approxEqualVectors(polyline(0.25), rev_poly(0.75)), failedUnitTest());
        acb = new SubRangedPath(new Bezier([a, c, b]), 0.5, 0.75);
        assert(approxEqualVectors(acb.dpdt(0), Vector3(-1, 0, 1)/4), failedUnitTest());
        assert(approxEqualVectors(acb.d2pdt2(0), Vector3(2,0,2)/16), failedUnitTest());
        auto r_acb = new ReversedPath(new Bezier([a, c, b]));
        assert(approxEqualVectors(r_acb.dpdt(0.5), -Vector3(-1, 0, 1)), failedUnitTest());
        assert(approxEqualVectors(r_acb.d2pdt2(0.5), Vector3(2,0,2)), failedUnitTest());

        a = Vector3([2.0, 0.0, 0.0]);
        b = Vector3([0.0, 2.0, 0.0]);
        c = Vector3([0.0, 0.0, 0.0]);
        abc = new Arc(a, b, c);
        auto abc_rotated = new RotatedAboutZAxisPath(abc, PI/4);
        // writeln("abc_rotated(1.0)=", abc_rotated(1.0));
        assert(approxEqualVectors(abc_rotated(1.0), sqrt(2.0)*Vector3(-1, 1, 0)),
               failedUnitTest());
        return 0;
    }
} // end mofifiedpath_test
