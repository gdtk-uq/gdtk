/**
 * nurbssurface.d
 * Author: Reece O.
 * Date: 2021-03-09
 * TODO Add error messages to nurbs_utils
 */

module geom.surface.nurbssurface;

import std.conv;
import std.stdio;
import std.format;

import geom.elements;
import geom.surface.parametricsurface;
import geom.misc.nurbs_utils;

class NURBSSurface : ParametricSurface {
public:

    this(const double[4][][] Pw, const double[] U, int p, const double[] V, int q)
    {
        // Test Pw is minimum viable
	    if (Pw.length < 2 || Pw[0].length < 2) {
	        string errMsg = "NURBSSurface() A NURBS surface requires at least two control points in each direction.\n";
	        errMsg ~= format("Supplied number of control points in direction u: %s\n", Pw.length);
	        errMsg ~= format("Supplied number of control points in direction v: %s", Pw[0].length);
	        throw new Error(text(errMsg));
	    }
	    // Test that p is valid
	    if (p > Pw.length-1) {
	        string errMsg = "NURBSSurface() Degree of NURBS surface in direction u is not compatible with given number of control points.\n";
	        errMsg ~= format("Supplied degree in u direction: %s\n", p);
	        errMsg ~= format("Maximum allowed degree in u direction: %s", Pw.length-1);
	        throw new Error(text(errMsg));
	    }
	    // Test that q is valid
	    if (q > Pw[0].length-1) {
	        string errMsg = "NURBSSurface() Degree of NURBS surface in direction v is not compatible with given number of control points.\n";
	        errMsg ~= format("Supplied degree in u direction: %s\n", q);
	        errMsg ~= format("Maximum allowed degree in u direction: %s", Pw[0].length-1);
	        throw new Error(text(errMsg));
	    }
	    // Test U knot vector is valid
        this.m_a = to!int(Pw.length-1);
	    if (U.length != m_a+p+2) {
            string errMsg = "NURBS() Knot vector in direction u is not compatible with given number control points and surface degree.\n";
            errMsg ~= format("Supplied number of knots in direction u: %s\n", U.length);
            errMsg ~= format("Required number of knots in direction u: %s", m_a+p+2);
            throw new Error(text(errMsg));
        }
        // Test V knot vector is valid
        this.m_b = to!int(Pw[0].length-1);
	    if (V.length != m_b+q+2) {
            string errMsg = "NURBS() Knot vector in direction u is not compatible with given number control points and surface degree.\n";
            errMsg ~= format("Supplied number of knots in direction v: %s\n", V.length);
            errMsg ~= format("Required number of knots in direction v: %s", m_b+q+2);
            throw new Error(text(errMsg));
        }

        this.mPw = Pw.dup;
        this.mU = U.dup;
        this.m_p = p;
        this.mV = V.dup;
        this.m_q = q;
        mNu.length = p + 1;
        mNws_u = NURBSWorkspace(p);
        mNv.length = q + 1;
        mNws_v = NURBSWorkspace(q);
	    temp.length = m_q+1;
    }
    this(ref const NURBSSurface other)
    {
        this.mPw = other.mPw.dup;
        this.mU = other.mU.dup;
        this.m_p = other.m_p;
        this.m_a = other.m_a;
        this.mV = other.mV.dup;
        this.m_q = other.m_q;
        this.m_b = other.m_p;
        mNu.length = other.m_p + 1;
        mNws_u = NURBSWorkspace(other.m_p);
        mNv.length = other.m_q + 1;
        mNws_v = NURBSWorkspace(other.m_q);
	    temp.length = m_q+1;
    }
    override NURBSSurface dup() const
    {
        return new NURBSSurface(mPw, mU, m_p, mV, m_q);
    }
    override Vector3 opCall(double u, double v) const
    {
        return deBoor(u, v);
    }
    override string toString() const
    {
        return "NURBSSurface(Pw=" ~ to!string(mPw) ~ ", U=" ~ to!string(mU) ~ ", p=" ~ to!string(m_p) ~ ", V=" ~ to!string(mV) ~ ", q=" ~ to!string(m_q) ~ ")";
    }

private:
    const double[4][][] mPw; // array of weighted control points

    const double[] mU; // knot vector in direction u
    int m_p; // degree of curve in direction u
    int m_a; // a+1 = number of control points in direction u

    const double[] mV; // knot vector in direction v
    int m_q; // degree of curve in direction v
    int m_b; // b+1 = number of control points in direction v

    static double[4] mSw;
    static double[3] mS;
    static double[] mNu;
    static double[] mNv;
    static double[4][] temp;
    static NURBSWorkspace mNws_u;
    static NURBSWorkspace mNws_v;

    Vector3 deBoor(double u, double v) const {
        // Returns the Cartesian coordinates of a point on a NURBS surface at a given set of parameter values
        // This is algorithm A4.3 from Piegl and Tiller (1997) - 'The NURBS Book'
        int uspan = findSpan(u, m_a, m_p, mU);
        basisFuns(uspan, u, m_p, mU, mNu, mNws_u);
        int vspan = findSpan(v, m_b, m_q, mV);
        basisFuns(vspan, v, m_q, mV, mNv, mNws_v);
        foreach (l; 0 .. m_q+1) {
            temp[l][] = 0.0;
            foreach (k; 0 .. m_p+1) temp[l][] += mNu[k]*mPw[uspan-m_p+k][vspan-m_q+l][];
        }
        mSw = 0.0;
        mS = 0.0;
        foreach (l; 0 .. m_q+1) mSw[] += mNv[l]*temp[l][];
        foreach (i; 0 .. 3) mS[i] = mSw[i]/mSw[3];
        return Vector3(mS);
    }
}

version(nurbssurface_test) {
    import util.msg_service;
    int main () {
    // This is example Ex4.3 from Piegl and Tiller (1997) - 'The NURBS Book'
    // Example only included 'activated' control points, so just set all others to 0 since they have no contribution
    double[4][][] Pw = [[[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]],
                        [[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]],
                        [[0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 1.0], [0.0, 6.0, 4.0, 2.0], [0.0, 2.0, 0.0, 1.0], [0.0, 0.0, 0.0, 0.0]],
                        [[0.0, 0.0, 0.0, 0.0], [4.0, 6.0, 8.0, 2.0], [12.0, 24.0, 12.0, 6.0], [4.0, 6.0, 0.0, 2.0], [0.0, 0.0, 0.0, 0.0]],
                        [[0.0, 0.0, 0.0, 0.0], [4.0, 2.0, 4.0, 1.0], [8.0, 6.0, 4.0, 2.0], [4.0, 2.0, 0.0, 1.0], [0.0, 0.0, 0.0, 0.0]],
                        [[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]],
                        [[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]],
                        [[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]]];
    int p = 2;
    int q = 2;
    double[] U = [0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 5.0, 5.0, 5.0];
    double[] V = [0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0];

    auto nsurf = new NURBSSurface(Pw, U, p, V, q);
    assert(approxEqualVectors(Vector3(2.0, 98.0/27.0, 68.0/27.0), nsurf.opCall(5.0/2.0, 1.0)), failedUnitTest());
    return 0;
    }
}
