/**
 * nurbs.d
 * Author: Reece O.
 * Date: 2021-02-24
 */

module geom.gpath.nurbs;

import std.conv;
import std.stdio;
import std.format;

import geom.elements;
import geom.gpath.path;
import geom.misc.nurbs_utils;

class NURBS : Path {
public:

    this(const double[4][] Pw, const double[] U, int p)
    {
	    // Test Pw is minimum viable
	    if (Pw.length < 2) {
	        string errMsg = "NURBS() A NURBS path requires at least two control points.\n";
	        errMsg ~= format("Supplied number of control points: %s", Pw.length);
	        throw new Error(text(errMsg));
	    }
	    // Test p is valid
	    if (p > Pw.length-1) {
	        string errMsg = "NURBS() Degree of NURBS path is not compatible with given number of control points.\n";
	        errMsg ~= format("Supplied path degree: %s\n", p);
	        errMsg ~= format("Maximum allowed path degree: %s", Pw.length-1);
	        throw new Error(text(errMsg));
	    }
        // Test knot vector is valid
        this.m_a = to!int(Pw.length-1);
	    if (U.length != m_a+p+2) {
            string errMsg = "NURBS() Knot vector is not compatible with given number control points and curve degree.\n";
            errMsg ~= format("Supplied number of knots: %s\n", U.length);
            errMsg ~= format("Required number of knots: %s", m_a+p+2);
            throw new Error(text(errMsg));
        }

        this.mPw = Pw.dup;
        this.mU = U.dup;
        this.m_p = p;
        mN.length = p + 1;
        mNws = NURBSWorkspace(p);

    }
    this(ref const NURBS other)
    {
        this.mPw = other.mPw.dup;
        this.mU = other.mU.dup;
        this.m_p = other.m_p;
        this.m_a = other.m_a;
        mN.length = other.m_p + 1;
        mNws = NURBSWorkspace(other.m_p);
    }
    override NURBS dup() const
    {
        return new NURBS(mPw, mU, m_p);
    }
    override Vector3 opCall(double u) const
    {
        return deBoor(u);
    }
    override string toString() const
    {
        return "NURBS(Pw=" ~ to!string(mPw) ~ ", U=" ~ to!string(mU) ~ ", p=" ~ to!string(m_p) ~ ")";
    }
    override string classString() const
    {
        return "NURBS";
    }

private:
    double[4][] mPw; // collection of weighted control points
    double[] mU; // knot vector
    int m_p; // degree of curve
    int m_a; // a+1 = number of control points
    static double[] mN;
    static double[4] mCw;
    static double[3] mC;
    static NURBSWorkspace mNws;

    Vector3 deBoor(double u) const {
        // Returns the Cartesian coordinates of a point on a NURBS curve at a given parameter value
        // This is algorithm A4.1 from Piegl and Tiller (1997) - 'The NURBS Book'
        int span = findSpan(u, m_a, m_p, mU);
        basisFuns(span, u, m_p, mU, mN, mNws);
        mCw[] = 0.0;
        mC[] = 0.0;
        foreach (i; 0 .. m_p+1) mCw[] += mN[i]*mPw[span-m_p+i][];
        foreach (i; 0 .. 3)  mC[i] = mCw[i]/mCw[3];
        return Vector3(mC);
    }

}

version(nurbs_test) {
    import util.msg_service;
    int main() {
        double[4][] Pw = [[-4.0, -4.0, 0.0, 1.0], [-2.0, 4.0, 0.0, 1.0], [2.0*5, -4.0*5, 0.0, 5], [4.0, 4.0, 0.0, 1.0], [3.778, 1.836, 2.933, 1.0], [4*2.772, -4*3.875, 4*1.736, 4.0]];
        int p = 2;
        double[] U = [0.0, 0.0, 0.0, 0.375, 0.5, 0.625, 1.0, 1.0, 1.0];
        auto ncurve = new NURBS(Pw, U, p);
        assert(approxEqualVectors(Vector3(3.782, 2.939, 0.435), ncurve(0.6)), failedUnitTest());
        return 0;
    }
}

