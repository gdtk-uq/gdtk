/** 
 * nurbsvolume.d
 * Author: Reece O. 
 * Date: 2021-03-09
 */

module geom.volume.nurbsvolume;

import std.conv;
import std.stdio;
import std.format;

import geom.elements;
import geom.volume.parametricvolume;
import nurbs_utils;

class NURBSVolume : ParametricVolume {
public:
    
    this(const double[4][][][] Pw, const double[] U, int p, const double[] V, int q, const double[] W, int r) 
    {
        // Test Pw is minimum viable
	    if (Pw.length < 2 || Pw[0].length < 2) {
	        string errMsg = "NURBSVolume() A NURBS volume requires at least two control points in each direction.\n";
	        errMsg ~= format("Supplied number of control points in direction u: %s\n", Pw.length);
	        errMsg ~= format("Supplied number of control points in direction v: %s", Pw[0].length);
	        throw new Error(text(errMsg));
	    }
	    // Test that p is valid
	    if (p > Pw.length-1) {
	        string errMsg = "NURBSVolume() Degree of NURBS volume in direction u is not compatible with given number of control points.\n";
	        errMsg ~= format("Supplied degree in u direction: %s\n", p);
	        errMsg ~= format("Maximum allowed degree in u direction: %s", Pw.length-1);
	        throw new Error(text(errMsg));
	    }
	    // Test that q is valid
	    if (q > Pw[0].length-1) {
	        string errMsg = "NURBSVolume() Degree of NURBS volume in direction v is not compatible with given number of control points.\n";
	        errMsg ~= format("Supplied degree in v direction: %s\n", q);
	        errMsg ~= format("Maximum allowed degree in v direction: %s", Pw[0].length-1);
	        throw new Error(text(errMsg));
	    }
	    // Test that r is valid
	    if (r > Pw[0][0].length-1) {
	        string errMsg = "NURBSVolume() Degree of NURBS volume in direction w is not compatible with given number of control points.\n";
	        errMsg ~= format("Supplied degree in w direction: %s\n", r);
	        errMsg ~= format("Maximum allowed degree in w direction: %s", Pw[0][0].length-1);
	        throw new Error(text(errMsg));
	    }
	    // Test U knot vector is valid
        this.m_a = to!int(Pw.length-1);
	    if (U.length != m_a+p+2) {
            string errMsg = "NURBSVolume() Knot vector in direction u is not compatible with given number control points and volume degree.\n";
            errMsg ~= format("Supplied number of knots in direction u: %s\n", U.length);
            errMsg ~= format("Required number of knots in direction u: %s", m_a+p+2);
            throw new Error(text(errMsg));
        }
        // Test V knot vector is valid
        this.m_b = to!int(Pw[0].length-1);
	    if (V.length != m_b+q+2) {
            string errMsg = "NURBSVolume() Knot vector in direction v is not compatible with given number control points and volume degree.\n";
            errMsg ~= format("Supplied number of knots in direction v: %s\n", V.length);
            errMsg ~= format("Required number of knots in direction v: %s", m_b+q+2);
            throw new Error(text(errMsg));
        }
        // Test W knot vector is valid
        this.m_c = to!int(Pw[0][0].length-1);
	    if (W.length != m_c+r+2) {
            string errMsg = "NURBSVolume() Knot vector in direction w is not compatible with given number control points and volume degree.\n";
            errMsg ~= format("Supplied number of knots in direction w: %s\n", W.length);
            errMsg ~= format("Required number of knots in direction w: %s", m_c+r+2);
            throw new Error(text(errMsg));
        }
        
        this.mPw = Pw.dup; 
        this.mU = U.dup; 
        this.m_p = p; 
        this.mV = V.dup; 
        this.m_q = q; 
        this.mW = W.dup; 
        this.m_r = r; 
        mNu.length = p + 1;
        mNws_u = NURBSWorkspace(p);
        mNv.length = q + 1;
        mNws_v = NURBSWorkspace(q);
        mNw.length = r + 1;
        mNws_w = NURBSWorkspace(r);
        temp1.length = m_q+1;
        foreach (i; 0 .. temp1.length) temp1[i].length = m_r+1;
	    temp2.length = m_r+1;
	    
    }
    this(ref const NURBSVolume other)
    {
        this.mPw = other.mPw.dup; 
        this.mU = other.mU.dup; 
        this.m_p = other.m_p; 
        this.m_a = other.m_a;
        this.mV = other.mV.dup; 
        this.m_q = other.m_q; 
        this.m_b = other.m_p;
        this.mW = other.mW.dup; 
        this.m_r = other.m_r; 
        this.m_c = other.m_c;
        mNu.length = other.m_p + 1;
        mNws_u = NURBSWorkspace(other.m_p);
        mNv.length = other.m_q + 1;
        mNws_v = NURBSWorkspace(other.m_q);
        mNw.length = other.m_r + 1;
        mNws_w = NURBSWorkspace(other.m_r);
        temp1.length = m_q+1;
        foreach (i; 0 .. temp1.length) temp1[i].length = m_r+1;
	    temp2.length = m_r+1;
    }
    override NURBSVolume dup() const
    {
        return new NURBSVolume(mPw, mU, m_p, mV, m_q, mW, m_r);
    }
    override Vector3 opCall(double u, double v, double w) const
    {
        return deBoor(u, v, w);
    }
    override string toString() const
    {
        return "NURBSVolume(Pw=" ~ to!string(mPw) ~ ", U=" ~ to!string(mU) ~ ", p=" ~ to!string(m_p) ~ ", V=" ~ to!string(mV) ~ ", q=" ~ to!string(m_q) ~ ", V=" ~ to!string(mW) ~ ", r=" ~ to!string(m_r) ~ ")";
    }

private:
    const double[4][][][] mPw; // array of weighted control points
    
    const double[] mU; // knot vector in direction u
    int m_p; // degree of volume in direction u
    int m_a; // a+1 = number of control points in direction u
    
    const double[] mV; // knot vector in direction v
    int m_q; // degree of volume in direction v
    int m_b; // b+1 = number of control points in direction v
    
    const double[] mW; // knot vector in direction w
    int m_r; // degree of volume in direction w
    int m_c; // c+1 = number of control points in direction w
    
    static double[4] mSw;
    static double[3] mS; 
    static double[] mNu;
    static double[] mNv;
    static double[] mNw;
    static double[4][][] temp1; //TODO change this
    static double[4][] temp2; //TODO change this
    static NURBSWorkspace mNws_u;
    static NURBSWorkspace mNws_v;
    static NURBSWorkspace mNws_w;
    
    Vector3 deBoor(double u, double v, double w) const {
        // Returns the Cartesian coordinates of a point on a NURBS volume at a given set of parameter values
        int uspan = FindSpan(u, m_a, m_p, mU);
        BasisFuns(uspan, u, m_p, mU, mNu, mNws_u);
        int vspan = FindSpan(v, m_b, m_q, mV);
        BasisFuns(vspan, v, m_q, mV, mNv, mNws_v);
        int wspan = FindSpan(w, m_c, m_r, mW);
        BasisFuns(wspan, w, m_r, mW, mNw, mNws_w);
        foreach (m; 0 .. m_r+1) {
            foreach (l; 0 .. m_q+1) {
                temp1[l][m][] = 0.0;
                foreach (k; 0 .. m_p+1) {
                    temp1[l][m][] += mNu[k] * mPw[uspan-m_p+k][vspan-m_q+l][wspan-m_r+m][];
                }
            }
        }
        foreach (m; 0 .. m_r+1) {
            temp2[m][] = 0.0;
            foreach (l; 0 .. m_q+1) temp2[m][] += mNv[l] * temp1[l][m][];
        }
        mSw = 0.0;
        mS = 0.0;
        foreach (m; 0 .. m_r+1) mSw[] += mNw[m]*temp2[m][];
        foreach (i; 0 .. 3) mS[i] = mSw[i]/mSw[3];
        return Vector3(mS);
    }
}

version(nurbsvolume_test) {
    import util.msg_service;
    import std.math;
    int main () {
        // solid cylinder point evaluation test
        double[4][][][] Pw = [[[[0.0, 0.0, -3.0, 1.0], [0.0, 0.0, -3.0, 1.0],           [0.0, 0.0, -3.0, 1.0], [0.0, 0.0, -3.0, 1.0],            [0.0, 0.0, -3.0, 1.0],  [0.0, 0.0, -3.0, 1.0],             [0.0, 0.0, -3.0, 1.0],  [0.0, 0.0, -3.0, 1.0],            [0.0, 0.0, -3.0, 1.0]],
                               [[0.0, 1.0, -3.0, 1.0], [1.0, 1.0, -3.0, 1.0/sqrt(2.0)], [1.0, 0.0, -3.0, 1.0], [1.0, -1.0, -3.0, 1.0/sqrt(2.0)], [0.0, -1.0, -3.0, 1.0], [-1.0, -1.0, -3.0, 1.0/sqrt(2.0)], [-1.0, 0.0, -3.0, 1.0], [-1.0, 1.0, -3.0, 1.0/sqrt(2.0)], [0.0, 1.0, -3.0, 1.0]]],
                              
                              [[[0.0, 0.0, -1.0, 1.0], [0.0, 0.0, -1.0, 1.0],           [0.0, 0.0, -1.0, 1.0], [0.0, 0.0, -1.0, 1.0],            [0.0, 0.0, -1.0, 1.0],  [0.0, 0.0, -1.0, 1.0],             [0.0, 0.0, -1.0, 1.0],  [0.0, 0.0, -1.0, 1.0],            [0.0, 0.0, -1.0, 1.0]],
                               [[0.0, 1.0, -1.0, 1.0], [1.0, 1.0, -1.0, 1.0/sqrt(2.0)], [1.0, 0.0, -1.0, 1.0], [1.0, -1.0, -1.0, 1.0/sqrt(2.0)], [0.0, -1.0, -1.0, 1.0], [-1.0, -1.0, -1.0, 1.0/sqrt(2.0)], [-1.0, 0.0, -1.0, 1.0], [-1.0, 1.0, -1.0, 1.0/sqrt(2.0)], [0.0, 1.0, -1.0, 1.0]]],
                               
                              [[[0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 1.0, 1.0],           [0.0, 0.0, 1.0, 1.0], [0.0, 0.0, 1.0, 1.0],            [0.0, 0.0, 1.0, 1.0],  [0.0, 0.0, 1.0, 1.0],             [0.0, 0.0, 1.0, 1.0],  [0.0, 0.0, 1.0, 1.0],            [0.0, 0.0, 1.0, 1.0]],
                               [[0.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0/sqrt(2.0)], [1.0, 0.0, 1.0, 1.0], [1.0, -1.0, 1.0, 1.0/sqrt(2.0)], [0.0, -1.0, 1.0, 1.0], [-1.0, -1.0, 1.0, 1.0/sqrt(2.0)], [-1.0, 0.0, 1.0, 1.0], [-1.0, 1.0, 1.0, 1.0/sqrt(2.0)], [0.0, 1.0, 1.0, 1.0]]],
                               
                              [[[0.0, 0.0, 3.0, 1.0], [0.0, 0.0, 3.0, 1.0],           [0.0, 0.0, 3.0, 1.0], [0.0, 0.0, 3.0, 1.0],            [0.0, 0.0, 3.0, 1.0],  [0.0, 0.0, 3.0, 1.0],             [0.0, 0.0, 3.0, 1.0],  [0.0, 0.0, 3.0, 1.0],            [0.0, 0.0, 3.0, 1.0]],
                               [[0.0, 1.0, 3.0, 1.0], [1.0, 1.0, 3.0, 1.0/sqrt(2.0)], [1.0, 0.0, 3.0, 1.0], [1.0, -1.0, 3.0, 1.0/sqrt(2.0)], [0.0, -1.0, 3.0, 1.0], [-1.0, -1.0, 3.0, 1.0/sqrt(2.0)], [-1.0, 0.0, 3.0, 1.0], [-1.0, 1.0, 3.0, 1.0/sqrt(2.0)], [0.0, 1.0, 3.0, 1.0]]]];
        
        
        int p = 3;
        int q = 1;
        int r = 2;
        double[] Z = [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0];
        double[] R = [0.0, 0.0, 1.0, 1.0];
        double[] THETA = [0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0];
        
        auto nvol = new NURBSVolume(Pw, Z, p, R, q, THETA, r);
        Vector3 P = Vector3(0.5, 0.0, 1.5);
        assert(approxEqualVectors(P, nvol(0.75, 0.5, 0.25)), failedUnitTest());

        return 0;
    }
}
