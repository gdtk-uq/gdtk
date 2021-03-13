/* nurbscurve.d
 * Author: Reece O. 
 * Date: 2021-02-24
 */

module geom.gpath.nurbscurve;

import std.conv;

import geom.elements;
import geom.gpath.path;
import nurbs_utils;

class NURBSCurve : Path {
public:
    const double[4][] Pw; // collection of weighted control points
    const double[] U; // knot vector
    int p; // degree of curve
    int a; // a+1 = number of control points
    
    this(in double[4][] Pw, in double[] U, in int p) 
    {
        this.Pw = Pw; this.U = U; this.p = p; this.a = to!int(U.length-p-2);   
    }
    this(ref const(NURBSCurve) other)
    {
        Pw = other.Pw; U = other.U; p = other.p; a = other.a;
    }
    override NURBSCurve dup() const
    {
        return new NURBSCurve(Pw, U, p);
    }
    override Vector3 opCall(double u) const
    {
        return deBoor(u);
    }
    override string toString() const
    {
        return "NURBS(Pw=" ~ to!string(Pw) ~ ", U=" ~ to!string(U) ~ ", p=" ~ to!string(p) ~ ")";
    }
    override string classString() const
    {
        return "NURBS";
    }
    
protected:
    Vector3 deBoor(double u) const {
        // Returns the Cartesian coordinates of a point on a NURBS curve at a given parameter value
        // This is algorithm A4.1 from Piegl and Tiller (1997) - 'The NURBS Book'
    
        int _span = FindSpan(u, a, p, U);
        double[] _N = BasisFuns(_span, u, p, U);
        double[] _Cw = [0.0, 0.0, 0.0, 0.0];
        double[] _C = [0.0, 0.0, 0.0];
        foreach (i; 0 .. p+1) _Cw[] += _N[i]*Pw[_span-p+i][];
        _C[] = _Cw[]/_Cw[_Cw.length-1];
        Vector3 _Cvec = Vector3(_C);
        return _Cvec;
    }
    
}

version(nurbscurve_test) {
    import util.msg_service;
    int main() {
        double[][] Pw = [[-4.0, -4.0, 0.0, 1.0], [-2.0, 4.0, 0.0, 1.0], [2.0*5, -4.0*5, 0.0, 5], [4.0, 4.0, 0.0, 1.0], [3.778, 1.836, 2.933, 1.0], [4*2.772, -4*3.875, 4*1.736, 4.0]];
        int p = 2;
        double[] U = [0.0, 0.0, 0.0, 0.375, 0.5, 0.625, 1.0, 1.0, 1.0];
        auto ncurve = new NURBSCurve(Pw, U, p);
        assert(approxEqualVectors(Vector3(3.782, 2.939, 0.435), ncurve.opCall(0.6)), failedUnitTest());
        return 0;
    }
}

