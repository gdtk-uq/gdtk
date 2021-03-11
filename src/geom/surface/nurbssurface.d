/* nurbssurface.d
Author: Reece O. 
Date: 2021-03-09
*/

module geom.surface.nurbssurface;

import std.conv;
import std.stdio;

import geom.elements;
import geom.surface.parametricsurface;
import nurbs_utils;

class NURBSSurface : ParametricSurface {
public:
    const double[][][] Pw; // array of weighted control points
    
    const double[] U; // knot vector in direction 1
    int p; // degree of curve in direction 1
    int a; // a+1 = number of control points in direction 1
    
    const double[] V; // knot vector in direction 2
    int q; // degree of curve in direction 2
    int b; // b+1 = number of control points in direction 2
    
    this(in double[][][] Pw, in double[] U, in int p, in double[] V, in int q) 
    {
        this.Pw = Pw; 
        this.U = U; this.p = p; this.a = to!int(U.length-p-2);
        this.V = V; this.q = q; this.b = to!int(V.length-q-2); 
    }
    this(ref const(NURBSSurface) other)
    {
        Pw = other.Pw; 
        U = other.U; p = other.p; a = other.a;
        V = other.V; q = other.q; b = other.b;
    }
    override NURBSSurface dup() const
    {
        return new NURBSSurface(Pw, U, p, V, q);
    }
    override Vector3 opCall(double u, double v) const
    {
        return deBoor(u, v);
    }
    override string toString() const
    {
        return "NURBSSurface(Pw=" ~ to!string(Pw) ~ ", U=" ~ to!string(U) ~ ", p=" ~ to!string(p) ~ ", V=" ~ to!string(V) ~ ", q=" ~ to!string(q) ~ ")";
    }

protected:
    Vector3 deBoor(double u, double v) const {
        // Returns the Cartesian coordinates of a point on a NURBS surface at a given set of parameter values
        // This is algorithm A4.3 from Piegl and Tiller (1997) - 'The NURBS Book'
    
        int _uspan = FindSpan(u, a, p, U);
        double[] _Nu = BasisFuns(_uspan, u, p, U);
        int _vspan = FindSpan(v, b, q, V);
        double[] _Nv = BasisFuns(_vspan, v, q, V);
        auto _temp = new double[4][q+1]; 
        foreach (l; 0 .. q+1) {
            _temp[l][] = [0.0, 0.0, 0.0, 0.0];
            foreach (k; 0 .. p+1) _temp[l][] += _Nu[k]*Pw[_uspan-p+k][_vspan-q+l][];
        }
        double[] _Sw = [0.0, 0.0, 0.0, 0.0];
        double[] _S = [0.0, 0.0, 0.0];
        foreach (l; 0 .. q+1) _Sw[] += _Nv[l]*_temp[l][];
        _S[] = _Sw[]/_Sw[_Sw.length-1];
        Vector3 _Svec = Vector3(_S);
        return _Svec;
    }
}

version(nurbssurface_test) {
    import util.msg_service;
    int main () {
    // This is example Ex4.3 from Piegl and Tiller (1997) - 'The NURBS Book'
    // Example only included 'activated' control points, so just set all others to 0 since they have no contribution
    double[][][] Pw = [[[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]],
                       [[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 4.0, 1.0], [4.0, 6.0, 8.0, 2.0], [4.0, 2.0, 4.0, 1.0]], 
                       [[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 6.0, 4.0, 2.0], [12.0, 24.0, 12.0, 6.0], [8.0, 6.0, 4.0, 2.0]],
                       [[0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 0.0, 1.0], [4.0, 6.0, 0.0, 2.0], [4.0, 2.0, 0.0, 1.0]]];
                   
    int p = 2;
    int q = 2;
    double[] U = [0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 3.0, 3.0];
    double[] V = [0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 5.0, 5.0, 5.0];
    auto nsurf = new NURBSSurface(Pw, U, p, V, q);
    assert(approxEqualVectors(Vector3(2.0, 98.0/27.0, 68.0/27.0), nsurf.opCall(1.0, 5.0/2.0)), failedUnitTest());
    return 0;
    }
}
