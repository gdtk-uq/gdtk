/*
 * Authors: Rowan G. and Peter J.
 * Date: 2019-04-02
 *
 * The Bezier triangle patch is not part of the ParametricSurface
 * family. Although the Bezier triangle is bivariate (in u and v),
 * the u,v ranges do not *independently* extend from 0 <= u,v <= 1.
 * Rather the u,v parameters are part of a triple (u, v, w) that
 * form barycentric coordinates. These have the constraint:
 *    u + v + w = 1.
 * We can treat the surface as bivariate because given u and v,
 * then w is implied:
 *    w = 1 - u - v
 *
 * We provide this Bezier triangle patch because it is useful as
 * a construction patch and for parameterisation of surfaces for
 * design optimization.
 *
 * The implementation here follows closesly the description in
 * Farin (2001) and Hansford (2002).
 *
 * References
 * ----------
 * Farin (2001)
 * Curves and Surfaces for CAGD: A Practical Guide.
 * Elsevier Science & Technology
 *
 * Hansford (2002)
 * Bezier Techniques, Ch.4 in Handbook of Computer Aided Geometric Design.
 * Editors: Farin et al.
 * Elsevier Science & Technology
 *
 */

module geom.surface.beziertrianglepatch;

import std.conv : text;
import std.range : iota;

import geom.elements;

class BezierTrianglePatch {
public:
    Vector3[] B; // control points
    int n; // degree of patch

    this(const Vector3[] _B, int _n)
    {
        n = _n;
        auto nCtrlPts = (n+1)*(n+2)/2;
        if (nCtrlPts != _B.length) {
            throw new Error(text("BezierTrianglePatch: incorrect number of control points supplied.",
                                 nCtrlPts, " expected, ", _B.length, " supplied."));
        }
        B.length = nCtrlPts;
        B[] = _B[];
        initialiseWorkingSpace();
    }

    void initialiseWorkingSpace()
    {
        _b.length = n + 1;
        foreach (i; 0 .. n + 1) {
            auto n_local = n - i;
            _b[i].length = (n_local + 1)*(n_local + 2) / 2;
        }
        // Copy B into _b[0].
        // Yes, we are duplicating the storage of those control points,
        // but it makes the implementation of de Casteljau's algorithm
        // much nicer. Typically, the number of control points will
        // be relatively small, so this shouldn't be a big storage penalty.
        _b[0][] = B[];
    }

    @nogc size_t toSingleIndex()(size_t i, size_t j, size_t k) const
    {
        auto n = i + j + k;
        auto m = n - i;
        return (m*(m+5)/2) - 2*j - k;
    }

    Vector3 opCall(double u, double v)
    {
        auto w = 1.0 - u - v;
        // We may allow u or v to get slightly larger than 1.0
        // if we are perturbing to get sensitivity.
        // In which case, we won't allow w to go negative.
        if (w < 0.0) w = 0.0;

        // de Casteljau algorithm as described in:
        // Farin (2001), pp. 310--311
        foreach (r; 1 .. n + 1) {
            // Loop over points via triangle indices.
            auto n_local = n - r;
            foreach (i; iota(n_local, -1, -1)) {
                foreach (j; iota(n_local-i, -1, -1)) {
                    auto k = n_local - (i+j);
                    auto idx = toSingleIndex!()(i, j, k);
                    auto uIdx = toSingleIndex!()(i+1, j, k);
                    auto vIdx = toSingleIndex!()(i, j+1, k);
                    auto wIdx = toSingleIndex!()(i, j, k+1);
                    // Do Vector3 operation component-wise to
                    // avoid allocation of temporaries.
                    _b[r][idx].refx = u*_b[r-1][uIdx].x + v*_b[r-1][vIdx].x + w*_b[r-1][wIdx].x;
                    _b[r][idx].refy = u*_b[r-1][uIdx].y + v*_b[r-1][vIdx].y + w*_b[r-1][wIdx].y;
                    _b[r][idx].refz = u*_b[r-1][uIdx].z + v*_b[r-1][vIdx].z + w*_b[r-1][wIdx].z;
                }
            }
        }
        return _b[n][0];
    }

private:
    Vector3[][] _b; // working space for de Casteljau algorithm

} // end class BezierTrianglePatch


version(beziertrianglepatch_test) {
    import util.msg_service;
    import std.stdio;
    int main()
    {
        // Example on 17.1 on p. 312 in Farin (2001)
        Vector3[] B = [Vector3(6, 0, 9),
                       Vector3(3, 3, 6), Vector3(3, 0, 0),
                       Vector3(0, 6, 0), Vector3(0, 3, 0), Vector3(0, 0, 0)];
        auto btp = new BezierTrianglePatch(B, 2);
        double u = 1./3;
        double v = 1./3;
        auto p = btp(u, v);
        // Farin gives results as (2, 2, 7/3);
        assert(approxEqualVectors(p, Vector3(2, 2, 7./3)), failedUnitTest());
        return 0;
    }

}
