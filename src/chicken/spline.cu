// spline.cu
// A simple cubic spline for chicken.
//
// Peter J.
// 2022-11-30 Port from D and build on the template functions for simple linear algebra.
//

#ifndef CUBICSPLINE_INCLUDED
#define CUBICSPLINE_INCLUDED

#include <cmath>
#include <sstream>
#include <iostream>

#include "number.cu"
#include "rsla.cu"

using namespace std;

template<size_t n>
struct CubicSpline {
    number xs[n+1];
    number ys[n+1];
    number a[n];
    number b[n];
    number c[n];

    CubicSpline(const number _xs[n+1], const number _ys[n+1])
    // Sets up the interpolatory cubic spline through the (x, y) points.
    //   xs : sequence of x-coordinates
    //   ys : sequence of y-coordinates
    //
    // There will be n+1 knots, n segments and n-1 interior knots.
    {
        for (int i=0; i < n+1; ++i) {
            xs[i] = _xs[i];
            ys[i] = _ys[i];
        }
        // for d, just use the y array
        number h[n];
        for (int i=0; i < n; ++i) { h[i] = xs[i+1] - xs[i]; }
        // Solve for the second-derivative at the interior knots.
        // The sigma[i], i=1...n-1 (inclusive) are the unknowns since
        // the natural end conditions are sigma[0]=sigma[n]=0.0
        constexpr size_t m = n-1; // number of unknowns (and constraint equations)
        number A[m][m+1]; rsla::set_all_zero<m,m+1>(A); // Augmented matrix
        for (int k=0; k < m; ++k) {
            size_t i = k + 1; // i-index as in the lecture notes
            if (k > 0) { A[k][k-1] = h[i-1]; }
            A[k][k] = 2.0*(h[i-1]+h[i]);
            if (k < m-1) { A[k][k+1] = h[i]; }
            A[k][m] = 6.0*((ys[i+1]-ys[i])/h[i] - (ys[i]-ys[i-1])/h[i-1]); // RHS elements
        }
        if (rsla::gaussJordanElimination<m,m+1>(A)) {
            throw runtime_error("Gauss-Jordan elimination failed while solving for sigma.");
        }
        // Put sigma values in place in the full array.
        number sigma[n+1]; sigma[0] = 0.0; sigma[n] = 0.0;
        for (int i=0; i < m; ++i) { sigma[i+1] = A[i][m]; }
        // Evaluate and store coefficients for the polynomial segments.
        for (int i=0; i < n; ++i) {
            a[i] = (sigma[i+1]-sigma[i])/(6.0*h[i]);
            b[i] = sigma[i]/2.0;
            c[i] = (ys[i+1]-ys[i])/h[i] - (2.0*sigma[i] + sigma[i+1])*h[i]/6.0;
        }
    } // end constructor

    string toString() const
    {
        stringstream ss;
        ss << "CubicSpline(n=" << n;
        ss << ", xs=" << rsla::toString<n+1>(xs);
        ss << ", ys=" << rsla::toString<n+1>(ys);
        ss << ", a=" << rsla::toString<n>(a);
        ss << ", b=" << rsla::toString<n>(b);
        ss << ", c=" << rsla::toString<n>(c);
        ss << ")";
        return ss.str();
    }

    number eval(number x) const
    // Evaluates the spline at point x by first searching for the
    // relevant segment and then evaluating the local polynomial piece.
    {
        int i = 0;
        if (x < xs[0]) {
            i = 0; // The point is off range to the left.
        } else {
            i = n - 1; // Start search from the right-hand end.
            while ((x < xs[i]) && (i > 0)) { i -= 1; }
        }
        // Now that we found the relevant segment, evaluate locally.
        double dx = x - xs[i];
        return ((a[i]*dx + b[i])*dx + c[i])*dx + ys[i];
    }
}; // end CubicSpline

#endif
