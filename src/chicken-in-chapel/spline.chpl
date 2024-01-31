// spline.chpl
// A simple cubic spline class for Chicken-in-Chapel.
//
// Peter J.
// 2024-02-01 Port from the C++/CUDA variant (which came from the D variant,
//            which came from the Python code for my mech3750 lectures in 2011.
//            This Chapel variant is moving back toward the Python code.
//

module Spline {

  param nSpline = 4;  // number of cubic-polynomial segments
  // [TODO] Make this more like the C++ template so that the one program
  // can have several splines of various sizes.

  record CubicSpline {
    var xs, ys: [{0..nSpline}]real;    // n+1 interpolation points
    var a, b, c: [{0..<nSpline}]real;  // cubic-polynomial coefficients
    //
    // Sets up the interpolatory cubic spline through the (x, y) points.
    //   xs : sequence of x-coordinates
    //   ys : sequence of y-coordinates
    //
    // There will be n+1 knots, n segments and n-1 interior knots.
    proc ref set(xs: []real, ys: []real) throws
    {
      const n = nSpline;
      if n != xs.size-1 {
        writeln("Incorrect size for xs, ys; nSpline=", nSpline);
        throw new owned Error();
      }
      this.xs = xs;
      this.ys = ys;
      // for polynomial coefficients d, just use the ys array
      const DSeg = {0..#n};
      var h: [DSeg]real;
      for i in DSeg do h[i] = xs[i+1] - xs[i];
      // Solve for the second-derivative at the interior knots.
      // The sigma[i], i=1...n-1 (inclusive) are the unknowns since
      // the natural end conditions are sigma[0]=sigma[n]=0.0
      const m = n-1; // number of unknowns (and constraint equations)
      // Set up a linear system and solve.
      // Presently, the solver needs both A and rhs to be 2D arrays.
      var A: [{0..#m,0..#m}]real;
      var rhs: [{0..#m,0..0}]real;
      foreach k in {0..#m} {
        var i = k + 1; // i-index as in the lecture notes
        if k > 0 then A[k,k-1] = h[i-1];
        A[k,k] = 2.0*(h[i-1]+h[i]);
        if k < m-1 then A[k,k+1] = h[i];
        rhs[k,0] = 6.0*((ys[i+1]-ys[i])/h[i] - (ys[i]-ys[i-1])/h[i-1]);
      }
      if gaussJordanElimination(A,rhs) != 0 {
        writeln("Gauss-Jordan elimination failed while solving for sigma.");
        throw new owned Error();
      }
      // Put sigma values in place in the full array.
      var sigma: [{0..#(n+1)}]real;
      sigma[0] = 0.0; sigma[n] = 0.0;
      for i in {0..#m} do sigma[i+1] = rhs[i,0];
      // Evaluate and store coefficients for the polynomial segments.
      foreach i in {0..#n} {
        a[i] = (sigma[i+1]-sigma[i])/(6.0*h[i]);
        b[i] = sigma[i]/2.0;
        c[i] = (ys[i+1]-ys[i])/h[i] - (2.0*sigma[i] + sigma[i+1])*h[i]/6.0;
      }
    } // end proc set()

    proc ref eval(x: real): real
    {
      var i = 0;
      if x < xs[0] then
        i = 0; // The point is off range to the left.
      else {
        i = nSpline-1; // Start the search from the right-hand end.
        while x < xs[i] && i > 0 do i -= 1;
      }
      // Now that the have found the relevant segment,
      // evaluate the local polynomial.
      var dx = x - xs[i];
      return ((a[i]*dx + b[i])*dx + c[i])*dx + ys[i];
    }

    proc ref x0(): real { return xs[0]; }
    proc ref y0(): real { return ys[0]; }
    proc ref xn(): real { return xs[nSpline]; }
    proc ref yn(): real { return ys[nSpline]; }
  } // end record CubicSpline

} // end module Spline
