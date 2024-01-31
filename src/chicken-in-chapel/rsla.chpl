// rsla.chpl
// Real simple linear algebra for the Chicken-in-Chapel flow solver.
//
// A linear equation solver using the text-book notation.
// Calculation of the flow gradients (for the viscous terms) requires
// the solution of many sets of 3x3 linear equations and the spline
// calculation requires a general linear solver.
//
// PJ 2024-01-31 Model the capability of the CUDA/C++ version in Chicken.

module Rsla {

  use Math;

  const D2x2: domain(2) = {0..#2, 0..#2};
  const D2x1: domain(1) = {0..#2};
  const D3x3: domain(2) = {0..#3, 0..#3};
  const D3x1: domain(1) = {0..#3};

  proc MVMult(const ref c: []real, const ref x: []real, ref result: []real)
  {
    select (c.domain, x.domain) {
      when (D2x2, D2x1) {
        result[0] = c[0,0]*x[0] + c[0,1]*x[1];
        result[1] = c[1,0]*x[0] + c[1,1]*x[1];
        return;
      }
      when (D3x3, D3x1) {
        result[0] = c[0,0]*x[0] + c[0,1]*x[1] + c[0,2]*x[2];
        result[1] = c[1,0]*x[0] + c[1,1]*x[1] + c[1,2]*x[2];
        result[2] = c[2,0]*x[0] + c[2,1]*x[1] + c[2,2]*x[2];
        return;
      }
      otherwise {
        // [TODO]
        return;
      }
    }
  } // end proc MVMult


  proc MInverse(const ref c: []real, ref cinv: []real,
                const very_small_value: real=1.0e-12): int
  {
    select c.domain {
      when D2x2 {
        // Compute inverse directly
        const det = c[0,0]*c[1,1] - c[0,1]*c[1,0];
        if abs(det) <= very_small_value then return -1; // singular
        const one_over_det = 1.0/det;
        cinv[0,0] =  c[1,1]*one_over_det; cinv[0,1] = -c[0,1]*one_over_det;
        cinv[1,0] = -c[1,0]*one_over_det; cinv[1,1] =  c[0,0]*one_over_det;
        return 0;
      }
      when D3x3 {
        // Compute inverse directly
        const det = c[0,0]*(c[1,1]*c[2,2] - c[1,2]*c[2,1])
              - c[0,1]*(c[1,0]*c[2,2] - c[1,2]*c[2,0])
              + c[0,2]*(c[1,0]*c[2,1] - c[1,1]*c[2,0]);
        if abs(det) <= very_small_value then return -1; // singular
        const one_over_det = 1.0/det;
        cinv[0,0] = (c[1,1]*c[2,2] - c[1,2]*c[2,1])*one_over_det;
        cinv[0,1] = (c[0,2]*c[2,1] - c[0,1]*c[2,2])*one_over_det;
        cinv[0,2] = (c[0,1]*c[1,2] - c[0,2]*c[1,1])*one_over_det;
        cinv[1,0] = (c[1,2]*c[2,0] - c[1,0]*c[2,2])*one_over_det;
        cinv[1,1] = (c[0,0]*c[2,2] - c[0,2]*c[2,0])*one_over_det;
        cinv[1,2] = (c[0,2]*c[1,0] - c[0,0]*c[1,2])*one_over_det;
        cinv[2,0] = (c[1,0]*c[2,1] - c[1,1]*c[2,0])*one_over_det;
        cinv[2,1] = (c[0,1]*c[2,0] - c[0,0]*c[2,1])*one_over_det;
        cinv[2,2] = (c[0,0]*c[1,1] - c[0,1]*c[1,0])*one_over_det;
        return 0;
      }
      otherwise {
        // [TODO] using the GaussJordan elimination function
        return -2; // not implemented
      }
    } // end select
  } // end MInverse

} // end module Rsla
