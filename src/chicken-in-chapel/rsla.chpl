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
        const cRows = c.dim(0);
        const cCols = c.dim(1);
        forall i in cRows {
          result[i] = 0.0;
          forall j in cCols do result[i] += c[i,j]*x[j];
        }
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
        var A = c; // Make a copy so that we don't mess with c.
        cinv = 0.0;
        for i in cinv.dim(0) do cinv[i,i] = 1.0;
        return gaussJordanElimination(A, cinv, very_small_value);
      }
    } // end select
  } // end MInverse

  // Perform Gauss-Jordan elimination on an augmented matrix [A|b]
  // such that the mutated matrix becomes [I|x]
  // where x is the solution vector(s) to A.x = b
  // Returns  0 on success
  //         -1 if matrix A is singular
  //         -2 if A and b are not compatible
  //         -3 if A is not square
  // Note that we expect b to be a 2-dimensional array.
  proc gaussJordanElimination(ref A: []real, ref b: []real,
                              very_small_value: real = 1.0e-16): int
  {
    var (rowsA, colsA) = A.dims();
    var (rowsB, colsB) = b.dims();
    if colsA != rowsB then return -2; // A and b not compatible
    if rowsA != colsA then return -3; // A not square
    for j in rowsA {
      // Select pivot.
      var p = j;
      for i in {(j+1)..rowsA.high} {
        if abs(A[i,j]) > abs(A[p,j]) then p = i;
      }
      if (abs(A[p,j]) < very_small_value) {
        return -1; // Matrix is essentially singular.
      }
      if (p != j) {
        forall k in colsA do (A[p,k], A[j,k]) = (A[j,k], A[p,k]);
        forall k in colsB do (b[p,k], b[j,k]) = (b[j,k], b[p,k]);
      }
      // Scale row j to get unity on the diagonal.
      const Ajj = A[j,j];
      forall k in colsA do A[j,k] /= Ajj;
      forall k in colsB do b[j,k] /= Ajj;
      // Do the elimination to get zeros in all off diagonal values in column j.
      for i in rowsA {
        if i != j {
          const Aij = A[i,j];
          for k in colsA do A[i,k] -= Aij * A[j,k];
          for k in colsB do b[i,k] -= Aij * b[j,k];
        }
      }
    } // end for j
    return 0; // success
  } // end proc gaussJordanElimination()

} // end module Rsla
