// matrix.d
module linalg;

import std.stdio;
import std.math;
import std.conv;
import std.file;
import std.format;
import std.string;
import number;

void prep_matrix(ref double[][] a, size_t ncols, size_t nrows) {
    a.length = nrows;
    foreach (ref row; a) {
        row.length = ncols;
        row[] = 0.0;
    }
} // end prep_matrix()

void prep_vector(ref double[] a, size_t len) {
    a.length = len;
    a[] = 0.0;
} // end prep_vector()

void matrix_multiply(double[][] a, double[][] b, ref double[][] c) {

    size_t a_nrows = a.length;
    size_t a_ncols = a[0].length;

    size_t b_nrows = b[0].length;
    size_t b_ncols = b.length;

    assert( (a_ncols == b_nrows), "incompatible matrices for matrix multiplication.");
    
    for (int i = 0; i < a_nrows; i++) {
        for (int j = 0; j < b_ncols; j++) {
            c[i][j] = 0;
            for (int k = 0; k < a_ncols; k++) {
                c[i][j] += a[i][k]*b[k][j];
	    }
	}
    }
    
} // end matrix_multiply()

void matrix_vector_multiply(double[][] a, double[] b, ref double[] c) {

    size_t a_nrows = a.length;
    size_t a_ncols = a[0].length;

    size_t b_len = b.length;

    assert( (a_ncols == b_len), "incompatible dimensions for matrix-vector multiplication.");
    c[] = 0.0;
    for (int i = 0; i < a_nrows; i++) {
        for (int j = 0; j < b_len; j++) {
            c[i] += a[i][j] * b[j];
	}
    }
    
} // end matrix_vector_multiply()


void matrix_inverse(ref double[][] a) {

    /*
      Simple Gauss-Jordan elimination routine that inverts a given matrix
      algorithm from:
          Numerical Recipes, Press et al., 3rd eition (2007), pg. 44
    */

    int i, icol, irow, j, k, l, ll;
    int n = to!int(a.length);
    double big, dum, pivinv;
    int[] indxc, indxr, ipiv;
    indxc.length = n;
    indxr.length = n;
    ipiv.length = n;
    
    for (j=0; j<n; j++) { ipiv[j] = 0; }
    for (i=0; i<n; i++) {
        big = 0.0;
        for (j=0; j<n; j++) {
            if (ipiv[j] != 1) {
                for (k=0; k<n; k++) {
                    if (ipiv[k] == 0) {
                        if (abs(a[j][k]) >= big) {
                            big = abs(a[j][k]);
                            irow = j;
                            icol = k;
                        }
                    }
                }
            }
        }
        ++(ipiv[icol]);
        if (irow != icol) {
            for (l=0; l<n; l++) {
                double tmp = a[irow][l];
                a[irow][l] = a[icol][l];
                a[icol][l] = tmp;
            }
        }
        indxr[i] = irow;
        indxc[i] = icol;
        assert( (a[icol][icol] != 0.0), "gaussj: Singular Matrix");
        pivinv = 1.0/a[icol][icol];
        a[icol][icol]=1.0;
        for (l=0; l<n; l++) { a[icol][l] *= pivinv; }
        for (ll=0; ll<n; ll++) {
            if (ll != icol) {
                dum = a[ll][icol];
                a[ll][icol]=0.0;
                for (l=0; l<n; l++) { a[ll][l] -= a[icol][l]*dum; }
            }
        }
    }
    for (l=n-1; l>=0; l--) {
        if (indxr[l] != indxc[l]) {
            for (k=0; k<n; k++) {
                double tmp = a[k][indxr[l]];
                a[k][indxr[l]] = a[k][indxc[l]];
                a[k][indxc[l]] = tmp;
            }
        }
    }
    
} // end matrix_inverse()

version(matrix_test) {
    import std.math;
    int main() {

        // inverse of a matrix
        double[][] a;
        int n = 3;
        prep_matrix(a, n, n);
        a[0][0] = 3.; a[0][1] = 0.; a[0][2] = 2.;
        a[1][0] = 2.; a[1][1] = 0.; a[1][2] = -2.;
        a[2][0] = 0.; a[2][1] = 1.; a[2][2] = 1.;

        matrix_inverse(a);
        
        double[][] sol;
        prep_matrix(sol, n, n);
        sol[0][0] = 0.2;  sol[0][1] = 0.2;  sol[0][2] = 0.;
        sol[1][0] = -0.2; sol[1][1] = 0.3;  sol[1][2] = 1.;
        sol[2][0] = 0.2;  sol[2][1] = -0.3; sol[2][2] = 0.;
        
        for (int i=0; i < n; i++) { assert(approxEqualNumbers(a[i], sol[i]), failedUnitTest()); }

        // matrix-matrix multiplication
        a[0][0] = 3.; a[0][1] = 0.; a[0][2] = 2.;
        a[1][0] = 2.; a[1][1] = 0.; a[1][2] = -2.;
        a[2][0] = 0.; a[2][1] = 1.; a[2][2] = 1.;

        double[][] b;
        prep_matrix(b, n, n);
        b[0][0] = 4.; b[0][1] = 7.; b[0][2] = -4.;
        b[1][0] = 8.; b[1][1] = 2.; b[1][2] = -3.;
        b[2][0] = 2.; b[2][1] = 1.; b[2][2] = 9.;

        double[][] c;
        prep_matrix(c, n, n);

        matrix_multiply(a, b, c);
        sol[0][0] = 16.; sol[0][1] = 23.; sol[0][2] = 6.;
        sol[1][0] = 4.;  sol[1][1] = 12.; sol[1][2] = -26.;
        sol[2][0] = 10.; sol[2][1] = 3.;  sol[2][2] = 6.;

        for (int i=0; i < n; i++) { assert(approxEqualNumbers(c[i], sol[i]), failedUnitTest()); }

        // matrix-vector multiplication
        double[] d;
        prep_vector(d, n);
        d[0] = 12.;
        d[1] = 16.;
        d[2] = -20.;

        double[] e;
        prep_vector(e, n);

        double[] sol_vec;
        prep_vector(sol_vec, n);
        sol_vec[0] = -4.;
        sol_vec[1] = 64.;
        sol_vec[2] = -4.;
        
        matrix_vector_multiply(a, d, e);
        
        assert(approxEqualNumbers(e, sol_vec), failedUnitTest());
        
        
        return 0;
    }
}
