/**
 * rsla: real simple linear algebra
 * A linear equation solver using the text-book notation.
 *
 * PJ 2016-01-16 factored out of the gradient calculation module for eilmer
 */

module nm.rsla;

import std.math;
import std.format;
import std.conv;
import std.stdio;
import ntypes.complex;

@nogc
T normInf(size_t N, size_t NDIM, size_t NDIM2, T)(ref T[NDIM2][NDIM] c)
// Return the max-absolute-row-sum of an augmented matrix c = [A|b]
// We are concerned with a measure of A only.
{
    T norm = 0.0;
    foreach (i; 0 .. N) {
        T rowsum = 0.0;
        foreach (j; 0 .. N) { rowsum += fabs(c[i][j]); }
        norm = fmax(rowsum, norm);
    }
    return norm;
} // end normInf()()

@nogc
int computeInverse(size_t N, size_t NDIM, size_t NDIM2, T)
    (ref T[NDIM2][NDIM] c, double very_small_value=1.0e-16)
// Perform Gauss-Jordan elimination on an augmented matrix.
// c = [A|b] such that the mutated matrix becomes [I|x]
// where x is the solution vector(s) to A.x = b
// When computing an inverse, the incoming data is assumed to be c=[A|I].
//
// Returns 0 normally, -1 if the matrix is essentially singular (zero pivot).
//
// We allow a larger storage than needed because the flow solver code will
// use these functions for both 2D and 3D simulations.
// When N < NDIM, the solver just uses a subset of the available storage
// and the rest remains untouched.
{
    assert(NDIM >= N, "Inadequate size of dimension for matrix");
    static if (N == 0) {
        assert(false, "Zero dimension linear system doesn't make sense.");
    }
    static if (N == 1) {
        T det = c[0][0];
        if (abs(det) <= very_small_value) return -1; // singular
        c[0][1] = 1.0/det; // inverse
        c[0][0] = 1.0; // identity
    }
    static if (N == 2) {
        T det = c[0][0]*c[1][1] - c[0][1]*c[1][0];
        if (abs(det) <= very_small_value) return -1; // singular
        // compute inverse directly
        T one_over_det = 1.0/det;
        c[0][2] =  c[1][1]*one_over_det; c[0][3] = -c[0][1]*one_over_det;
        c[1][2] = -c[1][0]*one_over_det; c[1][3] =  c[0][0]*one_over_det;
        // overwrite original elements with identity
        c[0][0] = 1.0; c[0][1] = 0.0;
        c[1][0] = 0.0; c[1][1] = 1.0;
    }
    static if (N == 3) {
        T det = c[0][0]*(c[1][1]*c[2][2] - c[1][2]*c[2][1])
            - c[0][1]*(c[1][0]*c[2][2] - c[1][2]*c[2][0])
            + c[0][2]*(c[1][0]*c[2][1] - c[1][1]*c[2][0]);
        if (abs(det) <= very_small_value) return -1; // singular
        // compute inverse directly
        T one_over_det = 1.0/det;
        c[0][3] = (c[1][1]*c[2][2] - c[1][2]*c[2][1])*one_over_det;
        c[0][4] = (c[0][2]*c[2][1] - c[0][1]*c[2][2])*one_over_det;
        c[0][5] = (c[0][1]*c[1][2] - c[0][2]*c[1][1])*one_over_det;
        c[1][3] = (c[1][2]*c[2][0] - c[1][0]*c[2][2])*one_over_det;
        c[1][4] = (c[0][0]*c[2][2] - c[0][2]*c[2][0])*one_over_det;
        c[1][5] = (c[0][2]*c[1][0] - c[0][0]*c[1][2])*one_over_det;
        c[2][3] = (c[1][0]*c[2][1] - c[1][1]*c[2][0])*one_over_det;
        c[2][4] = (c[0][1]*c[2][0] - c[0][0]*c[2][1])*one_over_det;
        c[2][5] = (c[0][0]*c[1][1] - c[0][1]*c[1][0])*one_over_det;
        // overwrite original elements with identity
        c[0][0] = 1.0; c[0][1] = 0.0; c[0][2] = 0.0;
        c[1][0] = 0.0; c[1][1] = 1.0; c[1][2] = 0.0;
        c[2][0] = 0.0; c[2][1] = 0.0; c[2][2] = 1.0;
    }
    static if (N > 3) {
        foreach(j; 0 .. N) {
            // Select pivot.
            size_t p = j;
            foreach(i; j+1 .. N) {
                if ( abs(c[i][j]) > abs(c[p][j]) ) p = i;
            }
            if (abs(c[p][j]) <= very_small_value) return -1; // singular
            if ( p != j ) { // Swap rows
                foreach(col; 0 .. 2*N) {
                    T tmp = c[p][col]; c[p][col] = c[j][col]; c[j][col] = tmp;
                }
            }
            // Scale row j to get unity on the diagonal.
            T cjj = c[j][j];
            foreach(col; 0 .. 2*N) c[j][col] /= cjj;
            // Do the elimination to get zeros in all off diagonal values in column j.
            foreach(i; 0 .. N) {
                if ( i == j ) continue;
                T cij = c[i][j];
                foreach(col; 0 .. 2*N) c[i][col] -= cij * c[j][col];
            }
        } // end foreach j
    } // end static if N > 3
    return 0; // success
} // end computeInverse()()

int computeInverseDebug(size_t N, size_t NDIM, size_t NDIM2, T)
    (ref T[NDIM2][NDIM] c, double very_small_value=1.0e-16)
// Perform Gauss-Jordan elimination on an augmented matrix.
// c = [A|b] such that the mutated matrix becomes [I|x]
// where x is the solution vector(s) to A.x = b
// When computing an inverse, the incoming data is assumed to be c=[A|I].
//
// Returns 0 normally, -1 if the matrix is essentially singular (zero pivot).
//
// This debug version gives more information.
{
    assert(NDIM >= N, "Inadequate size of dimension for matrix");
    T[2*N][N] csave;
    foreach(j; 0 .. N) {
        foreach(col; 0 .. 2*N) csave[j][col] = c[j][col];
    }
    writeln("At start of computeInverseDebug, c=", c);
    foreach(j; 0 .. N) {
        // Select pivot.
        size_t p = j;
        foreach(i; j+1 .. N) {
            if ( abs(c[i][j]) > abs(c[p][j]) ) p = i;
        }
        if (abs(c[p][j]) <= very_small_value) {
            string msg = format(" matrix is essentially singular j=%d p=%d", j, p) ~
                " \nc=" ~to!string(c) ~ " \ncsave=" ~ to!string(csave);
            writeln(msg);
            return -1; // singular
        }
        if ( p != j ) { // Swap rows
            foreach(col; 0 .. 2*N) {
                T tmp = c[p][col]; c[p][col] = c[j][col]; c[j][col] = tmp;
            }
        }
        // Scale row j to get unity on the diagonal.
        T cjj = c[j][j];
        foreach(col; 0 .. 2*N) c[j][col] /= cjj;
        // Do the elimination to get zeros in all off diagonal values in column j.
        foreach(i; 0 .. N) {
            if ( i == j ) continue;
            T cij = c[i][j];
            foreach(col; 0 .. 2*N) c[i][col] -= cij * c[j][col];
        }
    } // end foreach j
    return 0; // success
} // end computeInverseDebug()()

@nogc
void solveWithInverse(size_t N, size_t NDIM, size_t NDIM2, T)
    (ref T[NDIM2][NDIM] c, ref T[NDIM] rhs, ref T[NDIM] x)
// Multiply right-hand-side by the inverse part of the augmented matrix.
// Augmented matrix is assumed to be c=[I|Ainv]
{
    assert(NDIM >= N, "Inadequate size of dimension for matrix");
    static if (N == 0) {
        assert(false, "Zero dimension linear system doesn't make sense.");
    }
    static if (N == 1) {
        x[0] = c[0][1] * rhs[0];
    }
    static if (N == 2) {
        x[0] = c[0][2]*rhs[0] + c[0][3]*rhs[1];
        x[1] = c[1][2]*rhs[0] + c[1][3]*rhs[1];
    }
    static if (N == 3) {
        x[0] = c[0][3]*rhs[0] + c[0][4]*rhs[1] + c[0][5]*rhs[2];
        x[1] = c[1][3]*rhs[0] + c[1][4]*rhs[1] + c[1][5]*rhs[2];
        x[2] = c[2][3]*rhs[0] + c[2][4]*rhs[1] + c[2][5]*rhs[2];
    }
    static if (N > 3) {
        foreach(i; 0 .. N) {
            x[i] = 0.0;
            foreach(j; 0 .. N) {
                x[i] += c[i][N+j] * rhs[j];
            }
        }
    } // end static if N > 3
} // end solveWithInverse()()


version(rsla_test) {
    import util.msg_service;
    import nm.number;
    int main() {
        number[8][4] A = [[to!number(0.0),  to!number(2.0),  to!number(0.0),  to!number(1.0),
                           to!number(1.0), to!number(0.0), to!number(0.0), to!number(0.0)],
                          [to!number(2.0),  to!number(2.0),  to!number(3.0),  to!number(2.0),
                           to!number(0.0), to!number(1.0), to!number(0.0), to!number(0.0)],
                          [to!number(4.0), to!number(-3.0),  to!number(0.0),  to!number(1.0),
                           to!number(0.0), to!number(0.0), to!number(1.0), to!number(0.0)],
                          [to!number(6.0),  to!number(1.0), to!number(-6.0), to!number(-5.0),
                           to!number(0.0), to!number(0.0), to!number(0.0), to!number(1.0)]];
        assert(approxEqualNumbers(normInf!(4,4,8,number)(A), to!number(18.0)), failedUnitTest());
        computeInverse!(4,4,8,number)(A);
        number[4] b = [to!number(0.0), to!number(-2.0), to!number(-7.0), to!number(6.0)];
        number[4] x;
        solveWithInverse!(4,4,8,number)(A, b, x);
        assert(approxEqualNumbers(x[0], to!number(-0.5)) &&
               approxEqualNumbers(x[1], to!number(1.0)) &&
               approxEqualNumbers(x[2], to!number(1.0/3)) &&
               approxEqualNumbers(x[3], to!number(-2.0)),
               failedUnitTest());

        // Try same workspace with a smaller 2x2 system.
        x[0] = -0.5; x[1] = 1.0;
        A[0][0] = 0.0; A[0][1] = 2.0; A[0][2] = 1.0; A[0][3] = 0.0;
        A[1][0] = 2.0; A[1][1] = 2.0; A[1][2] = 0.0; A[1][3] = 1.0;
        b[0] = A[0][0]*x[0] + A[0][1]*x[1];
        b[1] = A[1][0]*x[0] + A[1][1]*x[1];
        assert(approxEqualNumbers(normInf!(2,4,8,number)(A), to!number(4.0)),
               failedUnitTest());
        computeInverse!(2,4,8,number)(A);
        x[0] = 0.0; x[1] = 0.0;
        solveWithInverse!(2,4,8,number)(A, b, x);
        assert(approxEqualNumbers(x[0], to!number(-0.5)) &&
               approxEqualNumbers(x[1], to!number(1.0)),
               failedUnitTest());

        // and again, with a 3x3 system.
        x[0] = -0.5; x[1] = 1.0; x[2] = 1.0/3;
        A[0][0] = 0.0; A[0][1] = 2.0; A[0][2] = 0.0;
        A[1][0] = 2.0; A[1][1] = 2.0; A[1][2] = 3.0;
        A[2][0] = 4.0; A[2][1] = -3.0; A[2][2] = 0.0;
        A[0][3] = 1.0; A[0][4] = 0.0; A[0][5] = 0.0;
        A[1][3] = 0.0; A[1][4] = 1.0; A[1][5] = 0.0;
        A[2][3] = 0.0; A[2][4] = 0.0; A[2][5] = 1.0;
        b[0] = A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2];
        b[1] = A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2];
        b[2] = A[2][0]*x[0] + A[2][1]*x[1] + A[2][2]*x[2];
        assert(approxEqualNumbers(normInf!(3,4,8,number)(A), to!number(7.0)),
               failedUnitTest());
        computeInverse!(3,4,8,number)(A);
        x[0] = 0.0; x[1] = 0.0; x[2] = 0.0;
        solveWithInverse!(3,4,8,number)(A, b, x);
        assert(approxEqualNumbers(x[0], to!number(-0.5)) &&
               approxEqualNumbers(x[1], to!number(1.0)) &&
               approxEqualNumbers(x[2], to!number(1.0/3)),
               failedUnitTest());

        return 0;
    }
}
