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

@nogc
int computeInverse(int N, int NDIM)
    (ref double[2*NDIM][NDIM] c, double very_small_value=1.0e-16)
// Perform Gauss-Jordan elimination on an augmented matrix.
// c = [A|b] such that the mutated matrix becomes [I|x]
// where x is the solution vector(s) to A.x = b
//
// Returns 0 normally, -1 if the matrix is essentially singular (zero pivot).
//
// We allow a larger storage than needed because the flow solver code will
// use these functions for both 2D and 3D simulations.
// When N < NDIM, the solver just uses a subset of the available storage
// and the rest remains untouched.
{
    assert(NDIM >= N, "Inadequate size of dimension for matrix");
    foreach(j; 0 .. N) {
	// Select pivot.
	size_t p = j;
	foreach(i; j+1 .. N) {
	    if ( abs(c[i][j]) > abs(c[p][j]) ) p = i;
	}
	if (abs(c[p][j]) <= very_small_value) return -1; // singular
	if ( p != j ) { // Swap rows
	    foreach(col; 0 .. 2*N) {
		double tmp = c[p][col]; c[p][col] = c[j][col]; c[j][col] = tmp;
	    }
	}
	// Scale row j to get unity on the diagonal.
	double cjj = c[j][j];
	foreach(col; 0 .. 2*N) c[j][col] /= cjj;
	// Do the elimination to get zeros in all off diagonal values in column j.
	foreach(i; 0 .. N) {
	    if ( i == j ) continue;
	    double cij = c[i][j];
	    foreach(col; 0 .. 2*N) c[i][col] -= cij * c[j][col]; 
	}
    } // end foreach j
    return 0; // success
} // end computeInverse()()

int computeInverseDebug(int N, int NDIM)
    (ref double[2*NDIM][NDIM] c, double very_small_value=1.0e-16)
// Perform Gauss-Jordan elimination on an augmented matrix.
// c = [A|b] such that the mutated matrix becomes [I|x]
// where x is the solution vector(s) to A.x = b
//
// Returns 0 normally, -1 if the matrix is essentially singular (zero pivot).
//
// This debug version gives more information.
{
    assert(NDIM >= N, "Inadequate size of dimension for matrix");
    double[2*N][N] csave;
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
		double tmp = c[p][col]; c[p][col] = c[j][col]; c[j][col] = tmp;
	    }
	}
	// Scale row j to get unity on the diagonal.
	double cjj = c[j][j];
	foreach(col; 0 .. 2*N) c[j][col] /= cjj;
	// Do the elimination to get zeros in all off diagonal values in column j.
	foreach(i; 0 .. N) {
	    if ( i == j ) continue;
	    double cij = c[i][j];
	    foreach(col; 0 .. 2*N) c[i][col] -= cij * c[j][col]; 
	}
    } // end foreach j
    return 0; // success
} // end computeInverseDebug()()

@nogc
void solveWithInverse(int N, int NDIM)
    (ref double[2*NDIM][NDIM] c, ref double[NDIM] rhs, ref double[NDIM] x)
// Multiply right-hand-side by the inverse part of the augmented matrix.
{
    assert(NDIM >= N, "Inadequate size of dimension for matrix");
    foreach(i; 0 .. N) {
	x[i] = 0.0;
	foreach(j; 0 .. N) {
	    x[i] += c[i][N+j] * rhs[j];
	}
    }
} // end solveWithInverse()()


version(rsla_test) {
    import util.msg_service;
    int main() {
	double[8][4] A = [[0.0,  2.0,  0.0,  1.0,  1.0, 0.0, 0.0, 0.0],
			  [2.0,  2.0,  3.0,  2.0,  0.0, 1.0, 0.0, 0.0],
			  [4.0, -3.0,  0.0,  1.0,  0.0, 0.0, 1.0, 0.0],
			  [6.0,  1.0, -6.0, -5.0,  0.0, 0.0, 0.0, 1.0]];
	computeInverse!(4,4)(A);
	double[4] b = [0.0, -2.0, -7.0, 6.0];
	double[4] x;
	solveWithInverse!(4,4)(A, b, x);
	assert(approxEqual(x[0], -0.5) && approxEqual(x[1], 1.0) &&
	       approxEqual(x[2], 1.0/3) && approxEqual(x[3], -2.0), failedUnitTest());

	return 0;
    }
}
