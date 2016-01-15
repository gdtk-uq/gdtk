/**
 * rsla: real simple linear algebra
 * A linear equation solver using the text-book notation.
 *
 * PJ 2016-01-16 factored out of the gradient calculation module for eilmer
 */

module nm.rsla;

import std.math;

@nogc
void computeInverse(int N)(ref double[2*N][N] c, double very_small_value=1.0e-16)
// Perform Gauss-Jordan elimination on an augmented matrix.
// c = [A|b] such that the mutated matrix becomes [I|x]
// where x is the solution vector(s) to A.x = b
{
    foreach(j; 0 .. N) {
	// Select pivot.
	size_t p = j;
	foreach(i; j+1 .. N) {
	    if ( abs(c[i][j]) > abs(c[p][j]) ) p = i;
	}
	assert(abs(c[p][j]) > very_small_value, "matrix is essentially singular");
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
} // end computeInverse()()

@nogc
void solveWithInverse(int N)(ref double[2*N][N] c, ref double[N] rhs, ref double[N] x)
// Multiply right-hand-side by the inverse part of the augmented matrix.
{
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
	computeInverse!4(A);
	double[4] b = [0.0, -2.0, -7.0, 6.0];
	double[4] x;
	solveWithInverse!4(A, b, x);
	assert(approxEqual(x[0], -0.5) && approxEqual(x[1], 1.0) &&
	       approxEqual(x[2], 1.0/3) && approxEqual(x[3], -2.0), failedUnitTest());

	return 0;
    }
}
