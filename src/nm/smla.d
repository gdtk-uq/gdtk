/**
 * A minimal sparse matrix linear algebra module.
 *
 * Author: Rowan J. Gollan
 * Date: 2015-12-18
 *
 */

module nm.smla;

import core.stdc.stdlib : exit;
import std.stdio;
import std.string;
import std.conv;
import std.math;
import std.algorithm : reduce;

import nm.bbla;

/**
 * References used in this module:
 *
 * Saad (2003)
 * Iterative Methods for Sparse Linear Systems, 2nd ed.
 * SIAM, Philaldelphia
 */

double dot(double[] a, double[] b)
{
    assert(a.length == b.length);
    double sum = 0.0;
    foreach (i; 0 .. a.length) sum += a[i]*b[i];
    return sum;
}

double norm2(double[] vec)
{
    auto ssquares = reduce!((a,b) => a + b * b)(0.0, vec);
    return sqrt(ssquares);
}

/**
 * This class is used to store and work with a sparse matrix.
 * The data is stored in Compressed Sparse Row format.
 * Usually a class is used to hide details of implementation
 * (in this case, storage). However, here we will expose those
 * items as public method data since the algorithms to manipulate
 * sparse matrices rely on the representation of the data.
 *
 */
class SMatrix {
public:
    double[] aa;
    size_t[] ja;
    size_t[] ia;

    this() {}

    this(double[] aa, size_t[] ja, size_t[] ia)
    {
	this.aa = aa.dup;
	this.ja = ja.dup;
	this.ia = ia.dup;
    }

    this(SMatrix other)
    {
	this(other.aa, other.ja, other.ia);
    }

    void addRow(double[] ai, size_t[] ji) 
    {
	if ( ia.length == 0 )
	    ia ~= 0;
	aa ~= ai[];
	ja ~= ji[];
	ia ~= aa.length; // Now ia is already ready for next addition.
    }

    void scaleRow(size_t row, double scaleFactor)
    {
	assert(row < ia.length-1);
	foreach (j; ia[row] .. ia[row+1]) {
	    aa[j] *= scaleFactor;
	}
    }

    const double opIndex(size_t row, size_t col)
    {
	// We need to search the given row to see if an entry
	// is present in the given column.
	foreach (j; ia[row] .. ia[row+1]) {
	    if ( ja[j] == col )
		return aa[j];
	}
	// else search failed, so we have a 0.0 entry
	return 0.0;
    }
    
    // Presently we allow assignment if and only if the element
    // is already non-zero. Otherwise, we'd have to shuffle all
    // the elements around.
    ref double opIndexAssign(double c, size_t row, size_t col) {
	foreach ( j; ia[row] .. ia[row+1] ) {
	    if ( ja[j] == col ) {
		aa[j] = c;
		return aa[j];
	    }
	}
	// If we get here, we tried to assign to a zero value.
	throw new Error("ERROR: Tried to assign a value to a zero entry in sparse matrix.");
    }

    override string toString() {
	string s = "SMatrix[\n";
	foreach (row; 0 .. ia.length-1) {
	    foreach (col; 0 .. ia.length-1) {
		s ~= to!string(this[row,col]);
		if ( col < ia.length-2 )
		    s ~= ", ";
	    }
	    s ~= "\n";
	}
	s ~= "]";
	return s;
    }
}

bool approxEqualMatrix(SMatrix a, SMatrix b)
{
    if ( a.aa.length != b.aa.length ) return false;
    if ( a.ia.length != b.ia.length ) return false;
    // Test equality in terms of non-zero positions in matrix
    if ( a.ja != b.ja ) return false;
    if ( a.ia != b.ia ) return false;
    // Then test individual non-zero elements.
    foreach ( i; 0 .. a.aa.length ) {
	if ( !approxEqual(a.aa[i], b.aa[i]) ) return false;
    }
    return true;
}

void multiply(SMatrix a, double[] b, double[] c)
in {
    assert(a.ia.length-1 == b.length);
    assert(b.length == c.length);
}
body {
    size_t k0, k1;
    foreach ( i; 0 .. a.ia.length-1 ) {
	k0 = a.ia[i];
	k1 = a.ia[i+1];
	c[i] = 0.0;
	foreach ( k; k0 .. k1 ) c[i] += a.aa[k] * b[a.ja[k]];
    }
}

/**
 * An ILU(0) decomposition.
 *
 * This algorithm changes the matrix a in place leaving
 * an L/U matrix. The algorithm is an IKJ variant that
 * is useful for row contiguous data storage as is
 * used here.
 *
 * This implements Algorithm 10.4 in Saad (2003)
 */

void decompILU0(SMatrix a)
{
    size_t n = a.ia.length-1;
    foreach ( i; 1 .. n ) { // Begin from 2nd row
	foreach ( k; a.ja[a.ia[i] .. a.ia[i+1]] ) {
	    // Only work on columns k < i-1
	    if ( k >= i ) break;
	    a[i,k] = a[i,k]/a[k,k];
	    foreach ( j; a.ja[a.ia[i] .. a.ia[i+1]] ) {
		//Only work with columns j >= k+1
		if ( j <= k ) continue;
		a[i,j] = a[i,j] - a[i,k]*a[k,j];
	    }
	}
    }
}

/**
 * An ILU(p) decomposition.
 *
 * This algorithm changes the matrix a in place leaving
 * an L/U matrix. The algorithm is an IKJ variant that
 * is useful for row contiguous data storage as is
 * used here.
 *
 * This implements Algorithm 10.5 in Saad (2003)
 */

SMatrix decompILUp(SMatrix s, int p)
{
    size_t n = s.ia.length-1;
    // TODO: use sparse matrix
    Matrix a = new Matrix(n, n); // LU factorisation
    Matrix lev = new Matrix(n, n); // fill levels

    // store input sparse matrix in full matrix format
    foreach(i; 0 .. n) {
	foreach(j; 0 .. n) {
	    a[i,j] = s[i,j];
	}
    }
    
    // assign initial fill levels
    foreach ( i; 0 .. n ) { 
	foreach ( j; 0 .. n ) { 
	    if (a[i,j] == 0.0) lev[i,j] = n-1;
	    else lev[i,j] = 0;
	}
    }
    
    // factorise matrix
    foreach ( i; 1 .. n ) { // Begin from 2nd row
	foreach ( k; 0 .. i ) {
	    if (lev[i,k] <= p && a[k,k] != 0.0) {
		a[i,k] = a[i,k]/a[k,k];
		foreach (j; k+1 .. n) {
		    a[i,j] = a[i,j] - a[i,k]*a[k,j];
		    if (lev[i,j] > lev[i,k] + lev[k,j] + 1) {
			lev[i,j] = lev[i,k] + lev[k,j] +1;
		    } // end if
		} // end for
	    } // end if
	} // end for
	foreach(k; 0 .. n) {
	    if (lev[i,k] > p) a[i,k] = 0.0;
	} // end for
    }

    // convert back in sparse matrix
    auto M = new SMatrix();
    foreach(i; 0..n) {
	bool first_nonzero_val_in_row = false;
	foreach(j; 0..n) {
	    if (a[i,j] != 0.0) {
		M.aa ~= a[i,j]; 
		M.ja ~= j;
		if (first_nonzero_val_in_row == false) {
		    M.ia ~= M.aa.length-1;
		    first_nonzero_val_in_row = true;
		}
	    }
	}
	if (i == n-1) M.ia ~= M.aa.length;
    }

    return M;
}

void solve(SMatrix LU, double[] b)
{
    int n = to!int(LU.ia.length-1);
    assert(b.length == n);
    // Forward elimination
    foreach ( i; 1 .. n ) {
	foreach ( j; LU.ja[LU.ia[i] .. LU.ia[i+1]] ) {
	    // Only work up to i
	    if ( j >= i ) break;
	    double multiplier = LU[i,j];
	    b[i] -= multiplier * b[j];
	}
    }
    // Back substitution
    b[n-1] /= LU[n-1,n-1];
    for ( int i = to!int(n-2); i >= 0; --i ) {
	double sum = b[i];
	foreach ( j; LU.ja[LU.ia[i] .. LU.ia[i+1]] ) {
	    // Only work with j >= i+1
	    if ( j <= i ) continue;
	    sum -= LU[i,j] * b[j];
	}
	b[i] = sum/LU[i,i];
    }
}

double[] gmres(SMatrix A, double[] b, double[] x0, int m)
in {
    assert(A.ia.length-1 == b.length);
    assert(b.length == x0.length);
    assert(m >= 1);
}
body {
    immutable double ZERO_TOL = 1.0e-15;
    // 0. Allocate working matrices and vectors
    size_t n = b.length;
    double[] Ax0, r0, v, w, g_old, g, x;
    Ax0.length = n;
    r0.length = n;
    v.length = n;
    w.length = n;
    x.length = n;
    g_old.length = m+1;
    g.length = m+1;
    auto H = new Matrix(m+1, m);
    H.zeros();
    auto Hold = new Matrix(m+1, m);
    auto Gamma = new Matrix(m+1, m+1);
    auto V = new Matrix(n, m+1);

    // 1. Compute r0, beta, v1
    multiply(A, x0, Ax0);
    foreach (i; 0 .. n) {
	r0[i] = b[i] - Ax0[i];
    }

    auto beta = norm2(r0);
    foreach (i; 0 .. n) {
	v[i] = r0[i]/beta;
	V[i,0] = v[i];
    }

    // 2. Do 'm' iterations of update
    foreach (j; 0 .. m) {
	multiply(A, v, w);
	foreach (i; 0 .. j+1) {
	    v = V.getColumn(i);
	    H[i,j] = dot(w, v);
	    foreach (k; 0 .. n) {
		w[k] -= H[i,j]*v[k]; 
	    }
	}
	H[j+1,j] = norm2(w);
	if ( H[j+1,j] <= ZERO_TOL ) {
	    m = j + 1;
	    break;
	}
	foreach (i; 0 .. n) {
	    v[i] = w[i]/H[j+1,j];
	    V[i,j+1] = v[i];
	}
    }

    // Use the plane rotations method to compute solution
    // and residual at each step.
    g_old[] = 0.0;
    g_old[0] = beta;
    copy(H, Hold);
    foreach (i; 0 .. m) {
	double c_i, s_i, denom;
	denom = sqrt(H[i,i]*H[i,i] + H[i+1,i]*H[i+1,i]);
	s_i = H[i+1,i]/denom; 
	c_i = H[i,i]/denom;
	Gamma.eye();
	Gamma[i,i] = c_i; Gamma[i,i+1] = s_i;
	Gamma[i+1,i] = -s_i; Gamma[i+1,i+1] = c_i;
	nm.bbla.dot(Gamma, Hold, H);
	nm.bbla.dot(Gamma, g_old, g);
	// Get Qold and g_old ready for next step
	g_old[] = g[];
	copy(H, Hold);
    }

    auto resid = fabs(g[m]);
    auto R = H.sliceDup(0, m, 0, m);
    auto gm = g[0 .. m];
    // At end H := R
    //        g := gm
    upperSolve(R, gm);
    auto Vm = V.sliceDup(0, n, 0, m);
    nm.bbla.dot(Vm, gm, x);

    foreach (i; 0 .. n) x[i] += x0[i];

    multiply(A, x, Ax0);
    foreach (i; 0 .. n) {
	r0[i] = b[i] - Ax0[i];
    }

    return x.dup();
}

double[] gmres2(SMatrix A, double[] b, double[] x0, int maxIters, double residTol)
in {
    assert(A.ia.length-1 == b.length);
    assert(b.length == x0.length);
    assert(maxIters >= 1);
}
body {
    double resid;
    // 0. Allocate working matrices and vectors
    size_t n = b.length;
    size_t m = maxIters;
    double[] Ax0, r0, v, w, x;
    Ax0.length = n;
    r0.length = n;
    v.length = n;
    w.length = n;
    x.length = n;
    auto V = new Matrix(n, m+1);
    auto H0 = new Matrix(m+1, m); H0.zeros();
    auto H1 = new Matrix(m+1, m); H1.zeros();
    auto Gamma = new Matrix(m+1, m+1); Gamma.eye();
    auto Q0 = new Matrix(m+1, m+1);
    auto Q1 = new Matrix(m+1, m+1);
    double[] g0, g1;
    g0.length = m+1; g0[] = 0.0;
    g1.length = m+1; g1[] = 0.0;
    double[] h, hR;
    h.length = m+1;
    hR.length = m+1;

    // 1. Compute r0, beta, v1
    multiply(A, x0, Ax0);
    foreach (i; 0 .. n) {
	r0[i] = b[i] - Ax0[i];
    }

    auto beta = norm2(r0);
    g0[0] = beta;
    foreach (i; 0 .. n) {
	v[i] = r0[i]/beta;
	V[i,0] = v[i];
    }

    // 2. Do 'm' iterations of update
    foreach (j; 0 .. m) {
	multiply(A, v, w);
	foreach (i; 0 .. j+1) {
	    v = V.getColumn(i);
	    H0[i,j] = dot(w, v);
	    foreach (k; 0 .. n) {
		w[k] -= H0[i,j]*v[k]; 
	    }
	}
	H0[j+1,j] = norm2(w);
	
	foreach (i; 0 .. n) {
	    v[i] = w[i]/H0[j+1,j];
	    V[i,j+1] = v[i];
	}

	// Build rotated Hessenberg progressively
	if ( j != 0 ) {
	    // Extract final column in H
	    foreach (i; 0 .. j+1) h[i] = H0[i,j];
	    // Rotate column by previous rotations (stored in Q0)
	    nm.bbla.dot(Q0, j+1, j+1, h, hR);
	    // Place column back in H
	    foreach (i; 0 .. j+1) H0[i,j] = hR[i];
	}
	// Now form new Gamma
	Gamma.eye();
	double c_j, s_j, denom;
	denom = sqrt(H0[j,j]*H0[j,j] + H0[j+1,j]*H0[j+1,j]);
	s_j = H0[j+1,j]/denom; 
	c_j = H0[j,j]/denom;
	Gamma[j,j] = c_j; Gamma[j,j+1] = s_j;
	Gamma[j+1,j] = -s_j; Gamma[j+1,j+1] = c_j;
	// Apply rotations
	nm.bbla.dot(Gamma, j+2, j+2, H0, j+1, H1);
	nm.bbla.dot(Gamma, j+2, j+2, g0, g1);
	// Accumulate Gamma rotations in Q.
	if ( j == 0 ) {
	    copy(Gamma, Q1);
	}
	else {
	    nm.bbla.dot(Gamma, j+2, j+2, Q0, j+2, Q1);
	}
	// Get residual
	resid = fabs(g1[j+1]);
	if ( resid <= residTol ) {
	    m = j+1;
	    break;
	}

	// Prepare for next step
	copy(H1, H0);
	g0[] = g1[];
	copy(Q1, Q0);
    }

    auto R = H1.sliceDup(0, m, 0, m);
    auto gm = g1[0 .. m];
    // At end H := R
    //        g := gm
    upperSolve(R, gm);
    auto Vm = V.sliceDup(0, n, 0, m);
    nm.bbla.dot(Vm, gm, x);

    foreach (i; 0 .. n) x[i] += x0[i];

    multiply(A, x, Ax0);
    foreach (i; 0 .. n) {
	r0[i] = b[i] - Ax0[i];
    }

    return x.dup();
}

struct GMRESWorkSpace {
    size_t n, m;
    double[] Ax0, r0, v, w, Pv, g0, g1, h, hR;
    Matrix V, H0, H1, Gamma, Q0, Q1;

    this(size_t n, int maxIters) {
	this.n = n;
	this.m = maxIters;
	// Allocate arrays
	this.Ax0.length = n;
	this.r0.length = n;
	this.v.length = n;
	this.w.length = n;
	this.Pv.length = n;
	this.g0.length = m+1; 
	this.g1.length = m+1; 
	this.h.length = m+1;
	this.hR.length = m+1;
	// Allocate matrices
	this.V = new Matrix(n, m+1);
	this.H0 = new Matrix(m+1, m); 
	this.H1 = new Matrix(m+1, m); 
	this.Gamma = new Matrix(m+1, m+1); 
	this.Q0 = new Matrix(m+1, m+1);
	this.Q1 = new Matrix(m+1, m+1);
    }
}

/**
 * A GMRES iterative solver with right preconditioning.
 * 
 * Params: 
 *    A =         coefficient matrix
 *    P =         pre-conditioning matrix in LU form
 *    b =         RHS vector
 *    x0 =        initial guess for solution
 *    x =         solution vector place in x
 *    maxIters =  maximum number of iterations
 *    residTol =  stop iterations when residual below this tolerance
 *    gws =       pre-allocated workspace
 */
void rpcGMRES(SMatrix A, SMatrix P, double[] b, double[] x0, double[] x,
	      int maxIters, double residTol, ref GMRESWorkSpace gws)
in {
    assert(A.ia.length-1 == b.length);
    assert(b.length == x0.length);
    assert(maxIters >= 1);
}
body {
    double resid;
    size_t n = b.length;
    size_t m = maxIters;
    // 0. Initialise the values in some of the pre-allocated storage
    gws.g0[] = 0.0;
    gws.g1[] = 0.0;
    gws.H0.zeros();
    gws.H1.zeros();
    gws.Gamma.eye();

    // 1. Compute r0, beta, v1
    multiply(A, x0, gws.Ax0);
    foreach (i; 0 .. n) {
	gws.r0[i] = b[i] - gws.Ax0[i];
    }

    auto beta = norm2(gws.r0);
    gws.g0[0] = beta;
    foreach (i; 0 .. n) {
	gws.v[i] = gws.r0[i]/beta;
	gws.V[i,0] = gws.v[i];
    }


    // 2. Do 'm' iterations of update
    foreach (j; 0 .. m) {
	gws.Pv[] = gws.v[];
	solve(P, gws.Pv);
	multiply(A, gws.Pv, gws.w);
	foreach (i; 0 .. j+1) {
	    foreach (k; 0 .. n ) gws.v[k] = gws.V[k,i]; // Extract column 'i'
	    gws.H0[i,j] = dot(gws.w, gws.v);
	    foreach (k; 0 .. n) gws.w[k] -= gws.H0[i,j]*gws.v[k]; 
	}
	gws.H0[j+1,j] = norm2(gws.w);
	
	foreach (k; 0 .. n) {
	    gws.v[k] = gws.w[k]/gws.H0[j+1,j];
	    gws.V[k,j+1] = gws.v[k];
	}

	// Build rotated Hessenberg progressively
	if ( j != 0 ) {
	    // Extract final column in H
	    foreach (i; 0 .. j+1) gws.h[i] = gws.H0[i,j];
	    // Rotate column by previous rotations (stored in Q0)
	    nm.bbla.dot(gws.Q0, j+1, j+1, gws.h, gws.hR);
	    // Place column back in H
	    foreach (i; 0 .. j+1) gws.H0[i,j] = gws.hR[i];
	}
	// Now form new Gamma
	gws.Gamma.eye();
	auto denom = sqrt(gws.H0[j,j]*gws.H0[j,j] + gws.H0[j+1,j]*gws.H0[j+1,j]);
	auto s_j = gws.H0[j+1,j]/denom; 
	auto c_j = gws.H0[j,j]/denom;
	gws.Gamma[j,j] = c_j; gws.Gamma[j,j+1] = s_j;
	gws.Gamma[j+1,j] = -s_j; gws.Gamma[j+1,j+1] = c_j;
	// Apply rotations
	nm.bbla.dot(gws.Gamma, j+2, j+2, gws.H0, j+1, gws.H1);
	nm.bbla.dot(gws.Gamma, j+2, j+2, gws.g0, gws.g1);
	// Accumulate Gamma rotations in Q.
	if ( j == 0 ) {
	    copy(gws.Gamma, gws.Q1);
	}
	else {
	    nm.bbla.dot(gws.Gamma, j+2, j+2, gws.Q0, j+2, gws.Q1);
	}
	// Get residual
	resid = fabs(gws.g1[j+1]);
	if ( resid <= residTol ) {
	    m = j+1;
	    break;
	}

	// Prepare for next step
	copy(gws.H1, gws.H0);
	gws.g0[] = gws.g1[];
	copy(gws.Q1, gws.Q0);
    }

    // At end H := R up to row m
    //        g := gm up to row m
    upperSolve(gws.H1, to!int(m), gws.g1);
    nm.bbla.dot(gws.V, n, m, gws.g1, x);
    solve(P, x);

    foreach (i; 0 .. n) x[i] += x0[i];
}

double[] fGMRES(SMatrix A, SMatrix P, double[] b, double[] x0, int mOuter, int mInner)
in {
    assert(A.ia.length-1 == b.length);
    assert(b.length == x0.length);
    assert(mOuter >= 1);
}
body {
    immutable double ZERO_TOL = 1.0e-15;
    // 0. Allocate working matrices and vectors
    size_t n = b.length;
    double[] Ax0, r0, v, w, g_old, g, x, z, x0_inner;
    Ax0.length = n;
    r0.length = n;
    v.length = n;
    z.length = n;
    w.length = n;
    x.length = n;
    x0_inner.length = n;
    x0_inner[] = 0.0;
    g_old.length = mOuter+1;
    g.length = mOuter+1;
    auto H = new Matrix(mOuter+1, mOuter);
    H.zeros();
    auto Hold = new Matrix(mOuter+1, mOuter);
    auto Gamma = new Matrix(mOuter+1, mOuter+1);
    auto V = new Matrix(n, mOuter+1);
    auto Z = new Matrix(n, mOuter+1);

    // 1. Compute r0, beta, v1
    multiply(A, x0, Ax0);
    foreach (i; 0 .. n) {
	r0[i] = b[i] - Ax0[i];
    }

    auto beta = norm2(r0);
    foreach (i; 0 .. n) {
	v[i] = r0[i]/beta;
	V[i,0] = v[i];
    }

    // 2. Do 'mOuter' iterations of update
    foreach (j; 0 .. mOuter) {
	z = gmres(P, v, x0_inner, mInner);
	// Save z vector in Z
	foreach (k; 0 .. n ) Z[k,j] = z[k];
	multiply(A, z, w);
	foreach (i; 0 .. j+1) {
	    v = V.getColumn(i);
	    H[i,j] = dot(w, v);
	    foreach (k; 0 .. n) {
		w[k] -= H[i,j]*v[k]; 
	    }
	}
	H[j+1,j] = norm2(w);
	if ( H[j+1,j] <= ZERO_TOL ) {
	    mOuter = j;
	    break;
	}
	foreach (i; 0 .. n) {
	    v[i] = w[i]/H[j+1,j];
	    V[i,j+1] = v[i];
	}
    }

    // Use the plane rotations method to compute solution
    // and residual at each step.
    g_old[] = 0.0;
    g_old[0] = beta;
    copy(H, Hold);
    foreach (i; 0 .. mOuter) {
	double c_i, s_i, denom;
	denom = sqrt(H[i,i]*H[i,i] + H[i+1,i]*H[i+1,i]);
	s_i = H[i+1,i]/denom; 
	c_i = H[i,i]/denom;
	Gamma.eye();
	Gamma[i,i] = c_i; Gamma[i,i+1] = s_i;
	Gamma[i+1,i] = -s_i; Gamma[i+1,i+1] = c_i;
	nm.bbla.dot(Gamma, Hold, H);
	nm.bbla.dot(Gamma, g_old, g);
	// Get Qold and g_old ready for next step
	g_old[] = g[];
	copy(H, Hold);
    }

    auto resid = fabs(g[mOuter]);
    auto R = H.sliceDup(0, mOuter, 0, mOuter);
    auto gm = g[0 .. mOuter];
    // At end H := R
    //        g := gm
    upperSolve(R, gm);
    auto Zm = Z.sliceDup(0, n, 0, mOuter);
    nm.bbla.dot(Zm, gm, x);

    foreach (i; 0 .. n) x[i] += x0[i];

    return x.dup();
}



version(smla_test) {
    import util.msg_service;
    int main() {
	// This example matrix is from Saad (2003), Sec. 3.4
	auto a = new SMatrix([1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12.],
			     [0, 3, 0, 1, 3, 0, 2, 3, 4, 2, 3, 4],
			     [0, 2, 5, 9, 11, 12]);
	// Test construction of matrix row by row
	auto b = new SMatrix();
	b.addRow([1., 2.], [0, 3]);
	b.addRow([3., 4., 5.], [0, 1, 3]);
	b.addRow([6., 7., 8., 9.], [0, 2, 3, 4]);
	b.addRow([10., 11.], [2, 3]);
	b.addRow([12.], [4]);
	assert(approxEqualMatrix(a, b), failedUnitTest());
	// Test matrix multiply
	double[] v = [1., 2., 3., 4., 5.];
	double[] c;
	c.length = v.length;
	multiply(a, v, c);
	double[] expected_c = [9., 31., 104., 74., 60.];
	foreach ( i; 0 .. c.length ) {
	    assert(approxEqual(c[i], expected_c[i]), failedUnitTest());
	}
	// Test decompILU0
	// This example matrix and decomposition is from Gerard and Wheatley, 6th ed, Sec. 2.4
	auto d = new SMatrix();
	d.addRow([-4., 2.], [0, 1]);
	d.addRow([1., -4., 1.], [0, 1, 2]);
	d.addRow([1., -4., 1.], [1, 2, 3]);
	d.addRow([1., -4., 1.], [2, 3, 4]);
	d.addRow([2., -4.], [3, 4]);
	decompILU0(d);
	auto dLU = new SMatrix();
	dLU.addRow([-4., 2.], [0, 1]);
	dLU.addRow([-0.25, -3.5, 1.], [0, 1, 2]);
	dLU.addRow([-0.2857, -3.7143, 1.], [1, 2, 3]);
	dLU.addRow([-0.2692, -3.7308, 1.], [2, 3, 4]);
	dLU.addRow([-0.5361, -3.4639], [3, 4]);
	assert(approxEqualMatrix(d, dLU), failedUnitTest());
	// Test solve.
	// Let's give the ILU(0) method a triangular matrix that is can solve exactly.
	// This example is taken from Faires and Burden (1993), Sec. 6.6, example 3. 
	auto e = new SMatrix();
	e.addRow([2., -1.], [0, 1]);
	e.addRow([-1., 2., -1.], [0, 1, 2]);
	e.addRow([-1., 2., -1.], [1, 2, 3]);
	e.addRow([-1., 2.], [2, 3]);
	decompILU0(e);
	double[] B = [1., 0., 0., 1.];
	solve(e, B);
	double[] B_exp = [1., 1., 1., 1.];
	foreach ( i; 0 .. B.length ) {
	    assert(approxEqual(B[i], B_exp[i]), failedUnitTest());
	}
	// Now let's see how we go at an approximate solve by using a non-triangular matrix.
	// This is example 2.2 from Gerard and Wheatley, 6th edition
	auto f = new SMatrix();
	f.addRow([3., 2., -1., 2.], [0, 1, 2, 3]);
	f.addRow([1., 4., 2.], [0, 1, 3]);
	f.addRow([2., 1., 2., -1.], [0, 1, 2, 3]);
	f.addRow([1., 1., -1., 3.], [0, 1, 2, 3]);
	decompILU0(f);
	double[] C = [2., 2., 0., 0.];
	solve(f, C);
	double[] C_exp = [0.333333, 0.666667, -1, -0.666667];
	foreach ( i; 0 .. C.length ) {
	    assert(approxEqual(C[i], C_exp[i]), failedUnitTest());
	}
	
	// Let's test the ILU(p) method
	SMatrix s = new SMatrix([1., 1., 4., 2., 4., 1., 2., 1., 8., 2., 4., 1., 3., 6., 2., 1.],
				[0, 1, 4, 1, 2, 4, 0, 1, 2, 3, 2, 3, 0, 1, 2, 4],
				[0, 3, 6, 10, 12, 16]);
	int p;
	SMatrix m = new SMatrix();
	// test for ILU(p=0)
	p = 0;
	m = decompILUp(s, p);
	auto sol0 = new SMatrix([1., 1., 4., 2., 4., 1., 2., -0.5, 10., 2., 0.4, 0.2, 3., 1.5, -0.4, -12.5],
				[0, 1, 4, 1, 2, 4, 0, 1, 2, 3, 2, 3, 0, 1, 2, 4],
				[0, 3, 6, 10, 12, 16]);
	assert(approxEqualMatrix(m, sol0), failedUnitTest());
	// test for ILU(p=2)
	p = 2;
	m = decompILUp(s, p);
	auto sol2 = new SMatrix([1., 1., 4., 2., 4., 1., 2., -0.5, 10., 2., -7.5, 0.4, 0.2, 3., 3., 1.5, -0.4, 4., -27.5],
				[0, 1, 4, 1, 2, 4, 0, 1, 2, 3, 4, 2, 3, 4, 0, 1, 2, 3, 4],
				[0, 3, 6, 11, 14, 19]);
	assert(approxEqualMatrix(m, sol2), failedUnitTest());
	
	// Test GMRES on Faires and Burden problem.

	auto g = new SMatrix();
	g.addRow([2., -1.], [0, 1]);
	g.addRow([-1., 2., -1.], [0, 1, 2]);
	g.addRow([-1., 2., -1.], [1, 2, 3]);
	g.addRow([-1., 2.], [2, 3]);
	double[] B1 = [1., 0., 0., 1.];
	double[] x0 = [1.2, 0.8, 0.9, 1.1];

	
	auto x = gmres(g, B1, x0, 4);
	foreach (i; 0 .. x.length) {
	    assert(approxEqual(x[i], B_exp[i]), failedUnitTest());
	}
	
	x = gmres2(g, B1, x0, 5, 1.0e-10);
	foreach (i; 0 .. x.length) {
	    assert(approxEqual(x[i], B_exp[i]), failedUnitTest());
	}

	// Test pre-conditioned GMRES on Gerald and Wheatley problem.
	// This time we expect the exact answer. Earlier we only used
	// an incomplete LU factorisation and so the result was
	// only approximate.
	auto h = new SMatrix();
	h.addRow([3., 2., -1., 2.], [0, 1, 2, 3]);
	h.addRow([1., 4., 2.], [0, 1, 3]);
	h.addRow([2., 1., 2., -1.], [0, 1, 2, 3]);
	h.addRow([1., 1., -1., 3.], [0, 1, 2, 3]);
	auto Ph = new SMatrix(h);
	decompILU0(Ph);
	double[] C1 = [2., 2., 0., 0.];
	x0 = [0.2, 0.5, -1.1, -0.6];
	int maxIters = 5;
	auto gws = GMRESWorkSpace(x0.length, maxIters);
	rpcGMRES(h, Ph, C1, x0, x, maxIters, 1.0e-15, gws);
	double[] C1_exp = [0.273, 0.773, -1.0, -0.682];

	foreach (i; 0 .. x.length) {
	    assert(approxEqual(x[i], C1_exp[i]), failedUnitTest());
	}

	auto Pi = new SMatrix();
	Pi.addRow([1.,], [0]);
	Pi.addRow([1.,], [1]);
	Pi.addRow([1.,], [2]);
	Pi.addRow([1.,], [3]);

	x = fGMRES(h, Pi, C1, x0, 5, 3);
	foreach (i; 0 .. x.length) {
	    assert(approxEqual(x[i], C1_exp[i]), failedUnitTest());
	}

	return 0;
    }
}
