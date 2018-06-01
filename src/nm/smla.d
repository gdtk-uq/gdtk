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
import std.array;
import core.memory;
import std.algorithm : reduce;
import nm.bbla;
import nm.complex;
import nm.number;

immutable double ESSENTIALLY_ZERO = 1.0e-50;

/**
 * References used in this module:
 *
 * Saad (2003)
 * Iterative Methods for Sparse Linear Systems, 2nd ed.
 * SIAM, Philaldelphia
 */

T dot(T)(T[] a, T[] b)
{
    assert(a.length == b.length);
    T sum = 0.0;
    foreach (i; 0 .. a.length) sum += a[i]*b[i];
    return sum;
}

T norm2(T)(T[] vec)
{
    auto ssquares = reduce!((a,b) => a + b * b)(to!T(0.0), vec);
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
class SMatrix(T) {
public:
    T[] aa;
    size_t[] ja;
    size_t[] ia;

    this() {}

    this(T[] aa, size_t[] ja, size_t[] ia)
    {
        this.aa = aa.dup;
        this.ja = ja.dup;
        this.ia = ia.dup;
    }

    this(SMatrix other)
    {
        this(other.aa, other.ja, other.ia);
    }

    void addRow(T[] ai, size_t[] ji) 
    {
        if ( ia.length == 0 )
            ia ~= 0;
        aa ~= ai[];
        ja ~= ji[];
        ia ~= aa.length; // Now ia is already ready for next addition.
    }

    void scaleRow(size_t row, T scaleFactor)
    {
        assert(row < ia.length-1);
        foreach (j; ia[row] .. ia[row+1]) {
            aa[j] *= scaleFactor;
        }
    }

    const T opIndex(size_t row, size_t col)
    {
        // We need to search the given row to see if an entry
        // is present in the given column.
        foreach (j; ia[row] .. ia[row+1]) {
            if ( ja[j] == col )
                return aa[j];
        }
        // else search failed, so we have a 0.0 entry
        return to!T(0.0);
    }
    
    // We now allow assingment to zero entries.
    ref T opIndexAssign(T c, size_t row, size_t col) {
        foreach ( j; ia[row] .. ia[row+1] ) {
            if ( ja[j] == col ) {
                aa[j] = c;
                return aa[j];
            }
        }

        // If we get here, we tried to assign to a zero value,
        // so to add the entry we will need to shuffle all the elements
        // WARNING: this code will happily add a zero to the sparse matrix;
        //          onus is on any user using this object to prevent
        //          unnecessary assingment of 0 values.
        
        // shuffle aa, ja, ia entries
        if ( col > ja[ia[row+1]-1] ) { // if the new entry occurs after all currently stored columns in the row
            aa.insertInPlace(ia[row+1], c);
            ja.insertInPlace(ia[row+1], col);
            foreach ( k; row+1 .. ia.length ) this.ia[k] += 1;
            return aa[ia[row+1]];
        }
        foreach ( j; ia[row] .. ia[row+1] ) { // else handle all other scenarios in a generic manner
            if ( col < ja[j] ) {
                aa.insertInPlace(j, c);
                ja.insertInPlace(j, col);
                foreach ( k; row+1 .. ia.length ) this.ia[k] += 1;
                return aa[j];
            }
        }
        // if we get here, something has gone wrong
        throw new Error("ERROR: An error occured while assigning a value to an entry in the sparse matrix.");
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
} // end class SMatrix

bool approxEqualMatrix(T)(SMatrix!T a, SMatrix!T b)
{
    if ( a.aa.length != b.aa.length ) return false;
    if ( a.ia.length != b.ia.length ) return false;
    // Test equality in terms of non-zero positions in matrix
    if ( a.ja != b.ja ) return false;
    if ( a.ia != b.ia ) return false;
    // Then test individual non-zero elements.
    foreach ( i; 0 .. a.aa.length ) {
        if ( !approxEqualNumbers(a.aa[i], b.aa[i]) ) return false;
    }
    return true;
}

void multiply(T)(SMatrix!T a, T[] b, T[] c)
in {
    assert(a.ia.length-1 == c.length);
    // Some faith is given that the user provides an appropriate length
    // input vector b, since in CRS we cannot easily determine the number of columns.
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

void decompILU0(T)(SMatrix!T a)
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

void decompILUp(T)(SMatrix!T a, int k)
{
    // NB. This pre-conditioner uses a sparse matrix, however some stored entries will be zero,
    // this is a result of entries starting as non-zero in the original matrix, and then becoming
    // zero during the factorisation.
    
    int n = to!int(a.ia.length)-1;
    int[][] lev; // fill levels
    lev.length = n;
    foreach ( i; 0..n) lev[i].length = n;

    // assign initial fill levels
    foreach ( i; 0 .. n ) { 
	foreach ( j; 0 .. n ) { 
            if (abs(a[i,j]) <= ESSENTIALLY_ZERO) lev[i][j] = n-1;
	    else lev[i][j] = 0;
	}
    }

    // symbolic phase
    foreach ( i; 1 .. n ) { // Begin from 2nd row
        foreach ( p; 0 .. i ) {
            if (lev[i][p] <= k) {
                foreach ( j ; p..n) {
                    if (lev[i][j] > ESSENTIALLY_ZERO && lev[i][j] > (lev[i][p]+lev[p][j]+1) ) {
                        lev[i][j] = lev[i][p]+lev[p][j]+1;
                    }
                } // end foreach
            } // end if
        } // end foreach
    } // end foreach
    
    // modify a matrix nonzero pattern
    foreach ( i; 0..n) {
        foreach ( j; 0..n) {
            if (lev[i][j] <= k && abs(a[i,j]) < ESSENTIALLY_ZERO) { a[i,j] = to!T(0.0); } 
        }
    }
    
    // clear the fill level matrix from memory
    destroy(lev);
    GC.minimize();
    
    // factorise phase (using ILU0 algorithm on new sparsity pattern)
    decompILU0!T(a);
}



void solve(T)(SMatrix!T LU, T[] b)
{
    int n = to!int(LU.ia.length-1);
    assert(b.length == n);
    // Forward elimination
    foreach ( i; 1 .. n ) {
        foreach ( j; LU.ja[LU.ia[i] .. LU.ia[i+1]] ) {
            // Only work up to i
            if ( j >= i ) break;
            T multiplier = LU[i,j];
            b[i] -= multiplier * b[j];
        }
    }
    // Back substitution
    b[n-1] /= LU[n-1,n-1];
    for ( int i = to!int(n-2); i >= 0; --i ) {
        T sum = b[i];
        foreach ( j; LU.ja[LU.ia[i] .. LU.ia[i+1]] ) {
            // Only work with j >= i+1
            if ( j <= i ) continue;
            sum -= LU[i,j] * b[j];
        }
        b[i] = sum/LU[i,i];
    }
}

T[] gmres(T)(SMatrix!T A, T[] b, T[] x0, int m)
in {
    assert(A.ia.length-1 == b.length);
    assert(b.length == x0.length);
    assert(m >= 1);
}
body {
    immutable double ZERO_TOL = 1.0e-15;
    // 0. Allocate working matrices and vectors
    size_t n = b.length;
    T[] Ax0, r0, v, w, g_old, g, x;
    Ax0.length = n;
    r0.length = n;
    v.length = n;
    w.length = n;
    x.length = n;
    g_old.length = m+1;
    g.length = m+1;
    auto H = new Matrix!T(m+1, m);
    H.zeros();
    auto Hold = new Matrix!T(m+1, m);
    auto Gamma = new Matrix!T(m+1, m+1);
    auto V = new Matrix!T(n, m+1);

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
    foreach( idx; 0..g_old.length) g_old[idx] = 0.0;
    g_old[0] = beta;
    copy(H, Hold);
    foreach (i; 0 .. m) {
        T c_i, s_i, denom;
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
    upperSolve!T(R, gm);
    auto Vm = V.sliceDup(0, n, 0, m);
    nm.bbla.dot!T(Vm, gm, x);

    foreach (i; 0 .. n) x[i] += x0[i];

    multiply!T(A, x, Ax0);
    foreach (i; 0 .. n) {
        r0[i] = b[i] - Ax0[i];
    }

    return x.dup();
}

T[] gmres2(T)(SMatrix!T A, T[] b, T[] x0, int maxIters, double residTol)
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
    T[] Ax0, r0, v, w, x;
    Ax0.length = n;
    r0.length = n;
    v.length = n;
    w.length = n;
    x.length = n;
    auto V = new Matrix!T(n, m+1);
    auto H0 = new Matrix!T(m+1, m); H0.zeros();
    auto H1 = new Matrix!T(m+1, m); H1.zeros();
    auto Gamma = new Matrix!T(m+1, m+1); Gamma.eye();
    auto Q0 = new Matrix!T(m+1, m+1);
    auto Q1 = new Matrix!T(m+1, m+1);
    T[] g0, g1;
    g0.length = m+1; foreach(idx; 0..g0.length) { g0[idx] = 0.0; }
    g1.length = m+1; foreach(idx; 0..g1.length) { g1[idx] = 0.0; }
    T[] h, hR;
    h.length = m+1;
    hR.length = m+1;

    // 1. Compute r0, beta, v1
    multiply!T(A, x0, Ax0);
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
            H0[i,j] = dot!T(w, v);
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
            nm.bbla.dot!T(Q0, j+1, j+1, h, hR);
            // Place column back in H
            foreach (i; 0 .. j+1) H0[i,j] = hR[i];
        }
        // Now form new Gamma
        Gamma.eye();
        T c_j, s_j, denom;
        denom = sqrt(H0[j,j]*H0[j,j] + H0[j+1,j]*H0[j+1,j]);
        s_j = H0[j+1,j]/denom; 
        c_j = H0[j,j]/denom;
        Gamma[j,j] = c_j; Gamma[j,j+1] = s_j;
        Gamma[j+1,j] = -s_j; Gamma[j+1,j+1] = c_j;
        // Apply rotations
        nm.bbla.dot!T(Gamma, j+2, j+2, H0, j+1, H1);
        nm.bbla.dot!T(Gamma, j+2, j+2, g0, g1);
        // Accumulate Gamma rotations in Q.
        if ( j == 0 ) {
            copy(Gamma, Q1);
        }
        else {
            nm.bbla.dot!T(Gamma, j+2, j+2, Q0, j+2, Q1);
        }
        // Get residual
        resid = fabs(g1[j+1]).re;
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
    upperSolve!T(R, gm);
    auto Vm = V.sliceDup(0, n, 0, m);
    nm.bbla.dot!T(Vm, gm, x);

    foreach (i; 0 .. n) x[i] += x0[i];

    multiply!T(A, x, Ax0);
    foreach (i; 0 .. n) {
        r0[i] = b[i] - Ax0[i];
    }

    return x.dup();
}

struct GMRESWorkSpace(T) {
    size_t n, m;
    T[] Ax0, r0, v, w, Pv, g0, g1, h, hR;
    Matrix!T V, H0, H1, Gamma, Q0, Q1;

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
        this.V = new Matrix!T(n, m+1);
        this.H0 = new Matrix!T(m+1, m); 
        this.H1 = new Matrix!T(m+1, m); 
        this.Gamma = new Matrix!T(m+1, m+1); 
        this.Q0 = new Matrix!T(m+1, m+1);
        this.Q1 = new Matrix!T(m+1, m+1);
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
void rpcGMRES(T)(SMatrix!T A, SMatrix!T P, T[] b, T[] x0, T[] x,
              int maxIters, double residTol, ref GMRESWorkSpace!T gws)
in {
    assert(A.ia.length-1 == b.length);
    assert(b.length == x0.length);
    assert(maxIters >= 1);
}
body {
    T resid;
    size_t n = b.length;
    size_t m = maxIters;
    // 0. Initialise the values in some of the pre-allocated storage
    foreach( idx; 0..gws.g0.length) gws.g0[idx] = 0.0;
    foreach( idx; 0..gws.g1.length) gws.g1[idx] = 0.0;
    gws.H0.zeros();
    gws.H1.zeros();
    gws.Gamma.eye();

    // 1. Compute r0, beta, v1
    multiply!T(A, x0, gws.Ax0);
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
        solve!T(P, gws.Pv);
        multiply!T(A, gws.Pv, gws.w);
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
            nm.bbla.dot!T(gws.Q0, j+1, j+1, gws.h, gws.hR);
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
        nm.bbla.dot!T(gws.Gamma, j+2, j+2, gws.H0, j+1, gws.H1);
        nm.bbla.dot!T(gws.Gamma, j+2, j+2, gws.g0, gws.g1);
        // Accumulate Gamma rotations in Q.
        if ( j == 0 ) {
            copy(gws.Gamma, gws.Q1);
        }
        else {
            nm.bbla.dot!T(gws.Gamma, j+2, j+2, gws.Q0, j+2, gws.Q1);
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
    upperSolve!T(gws.H1, to!int(m), gws.g1);
    nm.bbla.dot!T(gws.V, n, m, gws.g1, x);
    solve!T(P, x);

    foreach (i; 0 .. n) x[i] += x0[i];
}

T[] fGMRES(T)(SMatrix!T A, SMatrix!T P, T[] b, T[] x0, int mOuter, int mInner)
in {
    assert(A.ia.length-1 == b.length);
    assert(b.length == x0.length);
    assert(mOuter >= 1);
}
body {
    immutable double ZERO_TOL = 1.0e-15;
    // 0. Allocate working matrices and vectors
    size_t n = b.length;
    T[] Ax0, r0, v, w, g_old, g, x, z, x0_inner;
    Ax0.length = n;
    r0.length = n;
    v.length = n;
    z.length = n;
    w.length = n;
    x.length = n;
    x0_inner.length = n;
    foreach( idx; 0..x0_inner.length) x0_inner[idx] = 0.0;
    g_old.length = mOuter+1;
    g.length = mOuter+1;
    auto H = new Matrix!T(mOuter+1, mOuter);
    H.zeros();
    auto Hold = new Matrix!T(mOuter+1, mOuter);
    auto Gamma = new Matrix!T(mOuter+1, mOuter+1);
    auto V = new Matrix!T(n, mOuter+1);
    auto Z = new Matrix!T(n, mOuter+1);

    // 1. Compute r0, beta, v1
    multiply!T(A, x0, Ax0);
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
        z = gmres!T(P, v, x0_inner, mInner);
        // Save z vector in Z
        foreach (k; 0 .. n ) Z[k,j] = z[k];
        multiply!T(A, z, w);
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
    foreach( idx; 0..g_old.length) g_old[idx] = 0.0;
    g_old[0] = beta;
    copy(H, Hold);
    foreach (i; 0 .. mOuter) {
        T c_i, s_i, denom;
        denom = sqrt(H[i,i]*H[i,i] + H[i+1,i]*H[i+1,i]);
        s_i = H[i+1,i]/denom; 
        c_i = H[i,i]/denom;
        Gamma.eye();
        Gamma[i,i] = c_i; Gamma[i,i+1] = s_i;
        Gamma[i+1,i] = -s_i; Gamma[i+1,i+1] = c_i;
        nm.bbla.dot!T(Gamma, Hold, H);
        nm.bbla.dot!T(Gamma, g_old, g);
        // Get Qold and g_old ready for next step
        g_old[] = g[];
        copy(H, Hold);
    }

    auto resid = fabs(g[mOuter]);
    auto R = H.sliceDup(0, mOuter, 0, mOuter);
    auto gm = g[0 .. mOuter];
    // At end H := R
    //        g := gm
    upperSolve!T(R, gm);
    auto Zm = Z.sliceDup(0, n, 0, mOuter);
    nm.bbla.dot!T(Zm, gm, x);

    foreach (i; 0 .. n) x[i] += x0[i];

    return x.dup();
}



version(smla_test) {
    import util.msg_service;
    int main() {
        // This example matrix is from Saad (2003), Sec. 3.4
        auto a = new SMatrix!number([to!number(1.), to!number(2.), to!number(3.), to!number(4.),
                                     to!number(5.), to!number(6.), to!number(7.), to!number(8.),
                                     to!number(9.), to!number(10.), to!number(11.), to!number(12.)],
                                    [0, 3, 0, 1, 3, 0, 2, 3, 4, 2, 3, 4],
                                    [0, 2, 5, 9, 11, 12]);
        // Test construction of matrix row by row
        auto b = new SMatrix!number();
        b.addRow([to!number(1.), to!number(2.)], [0, 3]);
        b.addRow([to!number(3.), to!number(4.), to!number(5.)], [0, 1, 3]);
        b.addRow([to!number(6.), to!number(7.), to!number(8.), to!number(9.)], [0, 2, 3, 4]);
        b.addRow([to!number(10.), to!number(11.)], [2, 3]);
        b.addRow([to!number(12.)], [4]);
        assert(approxEqualMatrix(a, b), failedUnitTest());
        // Test matrix multiply
        number[] v = [to!number(1.), to!number(2.), to!number(3.), to!number(4.), to!number(5.)];
        number[] c;
        c.length = v.length;
        multiply(a, v, c);
        number[] expected_c = [to!number(9.), to!number(31.), to!number(104.), to!number(74.), to!number(60.)];
        foreach ( i; 0 .. c.length ) {
            assert(approxEqualNumbers(c[i], expected_c[i]), failedUnitTest());
        }
        // Test decompILU0
        // This example matrix and decomposition is from Gerard and Wheatley, 6th ed, Sec. 2.4
        auto d = new SMatrix!number();
        d.addRow([to!number(-4.), to!number(2.)], [0, 1]);
        d.addRow([to!number(1.), to!number(-4.), to!number(1.)], [0, 1, 2]);
        d.addRow([to!number(1.), to!number(-4.), to!number(1.)], [1, 2, 3]);
        d.addRow([to!number(1.), to!number(-4.), to!number(1.)], [2, 3, 4]);
        d.addRow([to!number(2.), to!number(-4.)], [3, 4]);
        decompILU0(d);
        auto dLU = new SMatrix!number();
        dLU.addRow([to!number(-4.), to!number(2.)], [0, 1]);
        dLU.addRow([to!number(-0.25), to!number(-3.5), to!number(1.)], [0, 1, 2]);
        dLU.addRow([to!number(-0.2857), to!number(-3.7143), to!number(1.)], [1, 2, 3]);
        dLU.addRow([to!number(-0.2692), to!number(-3.7308), to!number(1.)], [2, 3, 4]);
        dLU.addRow([to!number(-0.5361), to!number(-3.4639)], [3, 4]);
        assert(approxEqualMatrix(d, dLU), failedUnitTest());
        // Test solve.
        // Let's give the ILU(0) method a triangular matrix that is can solve exactly.
        // This example is taken from Faires and Burden (1993), Sec. 6.6, example 3. 
        auto e = new SMatrix!number();
        e.addRow([to!number(2.), to!number(-1.)], [0, 1]);
        e.addRow([to!number(-1.), to!number(2.), to!number(-1.)], [0, 1, 2]);
        e.addRow([to!number(-1.), to!number(2.), to!number(-1.)], [1, 2, 3]);
        e.addRow([to!number(-1.), to!number(2.)], [2, 3]);
        decompILU0(e);
        number[] B = [to!number(1.), to!number(0.), to!number(0.), to!number(1.)];
        solve(e, B);
        number[] B_exp = [to!number(1.), to!number(1.), to!number(1.), to!number(1.)];
        foreach ( i; 0 .. B.length ) {
            assert(approxEqualNumbers(B[i], B_exp[i]), failedUnitTest());
        }
        // Now let's see how we go at an approximate solve by using a non-triangular matrix.
        // This is example 2.2 from Gerard and Wheatley, 6th edition
        auto f = new SMatrix!number();
        f.addRow([to!number(3.), to!number(2.), to!number(-1.), to!number(2.)], [0, 1, 2, 3]);
        f.addRow([to!number(1.), to!number(4.), to!number(2.)], [0, 1, 3]);
        f.addRow([to!number(2.), to!number(1.), to!number(2.), to!number(-1.)], [0, 1, 2, 3]);
        f.addRow([to!number(1.), to!number(1.), to!number(-1.), to!number(3.)], [0, 1, 2, 3]);
        decompILU0(f);
        number[] C = [to!number(2.), to!number(2.), to!number(0.), to!number(0.)];
        solve(f, C);
        number[] C_exp = [to!number(0.333333), to!number(0.666667), to!number(-1), to!number(-0.666667)];
        foreach ( i; 0 .. C.length ) {
            assert(approxEqualNumbers(C[i], C_exp[i]), failedUnitTest());
        }
        
        // Let's test the ILU(p) method
        auto s = new SMatrix!number([to!number(1.), to!number(1.), to!number(4.), to!number(2.),
                                     to!number(4.), to!number(1.), to!number(2.), to!number(1.),
                                     to!number(8.), to!number(2.), to!number(4.), to!number(1.),
                                     to!number(3.), to!number(6.), to!number(2.), to!number(1.)],
                                    [0, 1, 4, 1, 2, 4, 0, 1, 2, 3, 2, 3, 0, 1, 2, 4],
                                    [0, 3, 6, 10, 12, 16]);
        int p;
	/*
        // test for ILU(p=0)
        p = 0;
        decompILUp(s, p);
        auto sol0 = new SMatrix!number([1., 1., 4., 2., 4., 1., 2., -0.5, 10., 2., 0.4, 0.2, 3., 1.5, -0.4, -12.5],
                                [0, 1, 4, 1, 2, 4, 0, 1, 2, 3, 2, 3, 0, 1, 2, 4],
                                [0, 3, 6, 10, 12, 16]);

        assert(approxEqualMatrix(s, sol0), failedUnitTest());
        
        // As a result of the note in decompILUp() we don't expect an exact match of the SMatrix classes
        foreach ( i; 0 .. 5) {
            foreach ( j; 0 .. 5) {
                assert(approxEqual(s[i,j], sol0[i,j]), failedUnitTest());
            }
        }
        */
        // test for ILU(p=2)
        s = new SMatrix!number([to!number(1.), to!number(1.), to!number(4.), to!number(2.),
                                to!number(4.), to!number(1.), to!number(2.), to!number(1.),
                                to!number(8.), to!number(2.), to!number(4.), to!number(1.),
                                to!number(3.), to!number(6.), to!number(2.), to!number(1.)],
                               [0, 1, 4, 1, 2, 4, 0, 1, 2, 3, 2, 3, 0, 1, 2, 4],
                               [0, 3, 6, 10, 12, 16]);
        p = 2;
        decompILUp(s, p);
        auto sol1 = new SMatrix!number([to!number(1.), to!number(1.), to!number(4.), to!number(2.),
                                        to!number(4.), to!number(1.), to!number(2.), to!number(-0.5),
                                        to!number(10.), to!number(2.), to!number(-7.5), to!number(0.4),
                                        to!number(0.2), to!number(3.), to!number(3.), to!number(1.5),
                                        to!number(-0.4), to!number(4.), to!number(-27.5)],
                                       [0, 1, 4, 1, 2, 4, 0, 1, 2, 3, 4, 2, 3, 4, 0, 1, 2, 3, 4],
                                       [0, 3, 6, 11, 14, 19]);
	assert(approxEqualMatrix!number(s, sol1), failedUnitTest());
        /*
        // As a result of the note in decompILUp() we don't expect an exact match of the SMatrix classes
        foreach ( i; 0 .. 5) {
            foreach ( j; 0 .. 5) {
                assert(approxEqual(s[i,j], sol1[i,j]), failedUnitTest());
            }
        }
        */
        // Test GMRES on Faires and Burden problem.

        auto g = new SMatrix!number();
        g.addRow([to!number(2.), to!number(-1.)], [0, 1]);
        g.addRow([to!number(-1.), to!number(2.), to!number(-1.)], [0, 1, 2]);
        g.addRow([to!number(-1.), to!number(2.), to!number(-1.)], [1, 2, 3]);
        g.addRow([to!number(-1.), to!number(2.)], [2, 3]);
        number[] B1 = [to!number(1.), to!number(0.), to!number(0.), to!number(1.)];
        number[] x0 = [to!number(1.2), to!number(0.8), to!number(0.9), to!number(1.1)];

        
        auto x = gmres!number(g, B1, x0, 4);
        foreach (i; 0 .. x.length) {
            assert(approxEqualNumbers(x[i], B_exp[i]), failedUnitTest());
        }
        
        x = gmres2!number(g, B1, x0, 5, 1.0e-10);
        foreach (i; 0 .. x.length) {
            assert(approxEqualNumbers(x[i], B_exp[i]), failedUnitTest());
        }

        // Test pre-conditioned GMRES on Gerald and Wheatley problem.
        // This time we expect the exact answer. Earlier we only used
        // an incomplete LU factorisation and so the result was
        // only approximate.
        auto h = new SMatrix!number();
        h.addRow([to!number(3.), to!number(2.), to!number(-1.), to!number(2.)], [0, 1, 2, 3]);
        h.addRow([to!number(1.), to!number(4.), to!number(2.)], [0, 1, 3]);
        h.addRow([to!number(2.), to!number(1.), to!number(2.), to!number(-1.)], [0, 1, 2, 3]);
        h.addRow([to!number(1.), to!number(1.), to!number(-1.), to!number(3.)], [0, 1, 2, 3]);
        auto Ph = new SMatrix!number(h);
        decompILU0!number(Ph);
        number[] C1 = [to!number(2.), to!number(2.), to!number(0.), to!number(0.)];
        x0 = [to!number(0.2), to!number(0.5), to!number(-1.1), to!number(-0.6)];
        int maxIters = 5;
        auto gws = GMRESWorkSpace!number(x0.length, maxIters);
        rpcGMRES!number(h, Ph, C1, x0, x, maxIters, 1.0e-15, gws);
        number[] C1_exp = [to!number(0.273), to!number(0.773), to!number(-1.0), to!number(-0.682)];

        foreach (i; 0 .. x.length) {
            assert(approxEqualNumbers(x[i], C1_exp[i]), failedUnitTest());
        }

        auto Pi = new SMatrix!number();
        Pi.addRow([to!number(1.),], [0]);
        Pi.addRow([to!number(1.),], [1]);
        Pi.addRow([to!number(1.),], [2]);
        Pi.addRow([to!number(1.),], [3]);

        x = fGMRES!number(h, Pi, C1, x0, 5, 3);
        foreach (i; 0 .. x.length) {
            assert(approxEqualNumbers(x[i], C1_exp[i]), failedUnitTest());
        }

        // This example tests the addition of values to zero-entries
        auto z = new SMatrix!number([to!number(1.), to!number(2.), to!number(3.), to!number(4.),
                                     to!number(5.), to!number(6.), to!number(7.), to!number(8.),
                                     to!number(9.), to!number(10.), to!number(11.), to!number(12.)],
                                    [0, 1, 2, 0, 1, 3, 0, 2, 3, 1, 2, 3],
                                    [0, 3, 6, 9, 12]);
        z[0,3] = to!number(99.0);
        z[1,2] = to!number(99.0);
        z[2,1] = to!number(99.0);
        z[3,0] = to!number(99.0);
        auto w = new SMatrix!number([to!number(1.), to!number(2.), to!number(3.), to!number(99.),
                                     to!number(4.), to!number(5.), to!number(99.), to!number(6.),
                                     to!number(7.), to!number(99.), to!number(8.), to!number(9.),
                                     to!number(99.), to!number(10.), to!number(11.), to!number(12.)],
                                    [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3],
                                    [0, 4, 8, 12, 16]);
        assert(approxEqualMatrix!number(z, w), failedUnitTest());
        return 0;
    }
}
