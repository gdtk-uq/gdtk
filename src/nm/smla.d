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
import std.algorithm : reduce, fill;
import std.algorithm.sorting : sort;
import std.algorithm.searching : maxElement;
import std.typecons : Tuple;
import nm.bbla;
import ntypes.complex;
import nm.number;
import nm.stmatrix;

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

    this(STMatrix!T stMat)
    {
        // Use the average number of non-zeros per row to reserve an array size.
        // This might help cut down on allocations while assembling a row.
        auto n_nzAvgPerRow = stMat.n_nzeros / stMat.n_rows;
        T[] ai;
        size_t[] ji;
        T[size_t] vi;
        Tuple!(size_t, size_t)[] keysToRemove;
        // We loop over all rows.
        // Every row should have at least one non-zero entry, otherwise there's an issue with our matrix.
        foreach (i; 0 .. stMat.n_rows) {
            foreach (c, v; stMat.val) {
                if (c[0] == i) {
                    vi[c[1]] = v;
                    keysToRemove ~= c;
                }
            }
            // Add the collected values as a new row in the sparse matrix.
            // We need to sort the entries in column-ascending order.

            ji = vi.keys;
            sort(ji);
            foreach (j; ji) ai ~= vi[j];
            addRow(ai, ji);
            // Remove the entries from the sparse triplet matrix that we just used
            // to save looking through them again.
            foreach (c; keysToRemove) {
                stMat.val.remove(c);
            }
            // Prepare temporary containers for next loop.
            ai.length = 0;
            ji.length = 0;
            keysToRemove.length = 0;
            vi.clear;
        }
    }

    void addRow(T[] ai, size_t[] ji)
    {
        if ( ia.length == 0 )
            ia ~= 0;
        aa ~= ai[];
        ja ~= ji[];
        ia ~= aa.length; // Now ia is already ready for next addition.
    }

    void dropRows(size_t nDrop)
    {
        if (nDrop > (ia.length-1)) {
            nDrop = ia.length-1;
        }
        size_t nvals = ia[$-1] - ia[$-(1+nDrop)];
        ja.length = ja.length - nvals;
        aa.length = aa.length - nvals;
        ia.length = ia.length - nDrop;
    }

    void scaleRow(size_t row, T scaleFactor)
    {
        assert(row < ia.length-1);
        foreach (j; ia[row] .. ia[row+1]) {
            aa[j] *= scaleFactor;
        }
    }

    void scaleCol(size_t col, T scaleFactor)
    {
        assert(col < ia.length-1);
        foreach (row; 0..ia.length-1) {
            foreach (offset, j; ja[ia[row]..ia[row+1]]) {
                if ( j == col ) {
                    aa[ia[row]+offset] *= scaleFactor;
                    break;
                }
            }
        }
    }

    void scale(T scaleFactor) {
        foreach (ref elem; aa) {
            elem *= scaleFactor;
        }
    }

    const T opIndex(size_t row, size_t col)
    {
        // We need to search the given row to see if an entry
        // is present in the given column.
        foreach (offset, j; ja[ia[row]..ia[row+1]]) {
            if ( j == col )
                return aa[ia[row]+offset];
        }
        // else search failed, so we have a 0.0 entry
        return to!T(0.0);
    }

    // We now allow assingment to zero entries.
    ref T opIndexAssign(T c, size_t row, size_t col) {
        foreach (offset, j; ja[ia[row]..ia[row+1]]) {
            if ( j == col ) {
                aa[ia[row]+offset] = c;
                return aa[ia[row]+offset];
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
        auto jmax = ja.maxElement;
        string s = "SMatrix[\n";
        foreach (row; 0 .. ia.length-1) {
            foreach (col; 0 .. jmax+1) {
                s ~= to!string(this[row,col]);
                if ( col < jmax )
                    s ~= ", ";
            }
            s ~= "\n";
        }
        s ~= "]";
        return s;
    }

    SMatrix!T transpose()
    {
	// Presently, assume square matrices
	auto tp = new SMatrix!T();
	tp.aa.length = aa.length;
	tp.ja.length = ja.length;
	tp.ia.length = ia.length;

	tp.ia[] = 0;

	// Treat original matrix a CSC
	// We'll count entries in the transposed rows
	// by looking at "columns" of original.
	foreach (row_idx; ja) {
	    tp.ia[row_idx]++;
	}

	// To get the offsets into ia, compute a running sum
	size_t sum = 0;
	foreach (row; 0 .. ia.length) {
	    size_t tmp = tp.ia[row];
	    tp.ia[row] = sum;
	    sum += tmp;
	}

	foreach (col; 0 .. ia.length-1) {
	    foreach (jj; ia[col] .. ia[col+1]) {
		auto row = ja[jj];
		auto dest = tp.ia[row];

		tp.ja[dest] = col;
		tp.aa[dest] = aa[jj];

		tp.ia[row]++;
	    }
	}

	size_t last = 0;
	foreach (row; 0 .. ia.length) {
	    auto tmp = tp.ia[row];
	    tp.ia[row] = last;
	    last = tmp;
	}

	return tp;

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

void multiply(T)(SMatrix!T a, T[] b, T[] c, size_t nrows=0, size_t ncols=0)
{
    size_t k0, k1;
    if (nrows == 0) { // just multiply using every row in matrix
        nrows = a.ia.length-1;
    }
    if (nrows > a.ia.length-1) { // don't fail, but reset to maximum available rows
        nrows = a.ia.length-1;
    }
    foreach (i; 0 .. nrows) {
        k0 = a.ia[i];
        k1 = a.ia[i+1];
        c[i] = 0.0;
        foreach (k; k0 .. k1) {
            if (ncols >= 1 && a.ja[k] >= ncols) continue;
            c[i] += a.aa[k] * b[a.ja[k]];
        }
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
    /*
    Solve a linear system using a prefactored sparse LU decomposition.

    Refactor by NNG, 18/08/23, Tangalooma QLD
    */

    // Note that ia has an extra element at the end, so the number of rows
    // n is actually ia.length-1. Meanwhile nn is the length of the under
    // -lying data array LU.aa, i.e. the number of nonzero entries in the matrix
    size_t n = LU.ia.length-1;
    size_t nn = LU.ja.length;
    assert(b.length == n);

    // Forward elimination
    foreach ( i; 1 .. n ) {
        foreach ( ii; LU.ia[i] .. LU.ia[i+1] ) {
            // Only work up to i
            size_t j = LU.ja[ii];
            if ( j >= i ) break;
            T multiplier = LU.aa[ii];
            b[i] -= multiplier * b[j];
        }
    }

    // Back substitution
    b[n-1] /= LU.aa[nn-1];
    for ( int i = to!int(n-2); i >= 0; --i ) {
        T sum = b[i];
        size_t diag_index;
        foreach ( ii; LU.ia[i] .. LU.ia[i+1] ) {
            size_t j = LU.ja[ii];
            // If we run across the diagonal element keep a note of where it is
            if ( j == i ) {
                diag_index = ii;
                continue;
            // Otherwise only work with j >= i+1
            } else if( j > i ) {
                sum -= LU.aa[ii] * b[j];
            }
        }
        b[i] = sum/LU.aa[diag_index];
    }
}


/**
 * This algorithm inverts the diagonal blocks of the sparse matrix A
 * the inverted blocks are stored in place, i.e. A = L + D^(-1) U
 * this can be useful for some of the iterative linear solvers such as Jacobi and Gauss-Seidel
 **/
void invert_block_diagonal(T)(SMatrix!T A, Matrix!T D, Matrix!T Dinv, size_t nblks, size_t blk_size) {

    foreach (k; 0..nblks) {
        // grab the current diagonal block
        foreach (i; 0..blk_size) {
            size_t idx = k*blk_size + i;
            foreach (j; 0..blk_size) {
                size_t jdx = k*blk_size + j;
                D[i,j] = A[idx,jdx];
            }
        }
        // invert the block matrix
        Dinv = inverse(D);
        // replace the diagonal block in the sparse matrix with its inverse
        foreach (i; 0..blk_size) {
            size_t idx = k*blk_size + i;
            foreach (j; 0..blk_size) {
                size_t jdx = k*blk_size + j;
                A[idx,jdx] = Dinv[i,j];
            }
        }
    }

} // end invert_block_diagonal()

/**
 * This algorithm multiplies only the block diagonal portion of the
 * sparse matrix (A) with a vector (b) and stores the result in (c)
 * this can be useful for some of the iterative linear solvers such as Jacobi and Gauss-Seidel
 **/
void multiply_block_diagonal(T)(SMatrix!T a, T[] b, T[] c, size_t nblks, size_t blk_size)
in {
    assert(a.ia.length-1 == c.length);
    // Some faith is given that the user provides an appropriate length
    // input vector b, since in CRS we cannot easily determine the number of columns.
}
do {
    size_t i0, i1, j0, j1;
    foreach (nb; 0..nblks) {
        i0 = nb*blk_size;
        i1 = nb*blk_size+blk_size;
        foreach ( i; i0..i1 ) {
            j0 = a.ia[i];
            j1 = a.ia[i+1];
            foreach ( j; a.ja[j0 .. j1] ) {
                // only work on entries that are in the diagonal block
                if ( j < i0 )
                    continue;
                if ( j > i1-1 )
                    break;
                c[i] += a[i,j] * b[j];
            }
        }
    }
 } // end multiply_block_diagonal()

/**
 * This algorithm multiplies only the block lower triangular portion of the
 * sparse matrix (A) with a vector (b) and stores the result in (c)
 * this is useful for some of the iterative linear solvers such as Jacobi and Gauss-Seidel
 **/
void multiply_block_lower_triangular(T)(SMatrix!T a, T[] b, T[] c, size_t nblks, size_t blk_size)
in {
    assert(a.ia.length-1 == c.length);
    // Some faith is given that the user provides an appropriate length
    // input vector b, since in CRS we cannot easily determine the number of columns.
}
do {
    size_t i0, i1, j0, j1;
    foreach (nb; 0..nblks) {
        i0 = nb*blk_size;
        i1 = nb*blk_size+blk_size;
        foreach ( i; i0..i1 ) {
            j0 = a.ia[i];
            j1 = a.ia[i+1];
            foreach ( j; a.ja[j0 .. j1] ) {
                // only work on entries where j < diagonal block
                if ( j >= i0 )
                    break;
                c[i] += a[i,j] * b[j];
            }
        }
    }
 } // end multiply_block_lower_triangular()

/**
 * This algorithm performs a block lower triangular solve
 * the routine assumes that the sparse matrix is in the form A = L + D^(-1) + U
 * this method can be used for iterative solvers such as Gauss-Seidel
 **/
void block_lower_triangular_solve(T)(SMatrix!T a, T[] b, T[] c, size_t nblks, size_t blk_size)
in {
    assert(a.ia.length-1 == c.length);
    // Some faith is given that the user provides an appropriate length
    // input vector b, since in CRS we cannot easily determine the number of columns.
}
do {
    size_t i0, i1, j0, j1;
    foreach (nb; 0..nblks) {
        i0 = nb*blk_size;
        i1 = nb*blk_size+blk_size;
        // move the up to date entries to the right hand side
        foreach ( i; i0..i1 ) {
            j0 = a.ia[i];
            j1 = a.ia[i+1];
            foreach ( j; a.ja[j0 .. j1] ) {
                // only work on entries where j < diagonal block
                if ( j >= i0 )
                    break;
                b[i] -= a[i,j] * c[j];
            }
        }
        // multiply by inverted diagonal block
        foreach ( i; i0..i1 ) {
            c[i] = 0.0;
            j0 = a.ia[i];
            j1 = a.ia[i+1];
            foreach ( j; a.ja[j0 .. j1] ) {
                // only work on entries that are in the diagonal block
                if ( j < i0 )
                    continue;
                if ( j > i1-1 )
                    break;
                c[i] += a[i,j] * b[j];
            }
        }
    }
 } // end block_lower_triangular_solve()

/**
 * This algorithm multiplies only the block upper triangular portion of the
 * sparse matrix (A) with a vector (b) and stores the result in (c)
 * this is useful for some of the iterative linear solvers such as Jacobi and Gauss-Seidel
 **/
void multiply_block_upper_triangular(T)(SMatrix!T a, T[] b, T[] c, size_t nblks, size_t blk_size)
in {
    assert(a.ia.length-1 == c.length);
    // Some faith is given that the user provides an appropriate length
    // input vector b, since in CRS we cannot easily determine the number of columns.
}
do {
    size_t i0, i1, j0, j1;
    foreach (nb; 0..nblks) {
        i0 = nb*blk_size;
        i1 = nb*blk_size+blk_size;
        foreach ( i; i0..i1 ) {
            j0 = a.ia[i];
            j1 = a.ia[i+1];
            foreach ( j; a.ja[j0 .. j1] ) {
                // only work on entries where j > diagonal block
                if ( j <= i1-1 )
                    continue;
                c[i] += a[i,j] * b[j];
            }
        }
    }
 } // multiply_block_upper_triangular()

/**
 * This algorithm performs a block upper triangular solve
 * the routine assumes that the sparse matrix is in the form A = L + D^(-1) + U
 * this method can be used for iterative solvers such as Gauss-Seidel
 **/
void block_upper_triangular_solve(T)(SMatrix!T a, T[] b, T[] c, size_t nblks, size_t blk_size)
in {
    assert(a.ia.length-1 == c.length);
    // Some faith is given that the user provides an appropriate length
    // input vector b, since in CRS we cannot easily determine the number of columns.
}
do {
    size_t i0, i1, j0, j1;
    for ( int nb = to!int(nblks)-1; nb >= 0; --nb ) {
        i0 = to!size_t(nb)*blk_size;
        i1 = to!size_t(nb)*blk_size+blk_size;
        // move the up to date entries to the right hand side
        foreach ( i; i0..i1 ) {
            j0 = a.ia[i];
            j1 = a.ia[i+1];
            foreach ( j; a.ja[j0 .. j1] ) {
                // only work on entries where j > diagonal block
                if ( j <= i1-1 )
                    continue;
                b[i] -= a[i,j] * c[j];
            }
        }
        // multiply by inverted diagonal block
        foreach ( i; i0..i1 ) {
            c[i] = 0.0;
            j0 = a.ia[i];
            j1 = a.ia[i+1];
            foreach ( j; a.ja[j0 .. j1] ) {
                // only work on entries that are in the diagonal block
                if ( j < i0 )
                    continue;
                if ( j > i1-1 )
                    break;
                c[i] += a[i,j] * b[j];
            }
        }
    }
 } // end block_upper_triangular_solve()

T[] gmres(T)(SMatrix!T A, T[] b, T[] x0, int m)
in {
    assert(A.ia.length-1 == b.length);
    assert(b.length == x0.length);
    assert(m >= 1);
}
do {
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
do {
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
 * A restarted GMRES iterative solver with right preconditioning.
 *
 * Params:
 *    A =           coefficient matrix
 *    P =           pre-conditioning matrix in LU form
 *    b =           RHS vector
 *    x0 =          initial guess for solution
 *    x =           solution vector place in x
 *    maxIters =    maximum number of iterations
 *    residTol =    stop iterations when residual below this tolerance
 *    gws =         pre-allocated workspace
 *    maxRestarts = maximum number of restarts
 */
void rpcGMRES(T)(SMatrix!T A, SMatrix!T P, T[] b, T[] x0, T[] x,
                 int maxIters, double residTol, ref GMRESWorkSpace!T gws, int maxRestarts=1)
in {
    assert(A.ia.length-1 == b.length);
    assert(b.length == x0.length);
    assert(maxIters >= 1);
}
do {
    T resid;
    size_t n = b.length;
    size_t m = maxIters;

    // Perform restarts
    for (size_t r = 0; r < maxRestarts; r++) {

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
        //writefln("beta: %.16e", beta);
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
            //writefln("restart: %d, iteration: %d, residual: %.8e", r, j, resid);
            if ( resid <= residTol ) {
                m = j+1;
                break;
            }

            // Prepare for next iteration
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

        if ( resid <= residTol || r+1 == maxRestarts ) {
            break;
        }

        // Else, prepare for a restart by setting the inital
        // guess to the current best estimate of the solution
        foreach (i; 0 .. n) x0[i] = x[i];
    }
}

T[] fGMRES(T)(SMatrix!T A, SMatrix!T P, T[] b, T[] x0, int mOuter, int mInner)
in {
    assert(A.ia.length-1 == b.length);
    assert(b.length == x0.length);
    assert(mOuter >= 1);
}
do {
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
        // Test dropping rows
        auto aa = new SMatrix!number(a);
        aa.dropRows(2);
        auto aaRef = new SMatrix!number([to!number(1.), to!number(2.), to!number(3.), to!number(4.),
                                         to!number(5.), to!number(6.), to!number(7.), to!number(8.),
                                         to!number(9.)],
                                        [0, 3, 0, 1, 3, 0, 2, 3, 4],
                                        [0, 2, 5, 9]);
        assert(approxEqualMatrix(aa, aaRef), failedUnitTest());
        // Test matrix multiply
        number[] v = [to!number(1.), to!number(2.), to!number(3.), to!number(4.), to!number(5.)];
        number[] c;
        c.length = v.length;
        multiply(a, v, c);
        number[] expected_c = [to!number(9.), to!number(31.), to!number(104.), to!number(74.), to!number(60.)];
        foreach ( i; 0 .. c.length ) {
            assert(approxEqualNumbers(c[i], expected_c[i]), failedUnitTest());
        }
	// Test transpose
	auto cc = new SMatrix!double();
	cc.addRow([1.0, 3.0], [0, 2]);
	cc.addRow([4.0], [2]);
	cc.addRow([7.0], [1]);
	auto ccT = cc.transpose();

	double[] ccT_aa_exp = [1.0, 7.0, 3.0, 4.0];
	size_t[] ccT_ja_exp = [0, 2, 0, 1];
	size_t[] ccT_ia_exp = [0, 1, 2, 4];
	foreach (i; 0 .. ccT.aa.length) {
	    assert(isClose(ccT.aa[i], ccT_aa_exp[i]));
	}
	foreach (i; 0 .. ccT.ja.length) {
	    assert(ccT.ja[i] == ccT_ja_exp[i]);
	}
	foreach (i; 0 .. ccT.ia.length) {
	    assert(ccT.ia[i] == ccT_ia_exp[i]);
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

        // Test pre-conditioned restarted GMRES on the linear system problem
        // outlined in Section III from:
        // Implementation of GMRES for chemically reacting flows by Maclean and White, 2012.
        auto mat_A = new SMatrix!number([to!number(1.), to!number(2.), to!number(-1.), to!number(3.), to!number(2.), to!number(-1.), to!number(-2.),
                                         to!number(2.), to!number(3.), to!number(-2.), to!number(-1.), to!number(2.), to!number(4.), to!number(2.),
                                         to!number(-2.), to!number(1.), to!number(5.), to!number(-1.), to!number(-1.), to!number(6.), to!number(-2.),
                                         to!number(-2.), to!number(3.), to!number(-1.), to!number(-1.), to!number(-5.), to!number(4.), to!number(3.),
                                         to!number(-2.), to!number(1.), to!number(2.), to!number(1.), to!number(-1.), to!number(3.), to!number(4.)],
                                        [0, 1, 5, 0, 1, 2, 6, 1, 2, 3, 7, 2, 3, 4, 8, 3, 4, 9, 0, 5, 1, 5, 6, 7, 2, 6, 7, 8, 3, 7, 8, 9, 4, 8, 9],
                                        [0, 3, 7, 11, 15, 18, 20, 24, 28, 32, 35]);
        auto pc_A = new SMatrix!number(mat_A);
        decompILU0!number(pc_A);
	number[] vec_b   = [to!number(1.), to!number(2.), to!number(3.), to!number(4.), to!number(5.),
                            to!number(6.), to!number(7.), to!number(8.), to!number(9.), to!number(10.)];
        number[] vec_x0  = [to!number(0.), to!number(0.), to!number(0.), to!number(0.), to!number(0.),
                            to!number(0.), to!number(0.), to!number(0.), to!number(0.), to!number(0.)];
        number[] vec_x   = [to!number(0.), to!number(0.), to!number(0.), to!number(0.), to!number(0.),
                            to!number(0.), to!number(0.), to!number(0.), to!number(0.), to!number(0.)];
        number[] vec_ref = [to!number(5.29051), to!number(-1.20438), to!number(4.15595), to!number(2.22681), to!number(0.0574555),
                            to!number(1.88175), to!number(3.65341), to!number(2.60547), to!number(6.66703), to!number(-2.48591)];
        int max_iterations = 5;
        int max_restarts   = 20;
        auto gws2 = GMRESWorkSpace!number(vec_x0.length, max_iterations);
        rpcGMRES!number(mat_A, pc_A, vec_b, vec_x0, vec_x, max_iterations, 1.0e-14, gws2, max_restarts);
        foreach (i; 0 .. vec_x.length) {
            assert(approxEqualNumbers(vec_x[i], vec_ref[i]), failedUnitTest());
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

        // Test initialisation from a sparse triplet matrix
        string fName = "test_data/b1_ss.mtx";
        auto matrix = readFromMatrixMarketFile!double(fName);
        auto testMat = new SMatrix!double(matrix);
        assert(approxEqualNumbers(testMat[4,0], -0.03599942, 1.0e-7), failedUnitTest());
        assert(approxEqualNumbers(testMat[6,6], 1.0, 1.0e-7), failedUnitTest());

        // We will test the suite of block sparse matrix methods by solving a linear system using two common iterative solution methods
        size_t blk_size = 2;
        size_t nblks = 5;
        auto A = new SMatrix!number([to!number(10.), to!number(2.),  to!number(-1.),
                                     to!number(3.),  to!number(20.), to!number(-1.), to!number(-2.),
                                     to!number(2.),  to!number(30.), to!number(-2.), to!number(-1.),
                                     to!number(2.),  to!number(40.), to!number(2.),  to!number(-2.),
                                     to!number(1.),  to!number(50.), to!number(-1.),
                                     to!number(-1.), to!number(60.),
                                     to!number(-2.), to!number(-2.), to!number(30.), to!number(-1.),
                                     to!number(-1.), to!number(-5.), to!number(40.), to!number(3.),
                                     to!number(-2.), to!number(1.),  to!number(20.), to!number(1.),
                                     to!number(-1.), to!number(3.),  to!number(40.),],
                                    [0, 1, 5, 0, 1, 2, 6, 1, 2, 3, 7, 2, 3, 4, 8, 3, 4, 9, 0, 5, 1, 5, 6, 7, 2, 6, 7, 8, 3, 7, 8, 9, 4, 8, 9],
                                    [0, 3, 7, 11, 15, 18, 20, 24, 28, 32, 35]);
        number[] R = [to!number(1.), to!number(2.0), to!number(3.0), to!number(4.0), to!number(5.0), to!number(6.0), to!number(7.0), to!number(8.0), to!number(9.0), to!number(10.0)];
        number[] sol = [to!number(0.0865856), to!number(0.117794), to!number(0.106302), to!number(0.111582), to!number(0.102159),
                        to!number(0.101443), to!number(0.254665), to!number(0.201483), to!number(0.440107), to!number(0.219546)];

        // test multiplying only the block diagonal of a sparse matrix by a vector
        number[] DdotR;
        DdotR.length = R.length;
        DdotR[] = to!number(0.0);
        multiply_block_diagonal(A, R, DdotR, nblks, blk_size);
        number[] s0 = [to!number(14.0), to!number(43.0), to!number(82.0), to!number(166.0), to!number(250.0),
                         to!number(360.0), to!number(202.0), to!number(285.0), to!number(190.0), to!number(427.0)];
        foreach (i; 0 .. s0.length) {
            assert(approxEqualNumbers(DdotR[i], s0[i]), failedUnitTest());
        }

        // test multiplying only the block lower triangular portion of a sparse matrix by a vector
        number[] LdotR;
        LdotR.length = R.length;
        LdotR[] = to!number(0.0);
        nm.smla.multiply_block_lower_triangular(A, R, LdotR, nblks, blk_size);
        number[] s1 = [to!number(0.0), to!number(0.0), to!number(4.0), to!number(0.0), to!number(4.0),
                         to!number(-1.0), to!number(-16.0), to!number(-3.0), to!number(0.0), to!number(-5.0)];
        foreach (i; 0 .. s1.length) {
            assert(approxEqualNumbers(LdotR[i], s1[i]), failedUnitTest());
        }

        // test multiplying only the block upper triangular portion of a sparse matrix by a vector
        number[] UdotR;
        UdotR.length = R.length;
        UdotR[] = to!number(0.0);
        nm.smla.multiply_block_upper_triangular(A, R, UdotR, nblks, blk_size);
        number[] s2 = [to!number(-6.0), to!number(-17.0), to!number(-8.0), to!number(-8.0), to!number(-10.0),
                         to!number(0.0), to!number(0.0), to!number(27.0), to!number(0.0), to!number(0.0)];
        foreach (i; 0 .. s2.length) {
            assert(approxEqualNumbers(UdotR[i], s2[i]), failedUnitTest());
        }

        // invert block diagonal in place => A = L + D^(-1) + U
        Matrix!number D; Matrix!number Dinv;
        D = new Matrix!number(blk_size,blk_size);
        Dinv = new Matrix!number(blk_size,blk_size);
        invert_block_diagonal(A, D, Dinv, nblks, blk_size,);

        // initialise data arrays
        number[] rhs;
        rhs.length = R.length;
        number[] x_curr;
        x_curr.length = R.length;
        x_curr[] = to!number(0.0);

        // perform 8 iterations of the Jacobi method
        x_curr[] = to!number(0.0);
        foreach (k; 0..8) {
            rhs[] = to!number(0.0);
            multiply_block_upper_triangular(A, x_curr, rhs, nblks, blk_size);
            multiply_block_lower_triangular(A, x_curr, rhs, nblks, blk_size);
            rhs[] = R[] - rhs[];
            x_curr[] = to!number(0.0);
            multiply_block_diagonal(A, rhs, x_curr, nblks, blk_size);
        }
        foreach (i; 0 .. sol.length) {
            assert(approxEqualNumbers(x_curr[i], sol[i]), failedUnitTest());
        }

        // perform 4 iterations of the Symmetric Gauss-Seidel method
        x_curr[] = to!number(0.0);
        foreach (k; 0..4) {
            // forward sweep
            rhs[] = to!number(0.0);
            multiply_block_upper_triangular(A, x_curr, rhs, nblks, blk_size);
            rhs[] = R[] - rhs[];
            block_lower_triangular_solve(A, rhs, x_curr, nblks, blk_size);

            // backward sweep
            rhs[] = to!number(0.0);
            multiply_block_lower_triangular(A, x_curr, rhs, nblks, blk_size);
            rhs[] = R[] - rhs[];
            block_upper_triangular_solve(A, rhs, x_curr, nblks, blk_size);
        }
        foreach (i; 0 .. sol.length) {
            assert(approxEqualNumbers(x_curr[i], sol[i]), failedUnitTest());
        }

        return 0;
    }
}
