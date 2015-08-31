/**
 * bbla.d
 * Bare-bones linear algebra functions.
 *
 * Author: Peter J.
 * Version: 2014-06-28, just enough to get our CFD code going.
 *          2015-05-30, added least-squares solver
 */

module bbla;

import std.conv;
import std.algorithm;
import std.math;
import std.exception;

class Matrix {
    size_t _nrows;
    size_t _ncols;
    double[][] _data;

    this(size_t n) {
	_nrows = n;
	_ncols = n;
	_data.length = n;
	foreach(i; 0 .. n) _data[i].length = n;
    }

    this(size_t nrows, size_t ncols) {
	_nrows = nrows;
	_ncols = ncols;
	_data.length = nrows;
	foreach(i; 0 .. nrows) _data[i].length = ncols;
    }

    this(in Matrix other) {
	this(other._nrows, other._ncols);
	foreach(row; 0 .. _nrows)
	    foreach(col; 0 .. _ncols)
		_data[row][col] = other._data[row][col];
    }

    this(in double[] vec, string orient="column") 
    {
	if ( orient == "column" ) {
	    this(vec.length, 1);
	    foreach(row; 0 .. _nrows) _data[row][0] = vec[row];
	} else {
	    this(1, vec.length);
	    foreach(col; 0 .. _ncols) _data[0][col] = vec[col];
	}
    }

    this(in float[] vec, string orient="column") {
	double[] my_vec;
	my_vec.length = vec.length;
	foreach(i; 0 .. my_vec.length) my_vec[i] = vec[i]; 
	this(my_vec, orient);
    }

    this(in int[] vec, string orient="column") {
	double[] my_vec;
	my_vec.length = vec.length;
	foreach(i; 0 .. my_vec.length) my_vec[i] = vec[i]; 
	this(my_vec, orient);
    }

    this(in double[][] other) {
	this(other.length, other[0].length);
	foreach(row; 0 .. _nrows)
	    foreach(col; 0 .. _ncols)
		_data[row][col] = other[row][col];
    }

    this(in float[][] other) {
	this(other.length, other[0].length);
	foreach(row; 0 .. _nrows)
	    foreach(col; 0 .. _ncols)
		_data[row][col] = other[row][col];
    }

    this(in int[][] other) {
	this(other.length, other[0].length);
	foreach(row; 0 .. _nrows)
	    foreach(col; 0 .. _ncols)
		_data[row][col] = other[row][col];
    }

    @property const size_t nrows() { return _nrows; }
    @property const size_t ncols() { return _ncols; }

    const double opIndex(size_t row, size_t col) {
	return _data[row][col];
    }

    ref double opIndexAssign(double c, size_t row, size_t col) {
	_data[row][col] = c;
	return _data[row][col];
    }

    ref double opIndexOpAssign(string op)(double c, size_t row, size_t col)
	if ( op == "+" || op == "-" || op == "*" || op == "/" )
    {
	static if ( op == "+" )
	    _data[row][col] += c;
	else if ( op == "-" )
	    _data[row][col] -= c;
	else if ( op == "*" )
	    _data[row][col] *= c;
	else if ( op == "/" )
	    _data[row][col] /= c;
	return _data[row][col];
    }

    Matrix opBinary(string op)(in Matrix rhs)
	if ( op == "+" || op == "-" )
    {
	enforce(_nrows == rhs._nrows && _ncols == rhs._ncols,
		"incompatible matrices");
	Matrix result = new Matrix(_nrows, _ncols);
	foreach(row; 0 .. _nrows) {
	    foreach(col; 0 .. _ncols) {
		static if ( op == "+" ) {
		    result._data[row][col] = _data[row][col] + rhs._data[row][col];
		} else if ( op == "-" ) {
		    result._data[row][col] = _data[row][col] - rhs._data[row][col];
		}
	    }
	}
	return result;
    }

    Matrix opBinary(string op)(double rhs)
	if ( op == "+" || op == "-" || op == "*" || op == "/" )
    {
	Matrix result = new Matrix(_nrows, _ncols);
	foreach(row; 0 .. _nrows) {
	    foreach(col; 0 .. _ncols) {
		static if ( op == "+" ) {
		    result._data[row][col] = _data[row][col] + rhs;
		} else if ( op == "-" ) {
		    result._data[row][col] = _data[row][col] - rhs;
		} else if ( op == "*" ) {
		    result._data[row][col] = _data[row][col] * rhs;
		} else if ( op == "/" ) {
		    result._data[row][col] = _data[row][col] / rhs;
		}
	    }
	}
	return result;
    }

    Matrix opBinaryRight(string op)(double lhs)
	if ( op == "+" || op == "-" || op == "*" )
    {
	Matrix result = new Matrix(_nrows, _ncols);
	foreach(row; 0 .. _nrows) {
	    foreach(col; 0 .. _ncols) {
		static if ( op == "+" ) {
		    result._data[row][col] = lhs + _data[row][col];
		} else if ( op == "-" ) {
		    result._data[row][col] = lhs - _data[row][col];
		} else if ( op == "*" ) {
		    result._data[row][col] = lhs * _data[row][col];
		}
	    }
	}
	return result;
    }

    override string toString() {
	string s = "Matrix[";
	foreach(row; 0 .. _nrows) {
	    s ~= "[";
	    foreach(col; 0 .. _ncols) {
		s ~= to!string(_data[row][col]);
		if ( col < _ncols-1 ) s ~= ",";
	    }
	    s ~= "]";
	    if ( row < _nrows-1 ) s ~= ",";
	}
	s ~= "]";
	return s;
    }

    void swapRows(size_t i1, size_t i2) {
	swap(_data[i1], _data[i2]);
    }

    double[] getColumn(size_t col) {
	double[] my_column;
	my_column.length = nrows;
	foreach(row; 0 .. nrows) my_column[row] = _data[row][col];
	return my_column;
    }

    double[] getRow(size_t row) {
	return _data[row].dup;
    }

    // Maybe there is a way of using built-in slices here.
    Matrix sliceDup(size_t row0, size_t row1, size_t col0, size_t col1)
    {
	enforce((row0 < row1) && (row1 <= this._nrows) &&
		(col0 < col1) && (col1 <= this._ncols),
		new Error("invalid subrange"));
	Matrix sub_matrix = new Matrix(row1-row0, col1-col0);
	foreach(row; row0 .. row1) {
	    foreach(col; col0 .. col1) {
		sub_matrix[row-row0, col-col0] = this._data[row][col];
	    }
	}
	return sub_matrix;
    }

    void sliceAssign(double c, size_t row0, size_t row1, size_t col0, size_t col1)
    {
	enforce((row0 < row1) && (row1 <= this._nrows) &&
		(col0 < col1) && (col1 <= this._ncols),
		new Error("invalid subrange"));
	foreach(row; row0 .. row1) {
	    foreach(col; col0 .. col1) {
		this._data[row][col] = c;
	    }
	}
    }

} // end class Matrix


bool approxEqualMatrix(in Matrix a, in Matrix b)
{
    if ( a.nrows != b.nrows ) return false;
    if ( a.ncols != b.ncols ) return false;
    bool is_equal = true;
    foreach(row; 0 .. a.nrows) {
	foreach(col; 0 .. a.ncols) {
	    if ( !approxEqual(a[row,col], b[row,col]) ) {
		is_equal = false;
		break;
	    }
	}
    }
    return is_equal;
}

Matrix zeros(size_t rows, size_t cols)
{
    Matrix my_matrix = new Matrix(rows, cols);
    foreach(row; 0 .. rows) {
	foreach(col; 0 .. cols) {
	    my_matrix[row,col] = 0.0;
	}
    }
    return my_matrix;
}

Matrix eye(size_t n)
{
    Matrix ident_matrix = new Matrix(n);
    foreach(row; 0 .. n) {
	foreach(col; 0 .. n) {
	    ident_matrix[row,col] = (row == col) ? 1.0 : 0.0;
	}
    }
    return ident_matrix;
}

Matrix transpose(in Matrix other)
{
    Matrix my_matrix = new Matrix(other.ncols, other.nrows);
    foreach(row; 0 .. other.nrows) {
	foreach(col; 0 .. other.ncols) {
	    my_matrix[col,row] = other[row,col];
	}
    }
    return my_matrix;
}

Matrix hstack(in Matrix[] matrixList)
{
    bool consistent = true;
    size_t nrows = matrixList[0].nrows;
    size_t ncols = 0;
    foreach(mat; matrixList) {
	ncols += mat.ncols;
	if ( nrows != mat.nrows ) consistent = false;
    }
    if ( !consistent ) {
	throw new Error("Matrices need to have the same number of rows");
    }
    Matrix result = new Matrix(nrows, ncols);
    size_t colStart = 0;
    foreach(mat; matrixList) {
	foreach(row; 0 .. mat.nrows) {
	    foreach(col; 0 .. mat.ncols) {
		result[row, colStart+col] = mat[row, col];
	    }
	}
	colStart += mat.ncols;
    }
    return result;
}

Matrix vstack(in Matrix[] matrixList)
{
    bool consistent = true;
    size_t ncols = matrixList[0].ncols;
    size_t nrows = 0;
    foreach(mat; matrixList) {
	nrows += mat.nrows;
	if ( ncols != mat.ncols ) consistent = false;
    }
    if ( !consistent ) {
	throw new Error("Matrices need to have the same number of columns");
    }
    Matrix result = new Matrix(nrows, ncols);
    size_t rowStart = 0;
    foreach(mat; matrixList) {
	foreach(row; 0 .. mat.nrows) {
	    foreach(col; 0 .. mat.ncols) {
		result[rowStart+row, col] = mat[row, col];
	    }
	}
	rowStart += mat.nrows;
    }
    return result;
}

Matrix dot(in Matrix a, in Matrix b)
{
    if ( a.ncols != b.nrows ) {
	throw new Exception("incompatible matrices for dot product");
    }
    size_t nrows = a.nrows;
    size_t ncols = b.ncols;
    Matrix c = zeros(nrows, ncols);
    foreach(row; 0 .. nrows) {
	foreach(col; 0 .. ncols) {
	    foreach(i; 0 .. a.ncols) {
		c[row,col] += a[row,i] * b[i,col];
	    }
	}
    }
    return c;
}

unittest {
    Matrix a = eye(3);
    assert(approxEqualMatrix(a, new Matrix([[1,0,0],[0,1,0],[0,0,1]])),
	   "eye(), indentity matrix");
    Matrix b = zeros(2,3);
    Matrix c = transpose(b);
    b[1,2] = 99.0;
    c[0,0] = 1.0; c[1,1] = 1.0;
    assert(approxEqualMatrix(b, new Matrix([[0,0,0],[0,0,99]])),
	   "zeros and assign");
    assert(approxEqualMatrix(c, new Matrix([[1,0],[0,1],[0,0]])),
	   "transpose and assign");

    Matrix e = new Matrix([1.0, 2.0, 3.0]);
    assert(approxEqualMatrix(e, new Matrix([[1],[2],[3]])),
	   "column vector entered as array");
    Matrix e2 = new Matrix([1, 2, 3], "row");
    assert(approxEqualMatrix(e2, new Matrix([[1,2,3]])),
	   "row vector entered as array of int");

    Matrix f = new Matrix([[1.0,2.0,3.0],[4.0,5.0,6.0]]);
    Matrix g = dot(f,c);
    assert(approxEqualMatrix(g, new Matrix([[1,2],[4,5]])),
	   "dot product");
}


/**
 * Perform Gauss-Jordan elimination on an augmented matrix.
 * c = [A|b] such that the mutated matrix becomes [I|x]
 * where x is the solution vector(s) to A.x = b
 */
void gaussJordanElimination(ref Matrix c, double very_small_value=1.0e-16)
{
    if (c.ncols < c.nrows) {
	throw new Exception("too few columns supplied");
    }
    foreach(j; 0 .. c.nrows) {
	// Select pivot.
	size_t p = j;
	foreach(i; j+1 .. c.nrows) {
	    if ( abs(c[i,j]) > abs(c[p,j]) ) p = i;
	}
	if ( abs(c[p,j]) < very_small_value ) {
	    throw new Exception("matrix is essentially singular");
	}
	if ( p != j ) c.swapRows(p,j);
	// Scale row j to get unity on the diagonal.
	double cjj = c[j,j];
	foreach(col; 0 .. c.ncols) c[j,col] /= cjj;
	// Do the elimination to get zeros in all off diagonal values in column j.
	foreach(i; 0 .. c.nrows) {
	    if ( i == j ) continue;
	    double cij = c[i,j];
	    foreach(col; 0 .. c.ncols) c[i,col] -= cij * c[j,col]; 
	}
    } // end foreach j
} // end gaussJordanElimination()

unittest {
    Matrix A = new Matrix([[0.0,  2.0,  0.0,  1.0],
			   [2.0,  2.0,  3.0,  2.0],
			   [4.0, -3.0,  0.0,  1.0],
			   [6.0,  1.0, -6.0, -5.0]]);
    Matrix b = new Matrix([0.0, -2.0, -7.0, 6.0], "column");
    Matrix Ab = hstack([A,b]);
    Matrix Aonly = Ab.sliceDup(0, 4, 0, 4);
    Matrix bonly = Ab.sliceDup(0, 4, 4, 5);
    gaussJordanElimination(Ab);
    assert(approxEqualMatrix(Ab, new Matrix([[1,0,0,0,-0.5],[0,1.0,0,0,1],
					     [0,0,1,0,1.0/3],[0,0,0,1,-2.0]])),
	   "Gauss-Jordan elimination with partial pivoting");
    double[] x = Ab.getColumn(4);
    Matrix new_rhs = dot(Aonly, new Matrix(x));
    assert(approxEqualMatrix(new_rhs, bonly), "check rhs");
    Matrix residual = new_rhs - b;
    assert(approxEqualMatrix(residual, new Matrix([0,0,0,0])), "zero residual");
}

/**
 * First stage of a linear equation solver that uses LU decomposition 
 * followed by a separate backsubstitution.
 *
 * Params:
 *     c: incoming matrix, outgoing LU matrices combined 
 *        as described in Section 2.2 of Gerald and Wheatley
 *        section 2.2 Elimination methods.
 *     very_small_value: used in our test for singularity
 *
 * Returns:
 *     A list of row-permutation pairs (the indices of the swapped rows).
 *     Since we allow partial-pivoting by swapping rows,
 *     we need to keep a record of the row permutations
 *     to be passed into the solve function.
 */
size_t[2][] decomp(ref Matrix c, double very_small_value=1.0e-16)
{
    if (c.ncols != c.nrows) {
	throw new Exception("require a square matrix");
    }
    size_t[2][] permuteList;

    foreach(j; 0 .. c.nrows) {
	// Select pivot.
	size_t p = j;
	foreach(i; j+1 .. c.nrows) {
	    if ( abs(c[i,j]) > abs(c[p,j]) ) p = i;
	}
	if ( abs(c[p,j]) < very_small_value ) {
	    throw new Exception("matrix is essentially singular");
	}
	if ( p != j ) {
	    c.swapRows(p,j);
	    permuteList ~= [p,j];
	}
	// Do the elimination to get zeros in column j, below the diagonal.
	// Don't disturb the previous multipliers stored in columns to the left.
	foreach(i; j+1 .. c.nrows) {
	    double multiplier = c[i,j]/c[j,j];
	    foreach(col; j .. c.ncols) c[i,col] -= multiplier * c[j,col];
	    c[i,j] = multiplier; // Save in the newly-zeroed spot.
	}
    } // end foreach j

    return permuteList;
} // end decomp()

/**
 * Second stage of the linear equation solver.
 * 
 * Params:
 *     c:  decomposed matrix from the first stage.
 *         This matrix may be reused for any number of RHS vectors.
 *     rhs: on input, one or more column vectors for the right-hand side
 *         on return, these are the solution vectors.
 *     permuteList: list of row-permutation pairs from the first stage.
 */
void solve(in Matrix c, ref Matrix rhs, in size_t[2][] permuteList)
{
    size_t nrows = c.nrows;
    if (rhs.ncols < 1 || rhs.nrows != nrows) {
	throw new Exception("invalid right-hand side");
    }
    // Get the right-hand side rows into final order.
    foreach(rowPair; permuteList) {
	rhs.swapRows(rowPair[0], rowPair[1]);
    }
    // Forward elimination, using the stored multipliers.
    foreach(i; 1 .. nrows) {
	foreach(j; 0 .. i) {
	    double multiplier = c[i,j];
	    foreach(col; 0 .. rhs.ncols) rhs[i,col] -= multiplier * rhs[j,col];
	}
    }
    // Back substitution to obtain the solution vector(s).
    foreach(col; 0 .. rhs.ncols) {
	rhs[nrows-1, col] /= c[nrows-1,nrows-1];
	for (int i = to!int(nrows-2); i >= 0; --i) {
	    double my_sum = rhs[i,col];
	    foreach(j; i+1 .. nrows) my_sum -= c[i,j] * rhs[j,col];
	    rhs[i,col] = my_sum/c[i,i];
	}
    }
} // end solve()

Matrix inverse(in Matrix a)
{
    auto n = a.nrows;
    if ( n != a.ncols && n == 0 ) {
	throw new Exception("matrix should be square and not empty");
    }
    auto c = new Matrix(a);
    auto perm = decomp(c);
    auto x = eye(n);
    solve(c, x, perm);
    return x;
}

unittest {
    auto A = new Matrix([[0.0,  2.0,  0.0,  1.0],
			 [2.0,  2.0,  3.0,  2.0],
			 [4.0, -3.0,  0.0,  1.0],
			 [6.0,  1.0, -6.0, -5.0]]);
    auto b = new Matrix([0.0, -2.0, -7.0, 6.0], "column");
    auto c = new Matrix(A);
    auto perm = decomp(c);
    auto x = new Matrix(b);
    solve(c, x, perm);
    auto residual = b - dot(A,x);
    assert(approxEqualMatrix(residual, new Matrix([0,0,0,0])), "zero residual");

    auto y = inverse(A);
    assert(approxEqualMatrix(dot(A,y), eye(4)), "inverse calculation");
}

/**
 * Solve an over-constrained linear system of constraint equations c.x=rhs
 * in a least-squares sense.
 *
 * Params:
 *     c: Matrix of coefficients with nrows >= ncols
 *     rhs: right-hand-side vector of size nrows,ncols2
 *
 * Returns:
 *     x: solution matrix of size ncols, ncols2
 *
 * Notes:
 *     ncols2 may be larger than 1 so that several right-hand-sides
 *     can be solved in one sitting.
 */
Matrix lsqsolve(const Matrix c, const Matrix rhs)
{
    size_t m = c.ncols; // number of unknowns
    size_t N = c.nrows; // number of linear constraint equations
    if (N < m) {
	throw new Error(text("too few constraints N=", N, " m=", m));
    }
    // Prepare the normal equations A.x = b
    Matrix a = zeros(m, m);
    Matrix x = zeros(m, rhs.ncols);
    foreach (k; 0 .. m) {
	foreach (j; 0 .. m) {
	    foreach (i; 0 .. N) { a[k,j] += c[i,k]*c[i,j]; }
	}
	foreach (j; 0 .. rhs.ncols) {
	    foreach (i; 0 .. N) { x[k,j] += c[i,k]*rhs[i,j]; }
	}
    }
    // Solve using decomposition and solve for the usual mxm linear system.
    auto perm = decomp(a);
    solve(a, x, perm);
    return x;
} // end lsqsolve()

unittest {
    auto A = new Matrix([[0.0,  2.0,  0.0,  1.0],
			 [2.0,  2.0,  3.0,  2.0],
			 [4.0,  4.0,  6.0,  4.0],
			 [4.0, -3.0,  0.0,  1.0],
			 [4.0, -3.0,  0.0,  1.0],
			 [6.0,  1.0, -6.0, -5.0]]);
    auto b = new Matrix([[ 0.0],
			 [-2.0],
			 [-4.0],
			 [-7.0],
			 [-7.0],
			 [ 6.0]]);
    auto xx = lsqsolve(A, b);
    auto expected_xx = new Matrix([-0.5, 1, 1.0/3, -2], "column");
    assert(approxEqualMatrix(xx, expected_xx), "least-squares solve");
}
