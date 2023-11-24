/**
 * bbla.d
 * Bare-bones linear algebra functions.
 *
 * Author: Peter J. and Rowan G.
 * Version: 2014-06-28, just enough to get our CFD code going.
 *          2015-05-30, added least-squares solver
 *          2018-04-21, added Crout's reduction for LU decomposition
 */

module nm.bbla;

import std.conv;
import std.algorithm;
import std.math;
import ntypes.complex;
import std.exception;

class Matrix(T) {
    size_t _nrows;
    size_t _ncols;
    T[] _data;

    this(size_t n) {
        _nrows = n;
        _ncols = n;
        _data.length = n*n;
    }

    this(size_t nrows, size_t ncols) {
        _nrows = nrows;
        _ncols = ncols;
        _data.length = nrows*ncols;
    }

    this(in Matrix other) {
        this(other._nrows, other._ncols);
        foreach(row; 0 .. _nrows)
            foreach(col; 0 .. _ncols)
                _data[index(row,col)] = other._data[index(row,col)];
    }

    this(in T[] vec, string orient="column")
    {
        if ( orient == "column" ) {
            this(vec.length, 1);
            foreach(row; 0 .. _nrows) _data[index(row,0)] = vec[row];
        } else {
            this(1, vec.length);
            foreach(col; 0 .. _ncols) _data[col] = vec[col];
        }
    }

    this(in float[] vec, string orient="column") {
        T[] my_vec;
        my_vec.length = vec.length;
        foreach(i; 0 .. my_vec.length) my_vec[i] = vec[i];
        this(my_vec, orient);
    }

    this(in int[] vec, string orient="column") {
        T[] my_vec;
        my_vec.length = vec.length;
        foreach(i; 0 .. my_vec.length) my_vec[i] = vec[i];
        this(my_vec, orient);
    }

    this(in T[][] other) {
        this(other.length, other[0].length);
        foreach(row; 0 .. _nrows)
            foreach(col; 0 .. _ncols)
                _data[index(row,col)] = other[row][col];
    }

    this(in float[][] other) {
        this(other.length, other[0].length);
        foreach(row; 0 .. _nrows)
            foreach(col; 0 .. _ncols)
                _data[index(row,col)] = other[row][col];
    }

    this(in int[][] other) {
        this(other.length, other[0].length);
        foreach(row; 0 .. _nrows)
            foreach(col; 0 .. _ncols)
                _data[index(row,col)] = other[row][col];
    }

    @property @nogc const size_t nrows() { return _nrows; }
    @property @nogc const size_t ncols() { return _ncols; }

    // Thinking that we get only one level of inlining,
    // let's do the index calculation explicitly for each
    // of the following index methods.

    @nogc
    size_t index(size_t row, size_t col) const
    in (row < _nrows, "row out of bounds")
    in (col < _ncols, "col out of bounds")
    {
        pragma(inline, true);
        return row*_ncols + col;
    }

    @nogc
    const T opIndex(size_t row, size_t col)
    in (row < _nrows, "row out of bounds")
    in (col < _ncols, "col out of bounds")
    {
        pragma(inline, true);
        return _data[row*_ncols + col];
    }

    @nogc
    ref T opIndexAssign(T c, size_t row, size_t col)
    in (row < _nrows, "row out of bounds")
    in (col < _ncols, "col out of bounds")
    {
        pragma(inline, true);
        size_t i = row*_ncols + col;
        _data[i] = c;
        return _data[i];
    }

    @nogc
    ref T opIndexOpAssign(string op)(T c, size_t row, size_t col)
        if ( op == "+" || op == "-" || op == "*" || op == "/" )
    in (row < _nrows, "row out of bounds")
    in (col < _ncols, "col out of bounds")
    {
        pragma(inline, true);
        size_t i = row*_ncols + col;
        static if ( op == "+" )
            _data[i] += c;
        else if ( op == "-" )
            _data[i] -= c;
        else if ( op == "*" )
            _data[i] *= c;
        else if ( op == "/" )
            _data[i] /= c;
        return _data[i];
    }

    // Element-by-element operations.
    Matrix!T opBinary(string op)(in Matrix rhs)
        if ( op == "+" || op == "-" || op == "*" || op == "/"  )
    {
        enforce(_nrows == rhs._nrows && _ncols == rhs._ncols,
                "incompatible matrices");
        Matrix!T result = new Matrix!T(_nrows, _ncols);
        foreach(row; 0 .. _nrows) {
            foreach(col; 0 .. _ncols) {
                size_t i = index(row,col);
                static if ( op == "+" ) {
                    result._data[i] = _data[i] + rhs._data[i];
                } else if ( op == "-" ) {
                    result._data[i] = _data[i] - rhs._data[i];
                } else if ( op == "*" ) {
                    result._data[i] = _data[i] * rhs._data[i];
                } else if ( op == "/" ) {
                    result._data[i] = _data[i] / rhs._data[i];
                }
            }
        }
        return result;
    }

    Matrix!T opBinary(string op)(double rhs)
        if ( op == "+" || op == "-" || op == "*" || op == "/" )
    {
        Matrix!T result = new Matrix!T(_nrows, _ncols);
        foreach(row; 0 .. _nrows) {
            foreach(col; 0 .. _ncols) {
                size_t i = index(row,col);
                static if ( op == "+" ) {
                    result._data[i] = _data[i] + rhs;
                } else if ( op == "-" ) {
                    result._data[i] = _data[i] - rhs;
                } else if ( op == "*" ) {
                    result._data[i] = _data[i] * rhs;
                } else if ( op == "/" ) {
                    result._data[i] = _data[i] / rhs;
                }
            }
        }
        return result;
    }

    version(complex_numbers) {
        Matrix!T opBinary(string op)(Complex!double rhs)
            if ( op == "+" || op == "-" || op == "*" || op == "/" )
        {
            Matrix!T result = new Matrix!T(_nrows, _ncols);
            foreach(row; 0 .. _nrows) {
                foreach(col; 0 .. _ncols) {
                    size_t i = index(row,col);
                    static if ( op == "+" ) {
                        result._data[i] = _data[i] + rhs;
                    } else if ( op == "-" ) {
                        result._data[i] = _data[i] - rhs;
                    } else if ( op == "*" ) {
                        result._data[i] = _data[i] * rhs;
                    } else if ( op == "/" ) {
                        result._data[i] = _data[i] / rhs;
                    }
                }
            }
            return result;
        }
    }

    Matrix!T opBinaryRight(string op)(double lhs)
        if ( op == "+" || op == "-" || op == "*" )
    {
        Matrix!T result = new Matrix!T(_nrows, _ncols);
        foreach(row; 0 .. _nrows) {
            foreach(col; 0 .. _ncols) {
                size_t i = index(row,col);
                static if ( op == "+" ) {
                    result._data[i] = lhs + _data[i];
                } else if ( op == "-" ) {
                    result._data[i] = lhs - _data[i];
                } else if ( op == "*" ) {
                    result._data[i] = lhs * _data[i];
                }
            }
        }
        return result;
    }

    version(complex_numbers) {
        Matrix!T opBinaryRight(string op)(Complex!double lhs)
            if ( op == "+" || op == "-" || op == "*" )
        {
            Matrix!T result = new Matrix!T(_nrows, _ncols);
            foreach(row; 0 .. _nrows) {
                foreach(col; 0 .. _ncols) {
                    size_t i = index(row,col);
                    static if ( op == "+" ) {
                        result._data[i] = lhs + _data[i];
                    } else if ( op == "-" ) {
                        result._data[i] = lhs - _data[i];
                    } else if ( op == "*" ) {
                        result._data[i] = lhs * _data[i];
                    }
                }
            }
            return result;
        }
    }

    @nogc
    void scale(double s)
    {
        foreach(row; 0 .. _nrows)
            foreach(col; 0 .. _ncols)
                _data[index(row,col)] *= s;
    }

    @nogc
    void add(in Matrix rhs)
    {
        foreach(row; 0 .. _nrows)
            foreach(col; 0 .. _ncols) {
                size_t i = index(row,col);
                _data[i] += rhs._data[i];
            }
    }

    override string toString() {
        string s = "Matrix["; // [TODO] add string form of type T here
        foreach(row; 0 .. _nrows) {
            s ~= "[";
            foreach(col; 0 .. _ncols) {
                s ~= to!string(_data[index(row,col)]);
                if ( col < _ncols-1 ) s ~= ",";
            }
            s ~= "]";
            if ( row < _nrows-1 ) s ~= ",";
        }
        s ~= "]";
        return s;
    }

    @nogc
    void swapRows(size_t i1, size_t i2) {
        foreach(col; 0 .. _ncols)
            swap(_data[index(i1,col)], _data[index(i2,col)]);
    }

    T[] getColumn(size_t col) {
        T[] my_column;
        my_column.length = nrows;
        foreach(row; 0 .. nrows) { my_column[row] = _data[index(row,col)]; }
        return my_column;
    }

    void getColumn(T[] my_column, size_t col) {
        assert (my_column.length == nrows, "Incorrect column length.");
        foreach(row; 0 .. nrows) { my_column[row] = _data[index(row,col)]; }
        return;
    }

    void setColumn(T[] my_column, size_t col) {
        assert (my_column.length == nrows, "Incorrect column length.");
        foreach(row; 0 .. nrows) { _data[index(row,col)] = my_column[row]; }
        return;
    }

    T[] getRow(size_t row) {
        T[] my_row;
        my_row.length = ncols;
        foreach(col; 0 .. _ncols) my_row[col] = _data[index(row,col)];
        return my_row;
    }

    // Maybe there is a way of using built-in slices here.
    Matrix!T sliceDup(size_t row0, size_t row1, size_t col0, size_t col1)
    {
        enforce((row0 < row1) && (row1 <= this._nrows) &&
                (col0 < col1) && (col1 <= this._ncols),
                new Error("invalid subrange"));
        Matrix sub_matrix = new Matrix(row1-row0, col1-col0);
        foreach(row; row0 .. row1) {
            foreach(col; col0 .. col1) {
                sub_matrix[row-row0, col-col0] = this._data[index(row,col)];
            }
        }
        return sub_matrix;
    }

    void sliceAssign(T c, size_t row0, size_t row1, size_t col0, size_t col1)
    {
        enforce((row0 < row1) && (row1 <= this._nrows) &&
                (col0 < col1) && (col1 <= this._ncols),
                new Error("invalid subrange"));
        foreach(row; row0 .. row1) {
            foreach(col; col0 .. col1) {
                this._data[index(row,col)] = c;
            }
        }
    }

    @nogc
    void zeros()
    {
        foreach(row; 0 .. _nrows)
            foreach(col; 0 .. _ncols)
                _data[index(row,col)] = to!T(0.0);
    }

    @nogc
    void ones()
    {
        foreach(row; 0 .. _nrows)
            foreach(col; 0 .. _ncols)
                _data[index(row,col)] = to!T(1.0);
    }

    @nogc
    void eye()
    {
        foreach(row; 0 .. _nrows) {
            foreach(col; 0 .. _ncols) {
                _data[index(row,col)] = to!T(0.0);
            }
            _data[index(row,row)] = to!T(1.0);
        }
    }
} // end class Matrix


bool approxEqualMatrix(T)(in Matrix!T a, in Matrix!T b)
{
    if (a.nrows != b.nrows) return false;
    if (a.ncols != b.ncols) return false;
    bool is_equal = true;
    foreach(row; 0 .. a.nrows) {
        foreach(col; 0 .. a.ncols) {
            if ( !approxEqualNumbers(a[row,col], b[row,col]) ) {
                is_equal = false;
                break;
            }
        }
    }
    return is_equal;
}

Matrix!T zeros(T)(size_t rows, size_t cols)
{
    Matrix!T my_matrix = new Matrix!T(rows, cols);
    foreach(row; 0 .. rows) {
        foreach(col; 0 .. cols) {
            my_matrix[row,col] = to!T(0.0);
        }
    }
    return my_matrix;
}

Matrix!T ones(T)(size_t rows, size_t cols)
{
    Matrix!T my_matrix = new Matrix!T(rows, cols);
    foreach(row; 0 .. rows) {
        foreach(col; 0 .. cols) {
            my_matrix[row,col] = to!T(1.0);
        }
    }
    return my_matrix;
}

Matrix!T eye(T)(size_t n)
{
    Matrix!T ident_matrix = new Matrix!T(n);
    foreach(row; 0 .. n) {
        foreach(col; 0 .. n) {
            ident_matrix[row,col] = (row == col) ? to!T(1.0) : to!T(0.0);
        }
    }
    return ident_matrix;
}

Matrix!T transpose(T)(in Matrix!T other)
{
    Matrix!T my_matrix = new Matrix!T(other.ncols, other.nrows);
    foreach(row; 0 .. other.nrows) {
        foreach(col; 0 .. other.ncols) {
            my_matrix[col,row] = other[row,col];
        }
    }
    return my_matrix;
}

Matrix!T hstack(T)(in Matrix!T[] matrixList)
{
    bool consistent = true;
    size_t nrows = matrixList[0].nrows;
    size_t ncols = 0;
    foreach(mat; matrixList) {
        ncols += mat.ncols;
        if (nrows != mat.nrows) consistent = false;
    }
    if (!consistent) {
        throw new Error("Matrices need to have the same number of rows");
    }
    Matrix!T result = new Matrix!T(nrows, ncols);
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

Matrix!T vstack(T)(in Matrix!T[] matrixList)
{
    bool consistent = true;
    size_t ncols = matrixList[0].ncols;
    size_t nrows = 0;
    foreach(mat; matrixList) {
        nrows += mat.nrows;
        if (ncols != mat.ncols) consistent = false;
    }
    if (!consistent) {
        throw new Error("Matrices need to have the same number of columns");
    }
    Matrix!T result = new Matrix!T(nrows, ncols);
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

Matrix!T dot(T)(in Matrix!T a, in Matrix!T b)
{
    if (a.ncols != b.nrows) {
        throw new Exception("incompatible matrices for dot product");
    }
    size_t nrows = a.nrows;
    size_t ncols = b.ncols;
    Matrix!T c = zeros!T(nrows, ncols);
    foreach(row; 0 .. nrows) {
        foreach(col; 0 .. ncols) {
            foreach(i; 0 .. a.ncols) {
                c[row,col] += a[row,i] * b[i,col];
            }
        }
    }
    return c;
}

@nogc
void dot(T)(in Matrix!T a, in Matrix!T b, ref Matrix!T c)
in {
    assert(a.ncols == b.nrows);
    assert(a.nrows == c.nrows);
    assert(b.ncols == c.ncols);
}
do {
    size_t nrows = a.nrows;
    size_t ncols = b.ncols;
    c.zeros();
    foreach(row; 0 .. nrows) {
        foreach(col; 0 .. ncols) {
            foreach(i; 0 .. a.ncols) {
                c[row,col] += a[row,i] * b[i,col];
            }
        }
    }
}

/**
 * A Dot product that only considers the upper left portion
 * of the matrices.
 *
 * a: matrix on left of dot
 * aRow: work up to (but not including) this row number of a
 * aCol: work up to (but not including) this column number of a
 * b: matrix on right of dot
 * bCol: work up to (but not including) this column number of b
 * c: matrix where result is placed
 *
 * Notes:
 *   1. ensure that all matrices are pre-allocated in size.
 *   2. matrix c is only changed where the the new reuls is
 *      computed based on the supplied row and column range.
 */
@nogc
void dot(T)(in Matrix!T a, size_t aRow, size_t aCol,
            in Matrix!T b, size_t bCol,
            ref Matrix!T c)
in {
    assert(aRow <= a.nrows);
    assert(aCol <= a.ncols);
    assert(aCol <= b.nrows);
    assert(bCol <= b.ncols);
    assert(aRow <= c.nrows);
    assert(bCol <= c.ncols);
}
do {
    foreach (row; 0 .. aRow) {
        foreach (col; 0 .. bCol) {
            c[row,col] = to!T(0.0);
            foreach (i; 0 .. aCol) {
                c[row,col] += a[row,i] * b[i,col];
            }
        }
    }
}

@nogc
void dot(T)(in Matrix!T a, T[] b, T[] c)
in {
    assert(a.ncols == b.length);
    assert(a.nrows == c.length);
}
do {
    size_t nrows = a.nrows;
    size_t ncols = a.ncols;
    foreach (idx; 0..c.length) c[idx] = 0.0;
    foreach (row; 0 .. nrows) {
        foreach (col; 0 .. ncols) {
            c[row] += a[row,col] * b[col];
        }
    }
}

@nogc
void dot(T)(in Matrix!T a, size_t aRow, size_t aCol, T[] b, T[] c)
in {
    assert(aRow <= a.nrows);
    assert(aCol <= a.ncols);
    assert(aCol <= b.length);
    assert(aRow <= c.length);
}
do {
    foreach (idx; 0..c.length) c[idx] = 0.0;
    foreach (row; 0 .. aRow) {
        foreach (col; 0 .. aCol) {
            c[row] += a[row,col] * b[col];
        }
    }
}

/*
    Some additional dot product routines that use flattened bare arrays,
    which somtimes are used as matrices in performance sensitive areas
    of the code.
    @author: NNG (July 2023)
*/
@nogc
void dot(T)(in T[] a, size_t a_stride, size_t aRow, size_t aCol, T[] b, T[] c)
{
    foreach (idx; 0 .. c.length) c[idx] = 0.0;
    foreach (row; 0 .. aRow) {
        foreach (col; 0 .. aCol) {
            c[row] += a[row*a_stride+col] * b[col];
        }
    }
}

@nogc
void transpose_and_dot(T)(in T[] a, size_t a_stride, size_t aRow, size_t aCol, T[] b, T[] c)
{
/*
    This routine does an inplace tranpose and dot, like this:
    [a1,  a2,  a3,  a4,  a5,  a6,  a7][q]   [a1*q + b1*r + c1*z]
    [b1,  b2,  b3,  b4,  b5,  b6,  b7][r] = [a2*q + b2*r + c2*z]
    [c1,  c2,  c3,  c4,  c5,  c6,  c7][z]   [a3*q + b3*r + c3*z]
                                            [a4*q + b4*r + c4*z]
                                            [a5*q + b5*r + c5*z]
                                            [a6*q + b6*r + c6*z]
                                            [a7*q + b7*r + c7*z]
    @author: NNG (July 2023)
*/
    foreach (idx; 0 .. c.length) c[idx] = 0.0;
    foreach (row; 0 .. aRow) {
        foreach (col; 0 .. aCol) {
            c[col] += a[row*a_stride+col] * b[row];
        }
    }
}

@nogc
void copy(T)(in Matrix!T src, ref Matrix!T tgt)
in {
    assert(src.nrows == tgt.nrows);
    assert(src.ncols == tgt.ncols);
}
do {
    foreach (row; 0 .. src.nrows) {
        foreach (col; 0 .. src.ncols) {
            tgt[row,col] = src[row,col];
        }
    }
}


version(bbla_test) {
    import util.msg_service;
    import std.conv;
    import ntypes.complex;
    import nm.number;
    int test_basic_operations() {
        Matrix!number a = eye!number(3);
        assert(approxEqualMatrix!number(a, new Matrix!number([[1,0,0],[0,1,0],[0,0,1]])),
               failedUnitTest());
        Matrix!number b = zeros!number(2,3);
        Matrix!number c = transpose!number(b);
        b[1,2] = to!number(99.0);
        c[0,0] = to!number(1.0); c[1,1] = to!number(1.0);
        assert(approxEqualMatrix!number(b, new Matrix!number([[0,0,0],[0,0,99]])),
               failedUnitTest());
        assert(approxEqualMatrix!number(c, new Matrix!number([[1,0],[0,1],[0,0]])),
               failedUnitTest());

        Matrix!number e = new Matrix!number([1.0, 2.0, 3.0]);
        assert(approxEqualMatrix!number(e, new Matrix!number([[1],[2],[3]])),
               failedUnitTest());
        Matrix!number e2 = new Matrix!number([1, 2, 3], "row");
        assert(approxEqualMatrix!number(e2, new Matrix!number([[1,2,3]])),
               failedUnitTest());

        Matrix!number f = new Matrix!number([[1.0,2.0,3.0],[4.0,5.0,6.0]]);
        Matrix!number g = dot!number(f,c);
        assert(approxEqualMatrix!number(g, new Matrix!number([[1,2],[4,5]])),
               failedUnitTest());
        g.zeros();
        assert(approxEqualMatrix!number(g, new Matrix!number([[0,0],[0,0]])).
               failedUnitTest());
        dot!number(f, c, g);
        assert(approxEqualMatrix!number(g, new Matrix!number([[1,2],[4,5]])),
               failedUnitTest());

        g.eye();
        assert(approxEqualMatrix!number(g, new Matrix!number([[1,0],[0,1]])),
               failedUnitTest());

        return 0;
    }
}


/**
 * Perform Gauss-Jordan elimination on an augmented matrix.
 * c = [A|b] such that the mutated matrix becomes [I|x]
 * where x is the solution vector(s) to A.x = b
 */
@nogc
void gaussJordanElimination(T)(ref Matrix!T c, double very_small_value=1.0e-16)
{
    if (c.ncols < c.nrows) {
        throw new Exception("too few columns supplied");
    }
    foreach(j; 0 .. c.nrows) {
        // Select pivot.
        size_t p = j;
        foreach(i; j+1 .. c.nrows) {
            if (abs(c[i,j]) > abs(c[p,j])) p = i;
        }
        if (abs(c[p,j]) < very_small_value) {
            throw new Exception("matrix is essentially singular");
        }
        if (p != j) { c.swapRows(p,j); }
        // Scale row j to get unity on the diagonal.
        T cjj = c[j,j];
        foreach (col; 0 .. c.ncols) { c[j,col] /= cjj; }
        // Do the elimination to get zeros in all off diagonal values in column j.
        foreach (i; 0 .. c.nrows) {
            if (i == j) continue;
            T cij = c[i,j];
            foreach (col; 0 .. c.ncols) { c[i,col] -= cij * c[j,col]; }
        }
    } // end foreach j
} // end gaussJordanElimination()

version(bbla_test) {
    int test_elimination() {
        Matrix!number A = new Matrix!number([[0.0,  2.0,  0.0,  1.0],
                                             [2.0,  2.0,  3.0,  2.0],
                                             [4.0, -3.0,  0.0,  1.0],
                                             [6.0,  1.0, -6.0, -5.0]]);
        Matrix!number b = new Matrix!number([0.0, -2.0, -7.0, 6.0], "column");
        Matrix!number Ab = hstack!number([A,b]);
        Matrix!number Aonly = Ab.sliceDup(0, 4, 0, 4);
        Matrix!number bonly = Ab.sliceDup(0, 4, 4, 5);
        gaussJordanElimination!number(Ab);
        assert(approxEqualMatrix!number(Ab, new Matrix!number([[1,0,0,0,-0.5],[0,1.0,0,0,1],
                                                               [0,0,1,0,1.0/3],[0,0,0,1,-2.0]])),
               failedUnitTest());
        number[] x = Ab.getColumn(4);
        Matrix!number new_rhs = dot!number(Aonly, new Matrix!number(x));
        assert(approxEqualMatrix!number(new_rhs, bonly), failedUnitTest());
        Matrix!number residual = new_rhs - b;
        assert(approxEqualMatrix!number(residual, new Matrix!number([0,0,0,0])), failedUnitTest());

        return 0;
    }
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
size_t[2][] decomp(T)(ref Matrix!T c, double very_small_value=1.0e-16)
{
    if (c.ncols != c.nrows) {
        throw new Exception("require a square matrix");
    }
    size_t[2][] permuteList;

    foreach (j; 0 .. c.nrows) {
        // Select pivot.
        size_t p = j;
        foreach (i; j+1 .. c.nrows) {
            if (abs(c[i,j]) > abs(c[p,j])) { p = i; }
        }
        if (abs(c[p,j]) < very_small_value) {
            throw new Exception("matrix is essentially singular");
        }
        if (p != j) {
            c.swapRows(p,j);
            permuteList ~= [p,j];
        }
        // Do the elimination to get zeros in column j, below the diagonal.
        // Don't disturb the previous multipliers stored in columns to the left.
        foreach (i; j+1 .. c.nrows) {
            T multiplier = c[i,j]/c[j,j];
            foreach (col; j .. c.ncols) { c[i,col] -= multiplier * c[j,col]; }
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
void solve(T)(in Matrix!T c, ref Matrix!T rhs, in size_t[2][] permuteList)
{
    size_t nrows = c.nrows;
    if (rhs.ncols < 1 || rhs.nrows != nrows) {
        throw new Exception("invalid right-hand side");
    }
    // Get the right-hand side rows into final order.
    foreach (rowPair; permuteList) {
        rhs.swapRows(rowPair[0], rowPair[1]);
    }
    // Forward elimination, using the stored multipliers.
    foreach (i; 1 .. nrows) {
        foreach (j; 0 .. i) {
            T multiplier = c[i,j];
            foreach (col; 0 .. rhs.ncols) { rhs[i,col] -= multiplier * rhs[j,col]; }
        }
    }
    // Back substitution to obtain the solution vector(s).
    foreach (col; 0 .. rhs.ncols) {
        rhs[nrows-1, col] /= c[nrows-1,nrows-1];
        for (int i = to!int(nrows-2); i >= 0; --i) {
            T my_sum = rhs[i,col];
            foreach (j; i+1 .. nrows) { my_sum -= c[i,j] * rhs[j,col]; }
            rhs[i,col] = my_sum/c[i,i];
        }
    }
} // end solve()

Matrix!T inverse(T)(in Matrix!T a)
{
    auto n = a.nrows;
    if (n != a.ncols && n == 0) {
        throw new Exception("matrix should be square and not empty");
    }
    auto c = new Matrix!T(a);
    auto perm = decomp!T(c);
    auto x = eye!T(n);
    solve!T(c, x, perm);
    return x;
}

void upperSolve(T)(in Matrix!T U, T[] b)
{
    int n = to!int(U.nrows);
    // Back subsitution
    b[n-1] /= U[n-1,n-1];
    for (int i = to!int(n-2); i >= 0; --i) {
        T sum = b[i];
        foreach (j; i+1 .. n) { sum -= U[i,j] * b[j]; }
        b[i] = sum/U[i,i];
    }
}

void upperSolve(T)(in Matrix!T U, int n, T[] b)
in {
    assert(n <= U.nrows);
    assert(n <= U.ncols);
    assert(n <= b.length);
}
do {
    // Back subsitution
    b[n-1] /= U[n-1,n-1];
    for (int i = to!int(n-2); i >= 0; --i) {
        T sum = b[i];
        foreach (j; i+1 .. n) { sum -= U[i,j] * b[j]; }
        b[i] = sum/U[i,i];
    }
}


version(bbla_test) {
    int test_decomp_and_inverse() {
        auto A = new Matrix!number([[0.0,  2.0,  0.0,  1.0],
                                    [2.0,  2.0,  3.0,  2.0],
                                    [4.0, -3.0,  0.0,  1.0],
                                    [6.0,  1.0, -6.0, -5.0]]);
        auto b = new Matrix!number([0.0, -2.0, -7.0, 6.0], "column");
        auto c = new Matrix!number(A);
        auto perm = decomp!number(c);
        auto x = new Matrix!number(b);
        solve!number(c, x, perm);
        auto residual = b - dot!number(A,x);
        assert(approxEqualMatrix!number(residual, new Matrix!number([0,0,0,0])),
               failedUnitTest());

        auto y = inverse!number(A);
        assert(approxEqualMatrix!number(dot!number(A,y), eye!number(4)), failedUnitTest());

        return 0;
    }
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
Matrix!T lsqsolve(T)(const Matrix!T c, const Matrix!T rhs)
{
    size_t m = c.ncols; // number of unknowns
    size_t N = c.nrows; // number of linear constraint equations
    if (N < m) {
        throw new Error(text("too few constraints N=", N, " m=", m));
    }
    // Prepare the normal equations A.x = b
    Matrix!T a = zeros!T(m, m);
    Matrix!T x = zeros!T(m, rhs.ncols);
    foreach (k; 0 .. m) {
        foreach (j; 0 .. m) {
            foreach (i; 0 .. N) { a[k,j] += c[i,k]*c[i,j]; }
        }
        foreach (j; 0 .. rhs.ncols) {
            foreach (i; 0 .. N) { x[k,j] += c[i,k]*rhs[i,j]; }
        }
    }
    // Solve using decomposition and solve for the usual mxm linear system.
    auto perm = decomp!T(a);
    solve!T(a, x, perm);
    return x;
} // end lsqsolve()()

version(bbla_test) {
    int test_lsqsolve() {
        auto A = new Matrix!number([[0.0,  2.0,  0.0,  1.0],
                                    [2.0,  2.0,  3.0,  2.0],
                                    [4.0,  4.0,  6.0,  4.0],
                                    [4.0, -3.0,  0.0,  1.0],
                                    [4.0, -3.0,  0.0,  1.0],
                                    [6.0,  1.0, -6.0, -5.0]]);
        auto b = new Matrix!number([[ 0.0],
                                    [-2.0],
                                    [-4.0],
                                    [-7.0],
                                    [-7.0],
                                    [ 6.0]]);
        auto xx = lsqsolve!number(A, b);
        auto expected_xx = new Matrix!number([-0.5, 1, 1.0/3, -2], "column");
        assert(approxEqualMatrix!number(xx, expected_xx), failedUnitTest());
        return 0;
    }
}

/**
 * Use Crout's reduction to form an LU decomposition.
 *
 * This code is adapted from Section 2.3 in Press et al. (2007).
 * The code in Numerical Recipes, 3rd edition is written for
 * C++ and makes use of classes. The implementation does not use
 * classes, just functions.
 *
 * Press et al. (2007)
 * Numerical Recipes: The Art of Scientific Computing, 3rd edition
 * Cambridge University Press, Cambridge, UK
 */

void LUDecomp(T)(Matrix!T a, ref int[] pivot, double verySmallValue=1.0e-40)
in {
    assert(a.nrows == a.ncols, "require a square matrix");
    assert(pivot.length == a.nrows, "pivot array wrongly sized");
}
do {
    int iMax;
    T largest, tmp;
    T[] vv;

    int n = to!int(a.nrows);
    vv.length = n;

    foreach (i; 0 .. n) {
        largest = 0.0;
        foreach (j; 0 .. n) {
            tmp = fabs(a[i,j]);
            if (tmp > largest) { largest = tmp; }
        }
        if (largest == 0.0)
            throw new Exception("Singular matrix in LUDecompCrout");
        vv[i] = 1.0/largest; // saves the scaling
    }

    foreach (k; 0 .. n) {
        largest = 0.0;
        foreach (i; k .. n) {
            tmp = vv[i]*fabs(a[i,k]);
            if (tmp > largest) {
                largest = tmp;
                iMax = i;
            }
        }
        if (k != iMax) { // if we need to interchange rows
            foreach (j; 0 .. n) {
                tmp = a[iMax,j];
                a[iMax,j] = a[k,j];
                a[k,j] = tmp;
            }
            vv[iMax] = vv[k];
        }
        pivot[k] = iMax;
        if (a[k,k] == 0.0) { a[k,k] = to!T(verySmallValue); }
        foreach (i; k+1 .. n) {
            a[i,k] /= a[k,k];
            tmp = a[i,k];
            foreach (j; k+1 .. n) {
                a[i,j] -= tmp*a[k,j];
            }
        }
    }
}

void LUSolve(T)(Matrix!T a, int[] pivot, T[] b, ref T[] x)
in {
    assert(a.nrows == a.ncols, "require a square matrix");
    assert(pivot.length == a.nrows, "pivot array wrongly sized");
    assert(b.length == a.nrows, "b array wrongly sized");
    assert(x.length == a.nrows, "x array wrongly sized");
}
do {
    int ii, ip, n;
    ii = 0;
    T sum;
    n = to!int(a.nrows);

    foreach (idx; 0..x.length) { x[idx] = b[idx]; }

    foreach (i; 0 .. n) {
        ip = pivot[i];
        sum = x[ip];
        x[ip] = x[i];
        if (ii != 0) {
            foreach (j; ii-1 .. i) {
                sum -= a[i,j]*x[j];
            }
        } else if (sum != 0.0) {
            ii = i + 1;
        }
        x[i] = sum;
    }

    for (int i = n-1; i >= 0; i--) {
        sum = x[i];
        foreach (j; i+1 .. n) { sum -= a[i,j]*x[j]; }
        x[i] = sum / a[i,i];
    }
}

version(bbla_test) {
    int test_Crout_decomp_and_solve() {
        auto A = new Matrix!number([[0.0,  2.0,  0.0,  1.0],
                                    [2.0,  2.0,  3.0,  2.0],
                                    [4.0, -3.0,  0.0,  1.0],
                                    [6.0,  1.0, -6.0, -5.0]]);
        number[] b = [to!number(0.0), to!number(-2.0), to!number(-7.0), to!number(6.0)];
        int[] pivot;
        pivot.length = 4;
        number[] x;
        x.length = 4;

        LUDecomp!number(A, pivot);
        LUSolve!number(A, pivot, b, x);

        // Reset A since it was converted to LU format.
        A = new Matrix!number([[0.0,  2.0,  0.0,  1.0],
                               [2.0,  2.0,  3.0,  2.0],
                               [4.0, -3.0,  0.0,  1.0],
                               [6.0,  1.0, -6.0, -5.0]]);

        number[] Ax;
        Ax.length = 4;
        dot(A, x, Ax);

        assert(approxEqualNumbers(b, Ax), failedUnitTest());

        return 0;
    }

    int main() {
        if ( test_basic_operations() != 0 )
            return 1;
        if ( test_elimination() != 0 )
            return 1;
        if ( test_decomp_and_inverse() != 0 )
            return 1;
        if ( test_lsqsolve() != 0 )
            return 1;
        if ( test_Crout_decomp_and_solve() != 0 )
            return 1;

        return 0;
    }


}
