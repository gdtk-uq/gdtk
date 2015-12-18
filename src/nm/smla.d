/**
 * A minimal sparse matrix linear algebra module.
 *
 * Author: Rowan J. Gollan
 * Date: 2015-12-18
 *
 */

module nm.smla;

import std.stdio;
import std.string;
import std.conv;
import std.math;

/**
 * This class is used to store and work with a sparse matrix.
 * The data is stored in Compressed Sparse Row format.
 * Usually a class is used to hide details of implementation
 * (in this case, storage). However, here we will expose those
 * items as public method data since the algorithms to manipulate
 * sparse matrices rely on the representation of the data.
 *
 * Reference:
 * Saad (2003)
 * Iterative Methods for Sparse Linear Systems, 2nd ed.
 * SIAM, Philaldelphia
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

    const double opIndex(size_t row, size_t col)
    {
	// We need to search the given row to see if an entry
	// is present in the given column.
	foreach ( j; ia[row] .. ia[row+1] ) {
	    if ( ja[j] == col )
		return aa[j];
	}
	// else search failed, so we have a 0.0 entry
	return 0.0;
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
{
    assert(a.ia.length-1 == b.length);
    assert(b.length == c.length);
    size_t k0, k1;
    foreach ( i; 0 .. a.ia.length-1 ) {
	k0 = a.ia[i];
	k1 = a.ia[i+1];
	c[i] = 0.0;
	foreach ( k; k0 .. k1 ) c[i] += a.aa[k] * b[a.ja[k]];
    }
}

version(smla_test) {
    import util.msg_service;
    int main() {
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
	double[] v = [1., 1., 1., 1., 1.];
	double[] c;
	c.length = v.length;
	multiply(a, v, c);
	double[] expected_c = [3., 12., 30., 21., 12.];
	foreach ( i; 0 .. c.length ) {
	    assert(approxEqual(c[i], expected_c[i]), failedUnitTest());
	}

	return 0;
    }
}
