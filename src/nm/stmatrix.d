/**
 * A simple container for a sparse triplet matrix.
 *
 * Author: Rowan J. Gollan
 * Date: 2020-03-25
 *
 */

module nm.stmatrix;

import core.stdc.stdlib : exit;
import std.stdio;
import std.string;
import std.conv;
import std.typecons : Tuple, tuple;
import ntypes.complex;
import nm.number;

class STMatrix(T) {
public:
    size_t n_rows, n_cols, n_nzeros;
    T[Tuple!(size_t, size_t)] val;

    this(size_t _n_rows, size_t _n_cols, size_t _n_nzeros)
    {
        n_rows = _n_rows;
        n_cols = _n_cols;
        n_nzeros = _n_nzeros;
    }
}

STMatrix!T readFromMatrixMarketFile(T)(string fName)
{
    auto f = File(fName, "r");
    // Keep throwing away comment lines until we find information about the size and number of entries.
    string line;
    while (true) {
        line = f.readln().strip();
        if (line[0] != '%') break;
    }
    // Read numer of rows, colums and non-zero entries
    auto tks = line.split();
    auto n_rows = to!size_t(tks[0]);
    auto n_cols = to!size_t(tks[1]);
    auto n_nzeros = to!size_t(tks[2]);
    auto matrix = new STMatrix!T(n_rows, n_cols, n_nzeros);
    // Now read all entries.
    // The MatrixMarket format allows duplicate entries for (i,j) BUT this code is specifically built for files from
    // the SuiteSparse collection. These matrices do NOT have duplicates, so we won't handle that case here.
    foreach (k; 0 .. n_nzeros) {
        line = f.readln().strip();
        tks = line.split();
        auto i = to!size_t(tks[0]) - 1; // Change 1- to 0-offset
        auto j = to!size_t(tks[1]) - 1; // Change 1- to 0-offset
        auto v = to!T(tks[2]);
        matrix.val[tuple(i,j)] = v;
    }
    f.close();
    return matrix;
}


version(stmatrix_test) {
    import util.msg_service;
    int main() {
        string fName = "test_data/b1_ss.mtx";
        auto matrix = readFromMatrixMarketFile!double(fName);
        assert(approxEqualNumbers(matrix.val[tuple!(size_t,size_t)(4,0)], -0.03599942, 1.0e-7), failedUnitTest());
        assert(approxEqualNumbers(matrix.val[tuple!(size_t,size_t)(6,6)], 1.0, 1.0e-7), failedUnitTest());
        return 0;
    }
}
