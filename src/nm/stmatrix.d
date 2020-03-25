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
import nm.complex;
import nm.number;

class STMatrix(T) {
public:
    size_t nrows, ncols, nzentries;
    size_t[] i;
    size_t[] j;
    T[] val;

    this(size_t _nrows, size_t _ncols, size_t _nzentries)
    {
        nrows = _nrows;
        ncols = _ncols;
        nzentries = _nzentries;
        i.length = nzentries;
        j.length = nzentries;
        val.length = nzentries;
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
    // Read nrows, ncols, nzentries
    auto tks = line.split();
    auto nrows = to!size_t(tks[0]);
    auto ncols = to!size_t(tks[1]);
    auto nzentries = to!size_t(tks[2]);
    auto matrix = new STMatrix!T(nrows, ncols, nzentries);
    // Now read all entries.
    // The MatrixMarket format allows duplicate entries for (i,j) BUT this code is specifically built for files from
    // the SuiteSparse collection. These matrices do NOT have duplicates, so we won't handle that case here.
    foreach (k; 0 .. nzentries) {
        line = f.readln().strip();
        tks = line.split();
        matrix.i[k] = to!size_t(tks[0]);
        matrix.j[k] = to!size_t(tks[1]);
        matrix.val[k] = to!T(tks[2]);
    }
    f.close();
    return matrix;
}


version(stmatrix_test) {
    import util.msg_service;
    int main() {
        string fName = "test_data/b1_ss.mtx";
        auto matrix = readFromMatrixMarketFile!double(fName);
        assert(matrix.i[0] == 5, failedUnitTest());
        assert(matrix.j[0] == 1, failedUnitTest());
        assert(approxEqualNumbers(matrix.val[0], -0.03599942, 1.0e-7), failedUnitTest());
        return 0;
    }
}
