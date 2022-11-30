// rsla.cu
// Real simple linear algebra for chicken flow solver.
// A linear equation solver using the text-book notation.
//
// PJ 2022-10-17 adapted the direct inverse calculations from the Dlang variant.
//   Calculation of the flow gradients (for the viscous terms) requires
//   the solution of 3x3 linear equations.
//   Would like to have used C++ array class but nvcc did not like them for __device__ fns.
//
// PJ 2022-11-29 Add templated Matrix class and functions.
//

#ifndef RSLA_INCLUDED
#define RSLA_INCLUDED

#include <cmath>
#include <sstream>
#include <iostream>

#include "number.cu"

using namespace std;

namespace rsla {

    __host__ __device__
    void MVMult(const number c[2][2], const number x[2], number result[2])
    {
        result[0] = c[0][0]*x[0] + c[0][1]*x[1];
        result[1] = c[1][0]*x[0] + c[1][1]*x[1];
    }

    __host__ __device__
    void MVMult(const number c[3][3], const number x[3], number result[3])
    {
        result[0] = c[0][0]*x[0] + c[0][1]*x[1] + c[0][2]*x[2];
        result[1] = c[1][0]*x[0] + c[1][1]*x[1] + c[1][2]*x[2];
        result[2] = c[2][0]*x[0] + c[2][1]*x[1] + c[2][2]*x[2];
    }

    __host__ __device__
    int MInverse(const number c[2][2], number d[2][2], number very_small_value=1.0e-12)
    {
        number det = c[0][0]*c[1][1] - c[0][1]*c[1][0];
        if (fabs(det) <= very_small_value) return -1; // singular
        // compute inverse directly
        number one_over_det = 1.0/det;
        d[0][0] =  c[1][1]*one_over_det; d[0][1] = -c[0][1]*one_over_det;
        d[1][0] = -c[1][0]*one_over_det; d[1][1] =  c[0][0]*one_over_det;
        return 0;
    }

    __host__ __device__
    int MInverse(const number c[3][3], number cinv[3][3], number very_small_value=1.0e-12)
    {
        number det = c[0][0]*(c[1][1]*c[2][2] - c[1][2]*c[2][1])
            - c[0][1]*(c[1][0]*c[2][2] - c[1][2]*c[2][0])
            + c[0][2]*(c[1][0]*c[2][1] - c[1][1]*c[2][0]);
        if (abs(det) <= very_small_value) return -1; // singular
        // compute inverse directly
        number one_over_det = 1.0/det;
        cinv[0][0] = (c[1][1]*c[2][2] - c[1][2]*c[2][1])*one_over_det;
        cinv[0][1] = (c[0][2]*c[2][1] - c[0][1]*c[2][2])*one_over_det;
        cinv[0][2] = (c[0][1]*c[1][2] - c[0][2]*c[1][1])*one_over_det;
        cinv[1][0] = (c[1][2]*c[2][0] - c[1][0]*c[2][2])*one_over_det;
        cinv[1][1] = (c[0][0]*c[2][2] - c[0][2]*c[2][0])*one_over_det;
        cinv[1][2] = (c[0][2]*c[1][0] - c[0][0]*c[1][2])*one_over_det;
        cinv[2][0] = (c[1][0]*c[2][1] - c[1][1]*c[2][0])*one_over_det;
        cinv[2][1] = (c[0][1]*c[2][0] - c[0][0]*c[2][1])*one_over_det;
        cinv[2][2] = (c[0][0]*c[1][1] - c[0][1]*c[1][0])*one_over_det;
        return 0;
    }

    // A more general matrix size can be implemented via templates,
    // in which we specify the size at compile time.

    template <size_t nrows, size_t ncols>
    void set_all_zero(number a[nrows][ncols])
    {
        for (int i=0; i < nrows; ++i) {
            for (int j=0; j < ncols; ++j) { a[i][j] = 0.0; }
        }
    }

    template <size_t nrows, size_t ncols>
    string toString(const number data[nrows][ncols])
    {
        stringstream ss;
        ss << "[";
        for (int i=0; i < nrows; ++i) {
            ss << "[";
            for (int j=0; j < ncols; ++j) {
                ss << data[i][j] << ((j < ncols-1) ? "," : "]");
            }
            ss << ((i < nrows-1) ? "," : "]");
        }
        return ss.str();
    }

    template <size_t nrows>
    string toString(const number data[nrows])
    {
        stringstream ss;
        ss << "[";
        for (int i=0; i < nrows; ++i) {
            ss << data[i];
            ss << ((i < nrows-1) ? "," : "]");
        }
        return ss.str();
    }

    template <size_t nrows, size_t ncols>
    bool approxEqual(const number data[nrows][ncols], const number other[nrows][ncols],
                     double tol=0.001)
    {
        bool isEqual = true; // Assume, then test elements.
        for (int i=0; i < nrows; ++i) {
            for (int j=0; j < ncols; ++j) {
                if (abs(data[i][j] - other[i][j]) > tol) isEqual = false;
            }
        }
        return isEqual;
    }

    template <size_t nrows>
    bool approxEqual(const number data[nrows][1], const number other[nrows],
                     double tol=0.001)
    {
        bool isEqual = true; // Assume, then test elements.
        for (int i=0; i < nrows; ++i) {
            if (abs(data[i][0] - other[i]) > tol) isEqual = false;
        }
        return isEqual;
    }

    template <size_t nrows, size_t ncols>
    int getColumn(number data[nrows][ncols], int j, number x[nrows])
    {
        if (j >= ncols) return -1;
        for (int i=0; i < nrows; ++i) { x[i] = data[i][j]; }
        return 0;
    }

    template <size_t nrows, size_t ncols>
    int getColumn(number data[nrows][ncols], int j, number x[nrows][1])
    {
        if (j >= ncols) return -1;
        for (int i=0; i < nrows; ++i) { x[i][0] = data[i][j]; }
        return 0;
    }

    template <size_t nrows1, size_t ncols1, size_t ncols2>
    void dot(const number a[nrows1][ncols1], const number b[ncols1][ncols2],
             number c[nrows1][ncols2])
    {
        for (int i=0; i < nrows1; ++i) {
            for (int j=0; j < ncols2; ++j) {
                number s = 0.0;
                for (int k=0; k < ncols1; ++k) { s += a[i][k]*b[k][j]; }
                c[i][j] = s;
            }
        }
    }

    template <size_t nrows1, size_t ncols1>
    void dot(const number a[nrows1][ncols1], const number b[ncols1],
             number c[nrows1])
    {
        for (int i=0; i < nrows1; ++i) {
            number s = 0.0;
            for (int k=0; k < ncols1; ++k) { s += a[i][k]*b[k]; }
            c[i] = s;
        }
    }

    template <size_t nrows1, size_t ncols1, size_t ncols2>
    void hstack(const number a[nrows1][ncols1], const number b[nrows1][ncols2],
                number c[nrows1][ncols1+ncols2])
    {
        for (int i=0; i < nrows1; ++i) {
            for (int j=0; j < ncols1; ++j) { c[i][j] = a[i][j]; }
            for (int j=0; j < ncols2; ++j) { c[i][ncols1+j] = b[i][j]; }
        }
    }

    template <size_t nrows1, size_t ncols1>
    void hstack(const number a[nrows1][ncols1], const number b[nrows1],
                number c[nrows1][ncols1+1])
    {
        for (int i=0; i < nrows1; ++i) {
            for (int j=0; j < ncols1; ++j) { c[i][j] = a[i][j]; }
            c[i][ncols1] = b[i];
        }
    }

    /**
     * Perform Gauss-Jordan elimination on an augmented matrix.
     * c = [A|b] such that the mutated matrix becomes [I|x]
     * where x is the solution vector(s) to A.x = b
     */
    template<size_t nrows, size_t ncols>
    int gaussJordanElimination(number c[nrows][ncols], double very_small_value=1.0e-16)
    {
        static_assert(ncols > nrows, "Too few columns supplied.");
        for(int j=0; j < nrows; ++j) {
            // Select pivot.
            size_t p = j;
            for(int i=j+1; i < nrows; ++i) {
                if (abs(c[i][j]) > abs(c[p][j])) p = i;
            }
            if (abs(c[p][j]) < very_small_value) {
                // Matrix is essentially singular.
                return -1;
            }
            if (p != j) { swap(c[p], c[j]); }
            // Scale row j to get unity on the diagonal.
            number cjj = c[j][j];
            for (int col=0; col < ncols; ++col) { c[j][col] /= cjj; }
            // Do the elimination to get zeros in all off diagonal values in column j.
            for (int i=0; i < nrows; ++i) {
                if (i == j) continue;
                number cij = c[i][j];
                for (int col=0; col < ncols; ++col) { c[i][col] -= cij * c[j][col]; }
            }
        } // end for j
        return 0; // success
    } // end gaussJordanElimination()

} // end namespace rsla

#endif
