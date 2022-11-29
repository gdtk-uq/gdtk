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
#include <array>

#include "number.cu"

using namespace std;


string toString(const number c[2][2])
{
    stringstream ss;
    ss << "[[" << c[0][0] << "," << c[0][1] << "],[" << c[1][0] << "," << c[1][1] << "]]";
    return ss.str();
}

string toString2(const number x[2])
{
    stringstream ss;
    ss << "[" << x[0] << "," << x[1] << "]";
    return ss.str();
}

string toString(const number c[3][3])
{
    stringstream ss;
    ss << "[[" << c[0][0] << "," << c[0][1] << "," << c[0][2] << "],["
       << c[1][0] << "," << c[1][1] << "," << c[1][2] << "],["
       << c[2][0] << "," << c[2][1] << "," << c[2][2] << "]]";
    return ss.str();
}

string toString3(const number x[3])
{
    stringstream ss;
    ss << "[" << x[0] << "," << x[1] << "," << x[2] << "]";
    return ss.str();
}

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
struct Matrix {
    array<array<number, ncols>, nrows> data;

    string toString() const
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

    bool approxEqual(const Matrix<nrows,ncols> other, double tol=0.001)
    {
        bool isEqual = true; // Assume, then test elements.
        for (int i=0; i < nrows; ++i) {
            for (int j=0; j < ncols; ++j) {
                if (abs(data[i][j] - other.data[i][j]) > tol) isEqual = false;
            }
        }
        return isEqual;
    }

    int getColumn(int j, Matrix<nrows,1>& x)
    {
        if (j >= ncols) return -1;
        for (int i=0; i < nrows; ++i) { x.data[i][0] = data[i][j]; }
        return 0;
    }

    int getColumn(int j, number x[])
    {
        if (j >= ncols) return -1;
        for (int i=0; i < nrows; ++i) { x[i] = data[i][j]; }
        return 0;
    }
}; // end Matrix

template <size_t nrows, size_t ncols>
ostream& operator<<(ostream& os, const Matrix<nrows,ncols> m)
{
    os << m.toString();
    return os;
}

template <size_t nrows1, size_t ncols1, size_t ncols2>
Matrix<nrows1, ncols2> dot(const Matrix<nrows1,ncols1>& a, const Matrix<ncols1,ncols2>& b)
{
    Matrix<nrows1,ncols2> c;
    for (int i=0; i < nrows1; ++i) {
        for (int j=0; j < ncols2; ++j) {
            number s = 0.0;
            for (int k=0; k < ncols1; ++k) { s += a.data[i][k]*b.data[k][j]; }
            c.data[i][j] = s;
        }
    }
    return c;
}

template <size_t nrows1, size_t ncols1, size_t ncols2>
Matrix<nrows1, ncols1+ncols2> hstack(const Matrix<nrows1,ncols1>& a, const Matrix<nrows1,ncols2>& b)
{
    Matrix<nrows1,ncols1+ncols2> c;
    for (int i=0; i < nrows1; ++i) {
        for (int j=0; j < ncols1; ++j) { c.data[i][j] = a.data[i][j]; }
        for (int j=0; j < ncols2; ++j) { c.data[i][ncols1+j] = b.data[i][j]; }
    }
    return c;
}

/**
 * Perform Gauss-Jordan elimination on an augmented matrix.
 * c = [A|b] such that the mutated matrix becomes [I|x]
 * where x is the solution vector(s) to A.x = b
 */
template<size_t nrows, size_t ncols>
int gaussJordanElimination(Matrix<nrows,ncols>& c, double very_small_value=1.0e-16)
{
    static_assert(ncols > nrows, "too few columns supplied");
    for(int j=0; j < nrows; ++j) {
        // Select pivot.
        size_t p = j;
        for(int i=j+1; i < nrows; ++i) {
            if (abs(c.data[i][j]) > abs(c.data[p][j])) p = i;
        }
        if (abs(c.data[p][j]) < very_small_value) {
            // Matrix is essentially singular.
            return -1;
        }
        if (p != j) { swap(c.data[p], c.data[j]); }
        // Scale row j to get unity on the diagonal.
        number cjj = c.data[j][j];
        for (int col=0; col < ncols; ++col) { c.data[j][col] /= cjj; }
        // Do the elimination to get zeros in all off diagonal values in column j.
        for (int i=0; i < nrows; ++i) {
            if (i == j) continue;
            number cij = c.data[i][j];
            for (int col=0; col < ncols; ++col) { c.data[i][col] -= cij * c.data[j][col]; }
        }
    } // end for j
    return 0; // success
} // end gaussJordanElimination()

#endif
