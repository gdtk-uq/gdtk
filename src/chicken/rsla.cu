// rsla.cu
// Real simple linear algebra for chicken flow solver.
// A linear equation solver using the text-book notation.
//
// PJ 2022-10-17 adapted the direct inverse calculations from the Dlang variant.
//   Calculation of the flow gradients (for the viscous terms) requires
//   the solution of 3x3 linear equations.
//   Would like to have used C++ array class but nvcc did like them for __device__ fns.

#ifndef RSLA_INCLUDED
#define RSLA_INCLUDED

#include <cmath>
#include <sstream>

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

#endif
