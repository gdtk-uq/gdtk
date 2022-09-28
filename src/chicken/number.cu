// number.cu
// Include file for chicken, to allow build for float or double numbers.
// PJ 2022-09-11

#ifndef NUMBER_INCLUDED
#define NUMBER_INCLUDED

#include <cmath>

# ifdef FLOAT_NUMBERS
typedef float number;
#else
typedef double number;
#endif

bool approxEquals(double a, double b, double e) {
    return fabs(a - b) < e;
}

#endif
