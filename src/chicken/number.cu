// number.cu
// Include file for chicken, to allow build for float or double numbers.
// PJ 2022-09-11

#ifndef NUMBER_INCLUDED
#define NUMBER_INCLUDED

#include <cmath>

# ifdef FLOAT_NUMBERS
typedef float number;
constexpr number zero = 0.0f;
constexpr number half = 0.5f;
constexpr number one_quarter = 0.25f;
constexpr number one_fifth = 0.2f;
constexpr number one_sixth = 1.0f/6.0f;
constexpr number four_sixth = 4.0f/6.0f;
constexpr number one_eighth = 0.125f;
constexpr number one_twelfth = 1.0f/12.0f;
constexpr number one = 1.0f;
constexpr number two = 2.0f;
constexpr number three = 3.0f;
constexpr number four = 4.0f;
constexpr number five = 5.0f;
constexpr number six = 6.0f;
constexpr number seven = 7.0f;
constexpr number eight = 8.0f;
#else
typedef double number;
constexpr number zero = 0.0;
constexpr number half = 0.5;
constexpr number one_quarter = 0.25;
constexpr number one_fifth = 0.2;
constexpr number one_sixth = 1.0/6.0;
constexpr number four_sixth = 4.0/6.0;
constexpr number one_eighth = 0.125;
constexpr number one_twelfth = 1.0/12.0;
constexpr number one = 1.0;
constexpr number two = 2.0;
constexpr number three = 3.0;
constexpr number four = 4.0;
constexpr number five = 5.0;
constexpr number six = 6.0;
constexpr number seven = 7.0;
constexpr number eight = 8.0;
#endif

bool approxEquals(double a, double b, double e) {
    return fabs(a - b) < e;
}

#endif
