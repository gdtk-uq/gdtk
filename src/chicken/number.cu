// number.cu
// Include file for chicken, to allow build for float or double numbers.
// PJ 2022-09-11

#ifndef NUMBER_INCLUDED
#define NUMBER_INCLUDED

# ifdef FLOAT_NUMBERS
typedef float number;
#else
typedef double number;
#endif

#endif
