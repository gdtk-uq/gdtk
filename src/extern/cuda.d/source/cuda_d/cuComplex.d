/*
  D bindings for CUDA.
  Authors:    Prasun Anand
  Copyright:  Copyright (c) 2017, Prasun Anand. All rights reserved.
  License:    BSD 3-Clause License
*/

module cuda_d.cuComplex;

import cuda_d.vector_types;

extern (C):

alias cuFloatComplex = float2_;

float cuCrealf (cuFloatComplex x);

float cuCimagf (cuFloatComplex x);

cuFloatComplex make_cuFloatComplex (float r, float i);

cuFloatComplex cuConjf (cuFloatComplex x);
cuFloatComplex cuCaddf (cuFloatComplex x, cuFloatComplex y);

cuFloatComplex cuCsubf (cuFloatComplex x, cuFloatComplex y);

cuFloatComplex cuCmulf (cuFloatComplex x, cuFloatComplex y);

cuFloatComplex cuCdivf (cuFloatComplex x, cuFloatComplex y);

float cuCabsf (cuFloatComplex x);

alias cuDoubleComplex = double2_;

double cuCreal (cuDoubleComplex x);

double cuCimag (cuDoubleComplex x);

cuDoubleComplex make_cuDoubleComplex (double r, double i);

cuDoubleComplex cuConj (cuDoubleComplex x);

cuDoubleComplex cuCadd (cuDoubleComplex x, cuDoubleComplex y);

cuDoubleComplex cuCsub (cuDoubleComplex x, cuDoubleComplex y);

cuDoubleComplex cuCmul (cuDoubleComplex x, cuDoubleComplex y);

cuDoubleComplex cuCdiv (cuDoubleComplex x, cuDoubleComplex y);

double cuCabs (cuDoubleComplex x);

alias cuComplex = float2_;

cuComplex make_cuComplex (float x, float y);

cuDoubleComplex cuComplexFloatToDouble (cuFloatComplex c);

cuFloatComplex cuComplexDoubleToFloat (cuDoubleComplex c);
