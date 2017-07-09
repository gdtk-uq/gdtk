/*
  D bindings for CUDA.
  Authors:    Prasun Anand
  Copyright:  Copyright (c) 2017, Prasun Anand. All rights reserved.
  License:    BSD 3-Clause License
*/

module cuda_d.cuda_profiler_api;

import cuda_d.cuda;

extern (C):

cudaError_t cudaProfilerInitialize (
    const(char)* configFile,
    const(char)* outputFile,
    cudaOutputMode_t outputMode);

cudaError_t cudaProfilerStart ();

cudaError_t cudaProfilerStop ();

