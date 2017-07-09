# cuda_d
D bindings for CUDA

[CUDA®](https://developer.nvidia.com/cuda-zone) is a parallel computing platform
and programming model invented by NVIDIA. It enables dramatic increases in computing
performance by harnessing the power of the graphics processing unit (GPU). With
millions of CUDA-enabled GPUs sold to date, software developers, scientists
and researchers are using GPU-accelerated computing for broad-ranging applications.

`cuda_d` helps you to use CUDA APIs in a D eco-system.

The current version provides bindings to:

1. CUDA
2. cuBLAS
2. cuBLASXT
4. cuRand
5. CUDA profiler

# Installation

1. Install `cuda` drivers.

2. Add `cuda_d` to `dependencies` and specify the `libs` in dub file. e.g.

`dub.json`:

```json
{
  "name": "cuda_d_example",
  "dependencies": {
    "cuda_d": "~>0.1.0"
  },
  "libs": [ "cuda", "cublas", "cudart" , "curand"],
  "description": "CUDA Example"
}
```

# Examples

The examples can be found [here](https://github.com/prasunanand/cuda_d_examples).

# LICENSE

This software is distributed under the [BSD 3-Clause License](LICENSE).

Copyright © 2017, Prasun Anand
