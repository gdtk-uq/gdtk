/*
  D bindings for CUDA.
  Authors:    Prasun Anand
  Copyright:  Copyright (c) 2017, Prasun Anand. All rights reserved.
  License:    BSD 3-Clause License
*/

extern (C):

enum CUDA_XT_DESCRIPTOR_VERSION = 0x01000000; // This is added to CUDART_VERSION

enum cudaXtCopyType_t
{
    LIB_XT_COPY_HOST_TO_DEVICE = 0,
    LIB_XT_COPY_DEVICE_TO_HOST = 1,
    LIB_XT_COPY_DEVICE_TO_DEVICE = 2
}

alias cudaLibXtCopyType = cudaXtCopyType_t;

enum libFormat_t
{
    LIB_FORMAT_CUFFT = 0,
    LIB_FORMAT_UNDEFINED = 1
}

alias libFormat = libFormat_t;

enum MAX_CUDA_DESCRIPTOR_GPUS = 64;

struct cudaXtDesc_t
{
    int version_; //descriptor version
    int nGPUs; //number of GPUs
    int[MAX_CUDA_DESCRIPTOR_GPUS] GPUs; //array of device IDs
    void*[MAX_CUDA_DESCRIPTOR_GPUS] data; //array of pointers to data, one per GPU
    size_t[MAX_CUDA_DESCRIPTOR_GPUS] size; //array of data sizes, one per GPU
    void* cudaXtState; //opaque CUDA utility structure
}

alias cudaXtDesc = cudaXtDesc_t;

struct cudaLibXtDesc_t
{
    int version_; //descriptor version
    cudaXtDesc* descriptor; //multi-GPU memory descriptor
    libFormat library; //which library recognizes the format
    int subFormat; //library specific enumerator of sub formats
    void* libDescriptor; //library specific descriptor e.g. FFT transform plan object
}

alias cudaLibXtDesc = cudaLibXtDesc_t;

