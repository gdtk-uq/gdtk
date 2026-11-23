/**
 * @file cuda_to_hip.h
 * @brief CUDA-to-HIP compatibility header for ROCm builds
 *
 * On ROCm this header maps CUDA runtime symbols to their HIP equivalents.
 * On CUDA it is a pass-through to the CUDA runtime header.
 *
 * This port was authored with assistance from Claude, an AI assistant by
 * Anthropic.
 */
#pragma once

#if defined(HIP) || defined(__HIP_PLATFORM_AMD__)

#include <hip/hip_runtime.h>

// Runtime API - memory
#define cudaMalloc                  hipMalloc
#define cudaFree                    hipFree
#define cudaMemcpy                  hipMemcpy
#define cudaMemcpyHostToDevice      hipMemcpyHostToDevice
#define cudaMemcpyDeviceToHost      hipMemcpyDeviceToHost

// Error handling
#define cudaError_t                 hipError_t
#define cudaSuccess                 hipSuccess
#define cudaGetLastError            hipGetLastError
#define cudaGetErrorString          hipGetErrorString

// Device management
#define cudaGetDeviceCount          hipGetDeviceCount

#else  // CUDA path

#include <cuda_runtime.h>

#endif
