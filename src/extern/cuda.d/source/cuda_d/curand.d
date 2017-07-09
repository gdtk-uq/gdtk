/*
  D bindings for CUDA.
  Authors:    Prasun Anand
  Copyright:  Copyright (c) 2017, Prasun Anand. All rights reserved.
  License:    BSD 3-Clause License
*/

module cuda_d.curand;

import cuda_d.cublas_api;

extern (C):

enum curandStatus
{
    CURAND_STATUS_SUCCESS = 0, ///< No errors
    CURAND_STATUS_VERSION_MISMATCH = 100, ///< Header file and linked library version do not match
    CURAND_STATUS_NOT_INITIALIZED = 101, ///< Generator not initialized
    CURAND_STATUS_ALLOCATION_FAILED = 102, ///< Memory allocation failed
    CURAND_STATUS_TYPE_ERROR = 103, ///< Generator is wrong type
    CURAND_STATUS_OUT_OF_RANGE = 104, ///< Argument out of range
    CURAND_STATUS_LENGTH_NOT_MULTIPLE = 105, ///< Length requested is not a multple of dimension
    CURAND_STATUS_DOUBLE_PRECISION_REQUIRED = 106, ///< GPU does not have double precision required by MRG32k3a
    CURAND_STATUS_LAUNCH_FAILURE = 201, ///< Kernel launch failure
    CURAND_STATUS_PREEXISTING_FAILURE = 202, ///< Preexisting failure on library entry
    CURAND_STATUS_INITIALIZATION_FAILED = 203, ///< Initialization of CUDA failed
    CURAND_STATUS_ARCH_MISMATCH = 204, ///< Architecture mismatch, GPU does not support requested feature
    CURAND_STATUS_INTERNAL_ERROR = 999 ///< Internal library error
}

alias curandStatus_t = curandStatus;

enum curandRngType
{
    CURAND_RNG_TEST = 0,
    CURAND_RNG_PSEUDO_DEFAULT = 100, ///< Default pseudorandom generator
    CURAND_RNG_PSEUDO_XORWOW = 101, ///< XORWOW pseudorandom generator
    CURAND_RNG_PSEUDO_MRG32K3A = 121, ///< MRG32k3a pseudorandom generator
    CURAND_RNG_PSEUDO_MTGP32 = 141, ///< Mersenne Twister MTGP32 pseudorandom generator
    CURAND_RNG_PSEUDO_MT19937 = 142, ///< Mersenne Twister MT19937 pseudorandom generator
    CURAND_RNG_PSEUDO_PHILOX4_32_10 = 161, ///< PHILOX-4x32-10 pseudorandom generator
    CURAND_RNG_QUASI_DEFAULT = 200, ///< Default quasirandom generator
    CURAND_RNG_QUASI_SOBOL32 = 201, ///< Sobol32 quasirandom generator
    CURAND_RNG_QUASI_SCRAMBLED_SOBOL32 = 202, ///< Scrambled Sobol32 quasirandom generator
    CURAND_RNG_QUASI_SOBOL64 = 203, ///< Sobol64 quasirandom generator
    CURAND_RNG_QUASI_SCRAMBLED_SOBOL64 = 204 ///< Scrambled Sobol64 quasirandom generator
}

alias curandRngType_t = curandRngType;

enum curandOrdering
{
    CURAND_ORDERING_PSEUDO_BEST = 100, ///< Best ordering for pseudorandom results
    CURAND_ORDERING_PSEUDO_DEFAULT = 101, ///< Specific default 4096 thread sequence for pseudorandom results
    CURAND_ORDERING_PSEUDO_SEEDED = 102, ///< Specific seeding pattern for fast lower quality pseudorandom results
    CURAND_ORDERING_QUASI_DEFAULT = 201 ///< Specific n-dimensional ordering for quasirandom results
}

/*
 * CURAND ordering of results in memory
 */
/** \cond UNHIDE_TYPEDEFS */
alias curandOrdering_t = curandOrdering;

enum curandDirectionVectorSet
{
    CURAND_DIRECTION_VECTORS_32_JOEKUO6 = 101, ///< Specific set of 32-bit direction vectors generated from polynomials recommended by S. Joe and F. Y. Kuo, for up to 20,000 dimensions
    CURAND_SCRAMBLED_DIRECTION_VECTORS_32_JOEKUO6 = 102, ///< Specific set of 32-bit direction vectors generated from polynomials recommended by S. Joe and F. Y. Kuo, for up to 20,000 dimensions, and scrambled
    CURAND_DIRECTION_VECTORS_64_JOEKUO6 = 103, ///< Specific set of 64-bit direction vectors generated from polynomials recommended by S. Joe and F. Y. Kuo, for up to 20,000 dimensions
    CURAND_SCRAMBLED_DIRECTION_VECTORS_64_JOEKUO6 = 104 ///< Specific set of 64-bit direction vectors generated from polynomials recommended by S. Joe and F. Y. Kuo, for up to 20,000 dimensions, and scrambled
}

alias curandDirectionVectorSet_t = curandDirectionVectorSet;

alias curandDirectionVectors32_t = uint[32];

alias curandDirectionVectors64_t = ulong[64];

struct curandGenerator_st;

alias curandGenerator_t = curandGenerator_st*;

alias curandDistribution_st = double;
alias curandDistribution_t = double*;
struct curandDistributionShift_st;
alias curandDistributionShift_t = curandDistributionShift_st*;

struct curandDistributionM2Shift_st;
alias curandDistributionM2Shift_t = curandDistributionM2Shift_st*;
struct curandHistogramM2_st;
alias curandHistogramM2_t = curandHistogramM2_st*;
alias curandHistogramM2K_st = uint;
alias curandHistogramM2K_t = uint*;
alias curandHistogramM2V_st = double;
alias curandHistogramM2V_t = double*;

struct curandDiscreteDistribution_st;
alias curandDiscreteDistribution_t = curandDiscreteDistribution_st*;

enum curandMethod
{
    CURAND_CHOOSE_BEST = 0, // choose best depends on args
    CURAND_ITR = 1,
    CURAND_KNUTH = 2,
    CURAND_HITR = 3,
    CURAND_M1 = 4,
    CURAND_M2 = 5,
    CURAND_BINARY_SEARCH = 6,
    CURAND_DISCRETE_GAUSS = 7,
    CURAND_REJECTION = 8,
    CURAND_DEVICE_API = 9,
    CURAND_FAST_REJECTION = 10,
    CURAND_3RD = 11,
    CURAND_DEFINITION = 12,
    CURAND_POISSON = 13
}

alias curandMethod_t = curandMethod;

curandStatus_t curandCreateGenerator (
    curandGenerator_t* generator,
    curandRngType_t rng_type);

curandStatus_t curandCreateGeneratorHost (
    curandGenerator_t* generator,
    curandRngType_t rng_type);


curandStatus_t curandDestroyGenerator (curandGenerator_t generator);


curandStatus_t curandGetVersion (int* version_);

curandStatus_t curandSetStream (
    curandGenerator_t generator,
    cudaStream_t stream);

curandStatus_t curandSetPseudoRandomGeneratorSeed (
    curandGenerator_t generator,
    ulong seed);

curandStatus_t curandSetGeneratorOffset (
    curandGenerator_t generator,
    ulong offset);

curandStatus_t curandSetGeneratorOrdering (
    curandGenerator_t generator,
    curandOrdering_t order);

curandStatus_t curandSetQuasiRandomGeneratorDimensions (
    curandGenerator_t generator,
    uint num_dimensions);

curandStatus_t curandGenerate (
    curandGenerator_t generator,
    uint* outputPtr,
    size_t num);

curandStatus_t curandGenerateLongLong (
    curandGenerator_t generator,
    ulong* outputPtr,
    size_t num);

curandStatus_t curandGenerateUniform (
    curandGenerator_t generator,
    float* outputPtr,
    size_t num);

curandStatus_t curandGenerateUniformDouble (
    curandGenerator_t generator,
    double* outputPtr,
    size_t num);

curandStatus_t curandGenerateNormal (
    curandGenerator_t generator,
    float* outputPtr,
    size_t n,
    float mean,
    float stddev);

curandStatus_t curandGenerateNormalDouble (
    curandGenerator_t generator,
    double* outputPtr,
    size_t n,
    double mean,
    double stddev);

curandStatus_t curandGenerateLogNormal (
    curandGenerator_t generator,
    float* outputPtr,
    size_t n,
    float mean,
    float stddev);

curandStatus_t curandGenerateLogNormalDouble (
    curandGenerator_t generator,
    double* outputPtr,
    size_t n,
    double mean,
    double stddev);

curandStatus_t curandCreatePoissonDistribution (
    double lambda,
    curandDiscreteDistribution_t* discrete_distribution);

curandStatus_t curandDestroyDistribution (
    curandDiscreteDistribution_t discrete_distribution);

curandStatus_t curandGeneratePoisson (
    curandGenerator_t generator,
    uint* outputPtr,
    size_t n,
    double lambda);

curandStatus_t curandGeneratePoissonMethod (
    curandGenerator_t generator,
    uint* outputPtr,
    size_t n,
    double lambda,
    curandMethod_t method);

curandStatus_t curandGenerateSeeds (curandGenerator_t generator);

curandStatus_t curandGetDirectionVectors32 (
    curandDirectionVectors32_t** vectors,
    curandDirectionVectorSet_t set);

curandStatus_t curandGetScrambleConstants32 (uint** constants);

curandStatus_t curandGetDirectionVectors64 (
    curandDirectionVectors64_t** vectors,
    curandDirectionVectorSet_t set);

curandStatus_t curandGetScrambleConstants64 (ulong** constants);
