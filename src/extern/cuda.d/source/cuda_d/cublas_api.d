/*
  D bindings for CUDA.
  Authors:    Prasun Anand
  Copyright:  Copyright (c) 2017, Prasun Anand. All rights reserved.
  License:    BSD 3-Clause License
*/

module cuda_d.cublas_api;

import cuda_d.cuComplex;

extern (C):

alias void* cudaStream_t;


struct half {
   ushort x;
}

struct half2 {
   uint x;
}

enum cublasStatus_t
{
    CUBLAS_STATUS_SUCCESS = 0,
    CUBLAS_STATUS_NOT_INITIALIZED = 1,
    CUBLAS_STATUS_ALLOC_FAILED = 3,
    CUBLAS_STATUS_INVALID_VALUE = 7,
    CUBLAS_STATUS_ARCH_MISMATCH = 8,
    CUBLAS_STATUS_MAPPING_ERROR = 11,
    CUBLAS_STATUS_EXECUTION_FAILED = 13,
    CUBLAS_STATUS_INTERNAL_ERROR = 14,
    CUBLAS_STATUS_NOT_SUPPORTED = 15,
    CUBLAS_STATUS_LICENSE_ERROR = 16
}

enum cublasFillMode_t
{
    CUBLAS_FILL_MODE_LOWER = 0,
    CUBLAS_FILL_MODE_UPPER = 1
}

enum cublasDiagType_t
{
    CUBLAS_DIAG_NON_UNIT = 0,
    CUBLAS_DIAG_UNIT = 1
}

enum cublasSideMode_t
{
    CUBLAS_SIDE_LEFT = 0,
    CUBLAS_SIDE_RIGHT = 1
}

enum cublasOperation_t
{
    CUBLAS_OP_N = 0,
    CUBLAS_OP_T = 1,
    CUBLAS_OP_C = 2
}

enum cublasPointerMode_t
{
    CUBLAS_POINTER_MODE_HOST = 0,
    CUBLAS_POINTER_MODE_DEVICE = 1
}

enum cublasAtomicsMode_t
{
    CUBLAS_ATOMICS_NOT_ALLOWED = 0,
    CUBLAS_ATOMICS_ALLOWED = 1
}

/* Used by cublasSgemmEx */
enum cublasDataType_t
{
    CUBLAS_DATA_FLOAT = 0,
    CUBLAS_DATA_DOUBLE = 1,
    CUBLAS_DATA_HALF = 2,
    CUBLAS_DATA_INT8 = 3
}

/* Opaque structure holding CUBLAS library context */
struct cublasContext;
alias cublasHandle_t = cublasContext*;

cublasStatus_t cublasCreate_v2 (cublasHandle_t* handle);
cublasStatus_t cublasDestroy_v2 (cublasHandle_t handle);
cublasStatus_t cublasGetVersion_v2 (cublasHandle_t handle, int* version_);
cublasStatus_t cublasSetStream_v2 (cublasHandle_t handle, cudaStream_t streamId);
cublasStatus_t cublasGetStream_v2 (cublasHandle_t handle, cudaStream_t* streamId);

cublasStatus_t cublasGetPointerMode_v2 (cublasHandle_t handle, cublasPointerMode_t* mode);
cublasStatus_t cublasSetPointerMode_v2 (cublasHandle_t handle, cublasPointerMode_t mode);

cublasStatus_t cublasGetAtomicsMode (cublasHandle_t handle, cublasAtomicsMode_t* mode);
cublasStatus_t cublasSetAtomicsMode (cublasHandle_t handle, cublasAtomicsMode_t mode);


cublasStatus_t cublasSetVector (
    int n,
    int elemSize,
    const(void)* x,
    int incx,
    void* devicePtr,
    int incy);

cublasStatus_t cublasGetVector (
    int n,
    int elemSize,
    const(void)* x,
    int incx,
    void* y,
    int incy);

cublasStatus_t cublasSetMatrix (
    int rows,
    int cols,
    int elemSize,
    const(void)* A,
    int lda,
    void* B,
    int ldb);

cublasStatus_t cublasGetMatrix (
    int rows,
    int cols,
    int elemSize,
    const(void)* A,
    int lda,
    void* B,
    int ldb);

cublasStatus_t cublasSetVectorAsync (
    int n,
    int elemSize,
    const(void)* hostPtr,
    int incx,
    void* devicePtr,
    int incy,
    cudaStream_t stream);

cublasStatus_t cublasGetVectorAsync (
    int n,
    int elemSize,
    const(void)* devicePtr,
    int incx,
    void* hostPtr,
    int incy,
    cudaStream_t stream);


cublasStatus_t cublasSetMatrixAsync (
    int rows,
    int cols,
    int elemSize,
    const(void)* A,
    int lda,
    void* B,
    int ldb,
    cudaStream_t stream);

cublasStatus_t cublasGetMatrixAsync (
    int rows,
    int cols,
    int elemSize,
    const(void)* A,
    int lda,
    void* B,
    int ldb,
    cudaStream_t stream);

void cublasXerbla (const(char)* srName, int info);
/* ---------------- CUBLAS BLAS1 functions ---------------- */
cublasStatus_t cublasSnrm2_v2 (
    cublasHandle_t handle,
    int n,
    const(float)* x,
    int incx,
    float* result); /* host or device pointer */

cublasStatus_t cublasDnrm2_v2 (
    cublasHandle_t handle,
    int n,
    const(double)* x,
    int incx,
    double* result); /* host or device pointer */

cublasStatus_t cublasScnrm2_v2 (
    cublasHandle_t handle,
    int n,
    const(cuComplex)* x,
    int incx,
    float* result); /* host or device pointer */

cublasStatus_t cublasDznrm2_v2 (
    cublasHandle_t handle,
    int n,
    const(cuDoubleComplex)* x,
    int incx,
    double* result); /* host or device pointer */

cublasStatus_t cublasSdot_v2 (
    cublasHandle_t handle,
    int n,
    const(float)* x,
    int incx,
    const(float)* y,
    int incy,
    float* result); /* host or device pointer */

cublasStatus_t cublasDdot_v2 (
    cublasHandle_t handle,
    int n,
    const(double)* x,
    int incx,
    const(double)* y,
    int incy,
    double* result); /* host or device pointer */

cublasStatus_t cublasCdotu_v2 (
    cublasHandle_t handle,
    int n,
    const(cuComplex)* x,
    int incx,
    const(cuComplex)* y,
    int incy,
    cuComplex* result); /* host or device pointer */

cublasStatus_t cublasCdotc_v2 (
    cublasHandle_t handle,
    int n,
    const(cuComplex)* x,
    int incx,
    const(cuComplex)* y,
    int incy,
    cuComplex* result); /* host or device pointer */

cublasStatus_t cublasZdotu_v2 (
    cublasHandle_t handle,
    int n,
    const(cuDoubleComplex)* x,
    int incx,
    const(cuDoubleComplex)* y,
    int incy,
    cuDoubleComplex* result); /* host or device pointer */

cublasStatus_t cublasZdotc_v2 (
    cublasHandle_t handle,
    int n,
    const(cuDoubleComplex)* x,
    int incx,
    const(cuDoubleComplex)* y,
    int incy,
    cuDoubleComplex* result); /* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasSscal_v2 (
    cublasHandle_t handle,
    int n,
    const(float)* alpha,
    float* x,
    int incx);

/* host or device pointer */
cublasStatus_t cublasDscal_v2 (
    cublasHandle_t handle,
    int n,
    const(double)* alpha,
    double* x,
    int incx);

/* host or device pointer */
cublasStatus_t cublasCscal_v2 (
    cublasHandle_t handle,
    int n,
    const(cuComplex)* alpha,
    cuComplex* x,
    int incx);

/* host or device pointer */
cublasStatus_t cublasCsscal_v2 (
    cublasHandle_t handle,
    int n,
    const(float)* alpha,
    cuComplex* x,
    int incx);

/* host or device pointer */
cublasStatus_t cublasZscal_v2 (
    cublasHandle_t handle,
    int n,
    const(cuDoubleComplex)* alpha,
    cuDoubleComplex* x,
    int incx);

/* host or device pointer */
cublasStatus_t cublasZdscal_v2 (
    cublasHandle_t handle,
    int n,
    const(double)* alpha,
    cuDoubleComplex* x,
    int incx);

/* host or device pointer */
cublasStatus_t cublasSaxpy_v2 (
    cublasHandle_t handle,
    int n,
    const(float)* alpha,
    const(float)* x,
    int incx,
    float* y,
    int incy);

/* host or device pointer */
cublasStatus_t cublasDaxpy_v2 (
    cublasHandle_t handle,
    int n,
    const(double)* alpha,
    const(double)* x,
    int incx,
    double* y,
    int incy);

/* host or device pointer */
cublasStatus_t cublasCaxpy_v2 (
    cublasHandle_t handle,
    int n,
    const(cuComplex)* alpha,
    const(cuComplex)* x,
    int incx,
    cuComplex* y,
    int incy);

/* host or device pointer */
cublasStatus_t cublasZaxpy_v2 (
    cublasHandle_t handle,
    int n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* x,
    int incx,
    cuDoubleComplex* y,
    int incy);

cublasStatus_t cublasScopy_v2 (
    cublasHandle_t handle,
    int n,
    const(float)* x,
    int incx,
    float* y,
    int incy);

cublasStatus_t cublasDcopy_v2 (
    cublasHandle_t handle,
    int n,
    const(double)* x,
    int incx,
    double* y,
    int incy);

cublasStatus_t cublasCcopy_v2 (
    cublasHandle_t handle,
    int n,
    const(cuComplex)* x,
    int incx,
    cuComplex* y,
    int incy);

cublasStatus_t cublasZcopy_v2 (
    cublasHandle_t handle,
    int n,
    const(cuDoubleComplex)* x,
    int incx,
    cuDoubleComplex* y,
    int incy);

cublasStatus_t cublasSswap_v2 (
    cublasHandle_t handle,
    int n,
    float* x,
    int incx,
    float* y,
    int incy);

cublasStatus_t cublasDswap_v2 (
    cublasHandle_t handle,
    int n,
    double* x,
    int incx,
    double* y,
    int incy);

cublasStatus_t cublasCswap_v2 (
    cublasHandle_t handle,
    int n,
    cuComplex* x,
    int incx,
    cuComplex* y,
    int incy);

cublasStatus_t cublasZswap_v2 (
    cublasHandle_t handle,
    int n,
    cuDoubleComplex* x,
    int incx,
    cuDoubleComplex* y,
    int incy);

cublasStatus_t cublasIsamax_v2 (
    cublasHandle_t handle,
    int n,
    const(float)* x,
    int incx,
    int* result); /* host or device pointer */

cublasStatus_t cublasIdamax_v2 (
    cublasHandle_t handle,
    int n,
    const(double)* x,
    int incx,
    int* result); /* host or device pointer */

cublasStatus_t cublasIcamax_v2 (
    cublasHandle_t handle,
    int n,
    const(cuComplex)* x,
    int incx,
    int* result); /* host or device pointer */

cublasStatus_t cublasIzamax_v2 (
    cublasHandle_t handle,
    int n,
    const(cuDoubleComplex)* x,
    int incx,
    int* result); /* host or device pointer */

cublasStatus_t cublasIsamin_v2 (
    cublasHandle_t handle,
    int n,
    const(float)* x,
    int incx,
    int* result); /* host or device pointer */

cublasStatus_t cublasIdamin_v2 (
    cublasHandle_t handle,
    int n,
    const(double)* x,
    int incx,
    int* result); /* host or device pointer */

cublasStatus_t cublasIcamin_v2 (
    cublasHandle_t handle,
    int n,
    const(cuComplex)* x,
    int incx,
    int* result); /* host or device pointer */

cublasStatus_t cublasIzamin_v2 (
    cublasHandle_t handle,
    int n,
    const(cuDoubleComplex)* x,
    int incx,
    int* result); /* host or device pointer */

cublasStatus_t cublasSasum_v2 (
    cublasHandle_t handle,
    int n,
    const(float)* x,
    int incx,
    float* result); /* host or device pointer */

cublasStatus_t cublasDasum_v2 (
    cublasHandle_t handle,
    int n,
    const(double)* x,
    int incx,
    double* result); /* host or device pointer */

cublasStatus_t cublasScasum_v2 (
    cublasHandle_t handle,
    int n,
    const(cuComplex)* x,
    int incx,
    float* result); /* host or device pointer */

cublasStatus_t cublasDzasum_v2 (
    cublasHandle_t handle,
    int n,
    const(cuDoubleComplex)* x,
    int incx,
    double* result); /* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasSrot_v2 (
    cublasHandle_t handle,
    int n,
    float* x,
    int incx,
    float* y,
    int incy,
    const(float)* c,
    const(float)* s); /* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasDrot_v2 (
    cublasHandle_t handle,
    int n,
    double* x,
    int incx,
    double* y,
    int incy,
    const(double)* c,
    const(double)* s); /* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasCrot_v2 (
    cublasHandle_t handle,
    int n,
    cuComplex* x,
    int incx,
    cuComplex* y,
    int incy,
    const(float)* c,
    const(cuComplex)* s); /* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasCsrot_v2 (
    cublasHandle_t handle,
    int n,
    cuComplex* x,
    int incx,
    cuComplex* y,
    int incy,
    const(float)* c,
    const(float)* s); /* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasZrot_v2 (
    cublasHandle_t handle,
    int n,
    cuDoubleComplex* x,
    int incx,
    cuDoubleComplex* y,
    int incy,
    const(double)* c,
    const(cuDoubleComplex)* s); /* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasZdrot_v2 (
    cublasHandle_t handle,
    int n,
    cuDoubleComplex* x,
    int incx,
    cuDoubleComplex* y,
    int incy,
    const(double)* c,
    const(double)* s); /* host or device pointer */

/* host or device pointer */
/* host or device pointer */
/* host or device pointer */
cublasStatus_t cublasSrotg_v2 (
    cublasHandle_t handle,
    float* a,
    float* b,
    float* c,
    float* s); /* host or device pointer */

/* host or device pointer */
/* host or device pointer */
/* host or device pointer */
cublasStatus_t cublasDrotg_v2 (
    cublasHandle_t handle,
    double* a,
    double* b,
    double* c,
    double* s); /* host or device pointer */

/* host or device pointer */
/* host or device pointer */
/* host or device pointer */
cublasStatus_t cublasCrotg_v2 (
    cublasHandle_t handle,
    cuComplex* a,
    cuComplex* b,
    float* c,
    cuComplex* s); /* host or device pointer */

/* host or device pointer */
/* host or device pointer */
/* host or device pointer */
cublasStatus_t cublasZrotg_v2 (
    cublasHandle_t handle,
    cuDoubleComplex* a,
    cuDoubleComplex* b,
    double* c,
    cuDoubleComplex* s); /* host or device pointer */

cublasStatus_t cublasSrotm_v2 (
    cublasHandle_t handle,
    int n,
    float* x,
    int incx,
    float* y,
    int incy,
    const(float)* param); /* host or device pointer */

cublasStatus_t cublasDrotm_v2 (
    cublasHandle_t handle,
    int n,
    double* x,
    int incx,
    double* y,
    int incy,
    const(double)* param); /* host or device pointer */

/* host or device pointer */
/* host or device pointer */
/* host or device pointer */
/* host or device pointer */
cublasStatus_t cublasSrotmg_v2 (
    cublasHandle_t handle,
    float* d1,
    float* d2,
    float* x1,
    const(float)* y1,
    float* param); /* host or device pointer */

/* host or device pointer */
/* host or device pointer */
/* host or device pointer */
/* host or device pointer */
cublasStatus_t cublasDrotmg_v2 (
    cublasHandle_t handle,
    double* d1,
    double* d2,
    double* x1,
    const(double)* y1,
    double* param
    );

/* host or device pointer */

/* --------------- CUBLAS BLAS2 functions  ---------------- */

/* GEMV */

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasSgemv_v2 (
    cublasHandle_t handle,
    cublasOperation_t trans,
    int m,
    int n,
    const(float)* alpha,
    const(float)* A,
    int lda,
    const(float)* x,
    int incx,
    const(float)* beta,
    float* y,
    int incy);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasDgemv_v2 (
    cublasHandle_t handle,
    cublasOperation_t trans,
    int m,
    int n,
    const(double)* alpha,
    const(double)* A,
    int lda,
    const(double)* x,
    int incx,
    const(double)* beta,
    double* y,
    int incy);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasCgemv_v2 (
    cublasHandle_t handle,
    cublasOperation_t trans,
    int m,
    int n,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* x,
    int incx,
    const(cuComplex)* beta,
    cuComplex* y,
    int incy);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasZgemv_v2 (
    cublasHandle_t handle,
    cublasOperation_t trans,
    int m,
    int n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* x,
    int incx,
    const(cuDoubleComplex)* beta,
    cuDoubleComplex* y,
    int incy);

/* GBMV */

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasSgbmv_v2 (
    cublasHandle_t handle,
    cublasOperation_t trans,
    int m,
    int n,
    int kl,
    int ku,
    const(float)* alpha,
    const(float)* A,
    int lda,
    const(float)* x,
    int incx,
    const(float)* beta,
    float* y,
    int incy);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasDgbmv_v2 (
    cublasHandle_t handle,
    cublasOperation_t trans,
    int m,
    int n,
    int kl,
    int ku,
    const(double)* alpha,
    const(double)* A,
    int lda,
    const(double)* x,
    int incx,
    const(double)* beta,
    double* y,
    int incy);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasCgbmv_v2 (
    cublasHandle_t handle,
    cublasOperation_t trans,
    int m,
    int n,
    int kl,
    int ku,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* x,
    int incx,
    const(cuComplex)* beta,
    cuComplex* y,
    int incy);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasZgbmv_v2 (
    cublasHandle_t handle,
    cublasOperation_t trans,
    int m,
    int n,
    int kl,
    int ku,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* x,
    int incx,
    const(cuDoubleComplex)* beta,
    cuDoubleComplex* y,
    int incy);

/* TRMV */
cublasStatus_t cublasStrmv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    const(float)* A,
    int lda,
    float* x,
    int incx);

cublasStatus_t cublasDtrmv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    const(double)* A,
    int lda,
    double* x,
    int incx);

cublasStatus_t cublasCtrmv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    const(cuComplex)* A,
    int lda,
    cuComplex* x,
    int incx);

cublasStatus_t cublasZtrmv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    const(cuDoubleComplex)* A,
    int lda,
    cuDoubleComplex* x,
    int incx);

/* TBMV */
cublasStatus_t cublasStbmv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    int k,
    const(float)* A,
    int lda,
    float* x,
    int incx);

cublasStatus_t cublasDtbmv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    int k,
    const(double)* A,
    int lda,
    double* x,
    int incx);

cublasStatus_t cublasCtbmv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    int k,
    const(cuComplex)* A,
    int lda,
    cuComplex* x,
    int incx);

cublasStatus_t cublasZtbmv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    int k,
    const(cuDoubleComplex)* A,
    int lda,
    cuDoubleComplex* x,
    int incx);

/* TPMV */
cublasStatus_t cublasStpmv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    const(float)* AP,
    float* x,
    int incx);

cublasStatus_t cublasDtpmv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    const(double)* AP,
    double* x,
    int incx);

cublasStatus_t cublasCtpmv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    const(cuComplex)* AP,
    cuComplex* x,
    int incx);

cublasStatus_t cublasZtpmv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    const(cuDoubleComplex)* AP,
    cuDoubleComplex* x,
    int incx);

/* TRSV */
cublasStatus_t cublasStrsv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    const(float)* A,
    int lda,
    float* x,
    int incx);

cublasStatus_t cublasDtrsv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    const(double)* A,
    int lda,
    double* x,
    int incx);

cublasStatus_t cublasCtrsv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    const(cuComplex)* A,
    int lda,
    cuComplex* x,
    int incx);

cublasStatus_t cublasZtrsv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    const(cuDoubleComplex)* A,
    int lda,
    cuDoubleComplex* x,
    int incx);

/* TPSV */
cublasStatus_t cublasStpsv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    const(float)* AP,
    float* x,
    int incx);

cublasStatus_t cublasDtpsv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    const(double)* AP,
    double* x,
    int incx);

cublasStatus_t cublasCtpsv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    const(cuComplex)* AP,
    cuComplex* x,
    int incx);

cublasStatus_t cublasZtpsv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    const(cuDoubleComplex)* AP,
    cuDoubleComplex* x,
    int incx);

/* TBSV */
cublasStatus_t cublasStbsv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    int k,
    const(float)* A,
    int lda,
    float* x,
    int incx);

cublasStatus_t cublasDtbsv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    int k,
    const(double)* A,
    int lda,
    double* x,
    int incx);

cublasStatus_t cublasCtbsv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    int k,
    const(cuComplex)* A,
    int lda,
    cuComplex* x,
    int incx);

cublasStatus_t cublasZtbsv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int n,
    int k,
    const(cuDoubleComplex)* A,
    int lda,
    cuDoubleComplex* x,
    int incx);

/* SYMV/HEMV */

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasSsymv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(float)* alpha,
    const(float)* A,
    int lda,
    const(float)* x,
    int incx,
    const(float)* beta,
    float* y,
    int incy);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasDsymv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(double)* alpha,
    const(double)* A,
    int lda,
    const(double)* x,
    int incx,
    const(double)* beta,
    double* y,
    int incy);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasCsymv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* x,
    int incx,
    const(cuComplex)* beta,
    cuComplex* y,
    int incy);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasZsymv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* x,
    int incx,
    const(cuDoubleComplex)* beta,
    cuDoubleComplex* y,
    int incy);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasChemv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* x,
    int incx,
    const(cuComplex)* beta,
    cuComplex* y,
    int incy);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasZhemv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* x,
    int incx,
    const(cuDoubleComplex)* beta,
    cuDoubleComplex* y,
    int incy);

/* SBMV/HBMV */

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasSsbmv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    int k,
    const(float)* alpha,
    const(float)* A,
    int lda,
    const(float)* x,
    int incx,
    const(float)* beta,
    float* y,
    int incy);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasDsbmv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    int k,
    const(double)* alpha,
    const(double)* A,
    int lda,
    const(double)* x,
    int incx,
    const(double)* beta,
    double* y,
    int incy);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasChbmv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    int k,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* x,
    int incx,
    const(cuComplex)* beta,
    cuComplex* y,
    int incy);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasZhbmv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    int k,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* x,
    int incx,
    const(cuDoubleComplex)* beta,
    cuDoubleComplex* y,
    int incy);

/* SPMV/HPMV */

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasSspmv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(float)* alpha,
    const(float)* AP,
    const(float)* x,
    int incx,
    const(float)* beta,
    float* y,
    int incy);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasDspmv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(double)* alpha,
    const(double)* AP,
    const(double)* x,
    int incx,
    const(double)* beta,
    double* y,
    int incy);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasChpmv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(cuComplex)* alpha,
    const(cuComplex)* AP,
    const(cuComplex)* x,
    int incx,
    const(cuComplex)* beta,
    cuComplex* y,
    int incy);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasZhpmv_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* AP,
    const(cuDoubleComplex)* x,
    int incx,
    const(cuDoubleComplex)* beta,
    cuDoubleComplex* y,
    int incy);

/* GER */

/* host or device pointer */
cublasStatus_t cublasSger_v2 (
    cublasHandle_t handle,
    int m,
    int n,
    const(float)* alpha,
    const(float)* x,
    int incx,
    const(float)* y,
    int incy,
    float* A,
    int lda);

/* host or device pointer */
cublasStatus_t cublasDger_v2 (
    cublasHandle_t handle,
    int m,
    int n,
    const(double)* alpha,
    const(double)* x,
    int incx,
    const(double)* y,
    int incy,
    double* A,
    int lda);

/* host or device pointer */
cublasStatus_t cublasCgeru_v2 (
    cublasHandle_t handle,
    int m,
    int n,
    const(cuComplex)* alpha,
    const(cuComplex)* x,
    int incx,
    const(cuComplex)* y,
    int incy,
    cuComplex* A,
    int lda);

/* host or device pointer */
cublasStatus_t cublasCgerc_v2 (
    cublasHandle_t handle,
    int m,
    int n,
    const(cuComplex)* alpha,
    const(cuComplex)* x,
    int incx,
    const(cuComplex)* y,
    int incy,
    cuComplex* A,
    int lda);

/* host or device pointer */
cublasStatus_t cublasZgeru_v2 (
    cublasHandle_t handle,
    int m,
    int n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* x,
    int incx,
    const(cuDoubleComplex)* y,
    int incy,
    cuDoubleComplex* A,
    int lda);

/* host or device pointer */
cublasStatus_t cublasZgerc_v2 (
    cublasHandle_t handle,
    int m,
    int n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* x,
    int incx,
    const(cuDoubleComplex)* y,
    int incy,
    cuDoubleComplex* A,
    int lda);

/* SYR/HER */

/* host or device pointer */
cublasStatus_t cublasSsyr_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(float)* alpha,
    const(float)* x,
    int incx,
    float* A,
    int lda);

/* host or device pointer */
cublasStatus_t cublasDsyr_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(double)* alpha,
    const(double)* x,
    int incx,
    double* A,
    int lda);

/* host or device pointer */
cublasStatus_t cublasCsyr_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(cuComplex)* alpha,
    const(cuComplex)* x,
    int incx,
    cuComplex* A,
    int lda);

/* host or device pointer */
cublasStatus_t cublasZsyr_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* x,
    int incx,
    cuDoubleComplex* A,
    int lda);

/* host or device pointer */
cublasStatus_t cublasCher_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(float)* alpha,
    const(cuComplex)* x,
    int incx,
    cuComplex* A,
    int lda);

/* host or device pointer */
cublasStatus_t cublasZher_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(double)* alpha,
    const(cuDoubleComplex)* x,
    int incx,
    cuDoubleComplex* A,
    int lda);

/* SPR/HPR */

/* host or device pointer */
cublasStatus_t cublasSspr_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(float)* alpha,
    const(float)* x,
    int incx,
    float* AP);

/* host or device pointer */
cublasStatus_t cublasDspr_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(double)* alpha,
    const(double)* x,
    int incx,
    double* AP);

/* host or device pointer */
cublasStatus_t cublasChpr_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(float)* alpha,
    const(cuComplex)* x,
    int incx,
    cuComplex* AP);

/* host or device pointer */
cublasStatus_t cublasZhpr_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(double)* alpha,
    const(cuDoubleComplex)* x,
    int incx,
    cuDoubleComplex* AP);

/* SYR2/HER2 */

/* host or device pointer */
cublasStatus_t cublasSsyr2_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(float)* alpha,
    const(float)* x,
    int incx,
    const(float)* y,
    int incy,
    float* A,
    int lda);

/* host or device pointer */
cublasStatus_t cublasDsyr2_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(double)* alpha,
    const(double)* x,
    int incx,
    const(double)* y,
    int incy,
    double* A,
    int lda);

/* host or device pointer */
cublasStatus_t cublasCsyr2_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(cuComplex)* alpha,
    const(cuComplex)* x,
    int incx,
    const(cuComplex)* y,
    int incy,
    cuComplex* A,
    int lda);

/* host or device pointer */
cublasStatus_t cublasZsyr2_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* x,
    int incx,
    const(cuDoubleComplex)* y,
    int incy,
    cuDoubleComplex* A,
    int lda);

/* host or device pointer */
cublasStatus_t cublasCher2_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(cuComplex)* alpha,
    const(cuComplex)* x,
    int incx,
    const(cuComplex)* y,
    int incy,
    cuComplex* A,
    int lda);

/* host or device pointer */
cublasStatus_t cublasZher2_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* x,
    int incx,
    const(cuDoubleComplex)* y,
    int incy,
    cuDoubleComplex* A,
    int lda);

/* SPR2/HPR2 */

/* host or device pointer */
cublasStatus_t cublasSspr2_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(float)* alpha,
    const(float)* x,
    int incx,
    const(float)* y,
    int incy,
    float* AP);

/* host or device pointer */
cublasStatus_t cublasDspr2_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(double)* alpha,
    const(double)* x,
    int incx,
    const(double)* y,
    int incy,
    double* AP);

/* host or device pointer */
cublasStatus_t cublasChpr2_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(cuComplex)* alpha,
    const(cuComplex)* x,
    int incx,
    const(cuComplex)* y,
    int incy,
    cuComplex* AP);

/* host or device pointer */
cublasStatus_t cublasZhpr2_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* x,
    int incx,
    const(cuDoubleComplex)* y,
    int incy,
    cuDoubleComplex* AP);

/* ---------------- CUBLAS BLAS3 functions ---------------- */

/* GEMM */

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasSgemm_v2 (
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int m,
    int n,
    int k,
    const(float)* alpha,
    const(float)* A,
    int lda,
    const(float)* B,
    int ldb,
    const(float)* beta,
    float* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasDgemm_v2 (
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int m,
    int n,
    int k,
    const(double)* alpha,
    const(double)* A,
    int lda,
    const(double)* B,
    int ldb,
    const(double)* beta,
    double* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasCgemm_v2 (
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int m,
    int n,
    int k,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* B,
    int ldb,
    const(cuComplex)* beta,
    cuComplex* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasZgemm_v2 (
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int m,
    int n,
    int k,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* B,
    int ldb,
    const(cuDoubleComplex)* beta,
    cuDoubleComplex* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasHgemm (
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int m,
    int n,
    int k,
    const(half)* alpha,
    const(half)* A,
    int lda,
    const(half)* B,
    int ldb,
    const(half)* beta,
    half* C,
    int ldc);

/* IO in FP16/FP32, computation in float */

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasSgemmEx (
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int m,
    int n,
    int k,
    const(float)* alpha,
    const(void)* A,
    cublasDataType_t Atype,
    int lda,
    const(void)* B,
    cublasDataType_t Btype,
    int ldb,
    const(float)* beta,
    void* C,
    cublasDataType_t Ctype,
    int ldc);

/* SYRK */

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasSsyrk_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int n,
    int k,
    const(float)* alpha,
    const(float)* A,
    int lda,
    const(float)* beta,
    float* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasDsyrk_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int n,
    int k,
    const(double)* alpha,
    const(double)* A,
    int lda,
    const(double)* beta,
    double* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasCsyrk_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int n,
    int k,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* beta,
    cuComplex* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasZsyrk_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int n,
    int k,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* beta,
    cuDoubleComplex* C,
    int ldc);

/* HERK */

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasCherk_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int n,
    int k,
    const(float)* alpha,
    const(cuComplex)* A,
    int lda,
    const(float)* beta,
    cuComplex* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasZherk_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int n,
    int k,
    const(double)* alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(double)* beta,
    cuDoubleComplex* C,
    int ldc);

/* SYR2K */

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasSsyr2k_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int n,
    int k,
    const(float)* alpha,
    const(float)* A,
    int lda,
    const(float)* B,
    int ldb,
    const(float)* beta,
    float* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasDsyr2k_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int n,
    int k,
    const(double)* alpha,
    const(double)* A,
    int lda,
    const(double)* B,
    int ldb,
    const(double)* beta,
    double* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasCsyr2k_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int n,
    int k,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* B,
    int ldb,
    const(cuComplex)* beta,
    cuComplex* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasZsyr2k_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int n,
    int k,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* B,
    int ldb,
    const(cuDoubleComplex)* beta,
    cuDoubleComplex* C,
    int ldc);

/* HER2K */

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasCher2k_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int n,
    int k,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* B,
    int ldb,
    const(float)* beta,
    cuComplex* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasZher2k_v2 (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int n,
    int k,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* B,
    int ldb,
    const(double)* beta,
    cuDoubleComplex* C,
    int ldc);

/* SYRKX : eXtended SYRK*/

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasSsyrkx (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int n,
    int k,
    const(float)* alpha,
    const(float)* A,
    int lda,
    const(float)* B,
    int ldb,
    const(float)* beta,
    float* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasDsyrkx (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int n,
    int k,
    const(double)* alpha,
    const(double)* A,
    int lda,
    const(double)* B,
    int ldb,
    const(double)* beta,
    double* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasCsyrkx (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int n,
    int k,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* B,
    int ldb,
    const(cuComplex)* beta,
    cuComplex* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasZsyrkx (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int n,
    int k,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* B,
    int ldb,
    const(cuDoubleComplex)* beta,
    cuDoubleComplex* C,
    int ldc);

/* HERKX : eXtended HERK */

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasCherkx (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int n,
    int k,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* B,
    int ldb,
    const(float)* beta,
    cuComplex* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasZherkx (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    int n,
    int k,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* B,
    int ldb,
    const(double)* beta,
    cuDoubleComplex* C,
    int ldc);

/* SYMM */

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasSsymm_v2 (
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    int m,
    int n,
    const(float)* alpha,
    const(float)* A,
    int lda,
    const(float)* B,
    int ldb,
    const(float)* beta,
    float* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasDsymm_v2 (
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    int m,
    int n,
    const(double)* alpha,
    const(double)* A,
    int lda,
    const(double)* B,
    int ldb,
    const(double)* beta,
    double* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasCsymm_v2 (
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    int m,
    int n,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* B,
    int ldb,
    const(cuComplex)* beta,
    cuComplex* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasZsymm_v2 (
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    int m,
    int n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* B,
    int ldb,
    const(cuDoubleComplex)* beta,
    cuDoubleComplex* C,
    int ldc);

/* HEMM */

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasChemm_v2 (
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    int m,
    int n,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* B,
    int ldb,
    const(cuComplex)* beta,
    cuComplex* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasZhemm_v2 (
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    int m,
    int n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* B,
    int ldb,
    const(cuDoubleComplex)* beta,
    cuDoubleComplex* C,
    int ldc);

/* TRSM */

/* host or device pointer */
cublasStatus_t cublasStrsm_v2 (
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int m,
    int n,
    const(float)* alpha,
    const(float)* A,
    int lda,
    float* B,
    int ldb);

/* host or device pointer */
cublasStatus_t cublasDtrsm_v2 (
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int m,
    int n,
    const(double)* alpha,
    const(double)* A,
    int lda,
    double* B,
    int ldb);

/* host or device pointer */
cublasStatus_t cublasCtrsm_v2 (
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int m,
    int n,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    int lda,
    cuComplex* B,
    int ldb);

/* host or device pointer */
cublasStatus_t cublasZtrsm_v2 (
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int m,
    int n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    int lda,
    cuDoubleComplex* B,
    int ldb);

/* TRMM */

/* host or device pointer */
cublasStatus_t cublasStrmm_v2 (
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int m,
    int n,
    const(float)* alpha,
    const(float)* A,
    int lda,
    const(float)* B,
    int ldb,
    float* C,
    int ldc);

/* host or device pointer */
cublasStatus_t cublasDtrmm_v2 (
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int m,
    int n,
    const(double)* alpha,
    const(double)* A,
    int lda,
    const(double)* B,
    int ldb,
    double* C,
    int ldc);

/* host or device pointer */
cublasStatus_t cublasCtrmm_v2 (
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int m,
    int n,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* B,
    int ldb,
    cuComplex* C,
    int ldc);

/* host or device pointer */
cublasStatus_t cublasZtrmm_v2 (
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int m,
    int n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* B,
    int ldb,
    cuDoubleComplex* C,
    int ldc);

/* BATCH GEMM */

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasSgemmBatched (
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int m,
    int n,
    int k,
    const(float)* alpha,
    const(float)** Aarray,
    int lda,
    const(float)** Barray,
    int ldb,
    const(float)* beta,
    float** Carray,
    int ldc,
    int batchCount);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasDgemmBatched (
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int m,
    int n,
    int k,
    const(double)* alpha,
    const(double)** Aarray,
    int lda,
    const(double)** Barray,
    int ldb,
    const(double)* beta,
    double** Carray,
    int ldc,
    int batchCount);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasCgemmBatched (
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int m,
    int n,
    int k,
    const(cuComplex)* alpha,
    const(cuComplex)** Aarray,
    int lda,
    const(cuComplex)** Barray,
    int ldb,
    const(cuComplex)* beta,
    cuComplex** Carray,
    int ldc,
    int batchCount);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasZgemmBatched (
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int m,
    int n,
    int k,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)** Aarray,
    int lda,
    const(cuDoubleComplex)** Barray,
    int ldb,
    const(cuDoubleComplex)* beta,
    cuDoubleComplex** Carray,
    int ldc,
    int batchCount);

/* ---------------- CUBLAS BLAS-like extension ---------------- */
/* GEAM */

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasSgeam (
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int m,
    int n,
    const(float)* alpha,
    const(float)* A,
    int lda,
    const(float)* beta,
    const(float)* B,
    int ldb,
    float* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasDgeam (
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int m,
    int n,
    const(double)* alpha,
    const(double)* A,
    int lda,
    const(double)* beta,
    const(double)* B,
    int ldb,
    double* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasCgeam (
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int m,
    int n,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* beta,
    const(cuComplex)* B,
    int ldb,
    cuComplex* C,
    int ldc);

/* host or device pointer */

/* host or device pointer */
cublasStatus_t cublasZgeam (
    cublasHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    int m,
    int n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* beta,
    const(cuDoubleComplex)* B,
    int ldb,
    cuDoubleComplex* C,
    int ldc);

/* Batched LU - GETRF*/

/*Device pointer*/

/*Device Pointer*/
/*Device Pointer*/
cublasStatus_t cublasSgetrfBatched (
    cublasHandle_t handle,
    int n,
    float** A,
    int lda,
    int* P,
    int* info,
    int batchSize);

/*Device pointer*/

/*Device Pointer*/
/*Device Pointer*/
cublasStatus_t cublasDgetrfBatched (
    cublasHandle_t handle,
    int n,
    double** A,
    int lda,
    int* P,
    int* info,
    int batchSize);

/*Device pointer*/

/*Device Pointer*/
/*Device Pointer*/
cublasStatus_t cublasCgetrfBatched (
    cublasHandle_t handle,
    int n,
    cuComplex** A,
    int lda,
    int* P,
    int* info,
    int batchSize);

/*Device pointer*/

/*Device Pointer*/
/*Device Pointer*/
cublasStatus_t cublasZgetrfBatched (
    cublasHandle_t handle,
    int n,
    cuDoubleComplex** A,
    int lda,
    int* P,
    int* info,
    int batchSize);

/* Batched inversion based on LU factorization from getrf */

/*Device pointer*/

/*Device pointer*/
/*Device pointer*/
cublasStatus_t cublasSgetriBatched (
    cublasHandle_t handle,
    int n,
    const(float)** A,
    int lda,
    const(int)* P,
    float** C,
    int ldc,
    int* info,
    int batchSize);

/*Device pointer*/

/*Device pointer*/
/*Device pointer*/
cublasStatus_t cublasDgetriBatched (
    cublasHandle_t handle,
    int n,
    const(double)** A,
    int lda,
    const(int)* P,
    double** C,
    int ldc,
    int* info,
    int batchSize);

/*Device pointer*/

/*Device pointer*/
/*Device pointer*/
cublasStatus_t cublasCgetriBatched (
    cublasHandle_t handle,
    int n,
    const(cuComplex)** A,
    int lda,
    const(int)* P,
    cuComplex** C,
    int ldc,
    int* info,
    int batchSize);

/*Device pointer*/

/*Device pointer*/
/*Device pointer*/
cublasStatus_t cublasZgetriBatched (
    cublasHandle_t handle,
    int n,
    const(cuDoubleComplex)** A,
    int lda,
    const(int)* P,
    cuDoubleComplex** C,
    int ldc,
    int* info,
    int batchSize);

/* Batched solver based on LU factorization from getrf */

cublasStatus_t cublasSgetrsBatched (
    cublasHandle_t handle,
    cublasOperation_t trans,
    int n,
    int nrhs,
    const(float)** Aarray,
    int lda,
    const(int)* devIpiv,
    float** Barray,
    int ldb,
    int* info,
    int batchSize);

cublasStatus_t cublasDgetrsBatched (
    cublasHandle_t handle,
    cublasOperation_t trans,
    int n,
    int nrhs,
    const(double)** Aarray,
    int lda,
    const(int)* devIpiv,
    double** Barray,
    int ldb,
    int* info,
    int batchSize);

cublasStatus_t cublasCgetrsBatched (
    cublasHandle_t handle,
    cublasOperation_t trans,
    int n,
    int nrhs,
    const(cuComplex)** Aarray,
    int lda,
    const(int)* devIpiv,
    cuComplex** Barray,
    int ldb,
    int* info,
    int batchSize);

cublasStatus_t cublasZgetrsBatched (
    cublasHandle_t handle,
    cublasOperation_t trans,
    int n,
    int nrhs,
    const(cuDoubleComplex)** Aarray,
    int lda,
    const(int)* devIpiv,
    cuDoubleComplex** Barray,
    int ldb,
    int* info,
    int batchSize);

/* TRSM - Batched Triangular Solver */

/*Host or Device Pointer*/
cublasStatus_t cublasStrsmBatched (
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int m,
    int n,
    const(float)* alpha,
    const(float)** A,
    int lda,
    float** B,
    int ldb,
    int batchCount);

/*Host or Device Pointer*/
cublasStatus_t cublasDtrsmBatched (
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int m,
    int n,
    const(double)* alpha,
    const(double)** A,
    int lda,
    double** B,
    int ldb,
    int batchCount);

/*Host or Device Pointer*/
cublasStatus_t cublasCtrsmBatched (
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int m,
    int n,
    const(cuComplex)* alpha,
    const(cuComplex)** A,
    int lda,
    cuComplex** B,
    int ldb,
    int batchCount);

/*Host or Device Pointer*/
cublasStatus_t cublasZtrsmBatched (
    cublasHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    int m,
    int n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)** A,
    int lda,
    cuDoubleComplex** B,
    int ldb,
    int batchCount);

/* Batched - MATINV*/

/*Device pointer*/

/*Device pointer*/

/*Device Pointer*/
cublasStatus_t cublasSmatinvBatched (
    cublasHandle_t handle,
    int n,
    const(float)** A,
    int lda,
    float** Ainv,
    int lda_inv,
    int* info,
    int batchSize);

/*Device pointer*/

/*Device pointer*/

/*Device Pointer*/
cublasStatus_t cublasDmatinvBatched (
    cublasHandle_t handle,
    int n,
    const(double)** A,
    int lda,
    double** Ainv,
    int lda_inv,
    int* info,
    int batchSize);

/*Device pointer*/

/*Device pointer*/

/*Device Pointer*/
cublasStatus_t cublasCmatinvBatched (
    cublasHandle_t handle,
    int n,
    const(cuComplex)** A,
    int lda,
    cuComplex** Ainv,
    int lda_inv,
    int* info,
    int batchSize);

/*Device pointer*/

/*Device pointer*/

/*Device Pointer*/
cublasStatus_t cublasZmatinvBatched (
    cublasHandle_t handle,
    int n,
    const(cuDoubleComplex)** A,
    int lda,
    cuDoubleComplex** Ainv,
    int lda_inv,
    int* info,
    int batchSize);

/* Batch QR Factorization */

/*Device pointer*/

/* Device pointer*/
cublasStatus_t cublasSgeqrfBatched (
    cublasHandle_t handle,
    int m,
    int n,
    float** Aarray,
    int lda,
    float** TauArray,
    int* info,
    int batchSize);

/*Device pointer*/

/* Device pointer*/
cublasStatus_t cublasDgeqrfBatched (
    cublasHandle_t handle,
    int m,
    int n,
    double** Aarray,
    int lda,
    double** TauArray,
    int* info,
    int batchSize);

/*Device pointer*/

/* Device pointer*/
cublasStatus_t cublasCgeqrfBatched (
    cublasHandle_t handle,
    int m,
    int n,
    cuComplex** Aarray,
    int lda,
    cuComplex** TauArray,
    int* info,
    int batchSize);

/*Device pointer*/

/* Device pointer*/
cublasStatus_t cublasZgeqrfBatched (
    cublasHandle_t handle,
    int m,
    int n,
    cuDoubleComplex** Aarray,
    int lda,
    cuDoubleComplex** TauArray,
    int* info,
    int batchSize);

/* Least Square Min only m >= n and Non-transpose supported */

/*Device pointer*/

/* Device pointer*/

/* Device pointer*/
cublasStatus_t cublasSgelsBatched (
    cublasHandle_t handle,
    cublasOperation_t trans,
    int m,
    int n,
    int nrhs,
    float** Aarray,
    int lda,
    float** Carray,
    int ldc,
    int* info,
    int* devInfoArray,
    int batchSize);

/*Device pointer*/

/* Device pointer*/

/* Device pointer*/
cublasStatus_t cublasDgelsBatched (
    cublasHandle_t handle,
    cublasOperation_t trans,
    int m,
    int n,
    int nrhs,
    double** Aarray,
    int lda,
    double** Carray,
    int ldc,
    int* info,
    int* devInfoArray,
    int batchSize);

/*Device pointer*/

/* Device pointer*/
cublasStatus_t cublasCgelsBatched (
    cublasHandle_t handle,
    cublasOperation_t trans,
    int m,
    int n,
    int nrhs,
    cuComplex** Aarray,
    int lda,
    cuComplex** Carray,
    int ldc,
    int* info,
    int* devInfoArray,
    int batchSize);

/*Device pointer*/

/* Device pointer*/
cublasStatus_t cublasZgelsBatched (
    cublasHandle_t handle,
    cublasOperation_t trans,
    int m,
    int n,
    int nrhs,
    cuDoubleComplex** Aarray,
    int lda,
    cuDoubleComplex** Carray,
    int ldc,
    int* info,
    int* devInfoArray,
    int batchSize);

/* DGMM */
cublasStatus_t cublasSdgmm (
    cublasHandle_t handle,
    cublasSideMode_t mode,
    int m,
    int n,
    const(float)* A,
    int lda,
    const(float)* x,
    int incx,
    float* C,
    int ldc);

cublasStatus_t cublasDdgmm (
    cublasHandle_t handle,
    cublasSideMode_t mode,
    int m,
    int n,
    const(double)* A,
    int lda,
    const(double)* x,
    int incx,
    double* C,
    int ldc);

cublasStatus_t cublasCdgmm (
    cublasHandle_t handle,
    cublasSideMode_t mode,
    int m,
    int n,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* x,
    int incx,
    cuComplex* C,
    int ldc);

cublasStatus_t cublasZdgmm (
    cublasHandle_t handle,
    cublasSideMode_t mode,
    int m,
    int n,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* x,
    int incx,
    cuDoubleComplex* C,
    int ldc);

/* TPTTR : Triangular Pack format to Triangular format */
cublasStatus_t cublasStpttr (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(float)* AP,
    float* A,
    int lda);

cublasStatus_t cublasDtpttr (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(double)* AP,
    double* A,
    int lda);

cublasStatus_t cublasCtpttr (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(cuComplex)* AP,
    cuComplex* A,
    int lda);

cublasStatus_t cublasZtpttr (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(cuDoubleComplex)* AP,
    cuDoubleComplex* A,
    int lda);

/* TRTTP : Triangular format to Triangular Pack format */
cublasStatus_t cublasStrttp (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(float)* A,
    int lda,
    float* AP);

cublasStatus_t cublasDtrttp (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(double)* A,
    int lda,
    double* AP);

cublasStatus_t cublasCtrttp (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(cuComplex)* A,
    int lda,
    cuComplex* AP);

cublasStatus_t cublasZtrttp (
    cublasHandle_t handle,
    cublasFillMode_t uplo,
    int n,
    const(cuDoubleComplex)* A,
    int lda,
    cuDoubleComplex* AP);
