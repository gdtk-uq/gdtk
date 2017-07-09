/*
  D bindings for CUDA.
  Authors:    Prasun Anand
  Copyright:  Copyright (c) 2017, Prasun Anand. All rights reserved.
  License:    BSD 3-Clause License
*/
module cuda_d.cublastXt;

import cuda_d.cublas_api;
import cuda_d.cuComplex;

extern (C):

struct cublasXtContext;
alias cublasXtHandle_t = cublasXtContext*;

cublasStatus_t cublasXtCreate (cublasXtHandle_t* handle);
cublasStatus_t cublasXtDestroy (cublasXtHandle_t handle);
cublasStatus_t cublasXtGetNumBoards (int nbDevices, int* deviceId, int* nbBoards);
cublasStatus_t cublasXtMaxBoards (int* nbGpuBoards);
/* This routine selects the Gpus that the user want to use for CUBLAS-XT */
cublasStatus_t cublasXtDeviceSelect (cublasXtHandle_t handle, int nbDevices, int* deviceId);

/* This routine allows to change the dimension of the tiles ( blockDim x blockDim ) */
cublasStatus_t cublasXtSetBlockDim (cublasXtHandle_t handle, int blockDim);
cublasStatus_t cublasXtGetBlockDim (cublasXtHandle_t handle, int* blockDim);

enum cublasXtPinnedMemMode_t
{
    CUBLASXT_PINNING_DISABLED = 0,
    CUBLASXT_PINNING_ENABLED = 1
}

cublasStatus_t cublasXtGetPinningMemMode (cublasXtHandle_t handle, cublasXtPinnedMemMode_t* mode);
cublasStatus_t cublasXtSetPinningMemMode (cublasXtHandle_t handle, cublasXtPinnedMemMode_t mode);

enum cublasXtOpType_t
{
    CUBLASXT_FLOAT = 0,
    CUBLASXT_DOUBLE = 1,
    CUBLASXT_COMPLEX = 2,
    CUBLASXT_DOUBLECOMPLEX = 3
}

enum cublasXtBlasOp_t
{
    CUBLASXT_GEMM = 0,
    CUBLASXT_SYRK = 1,
    CUBLASXT_HERK = 2,
    CUBLASXT_SYMM = 3,
    CUBLASXT_HEMM = 4,
    CUBLASXT_TRSM = 5,
    CUBLASXT_SYR2K = 6,
    CUBLASXT_HER2K = 7,

    CUBLASXT_SPMM = 8,
    CUBLASXT_SYRKX = 9,
    CUBLASXT_HERKX = 10,
    CUBLASXT_TRMM = 11,
    CUBLASXT_ROUTINE_MAX = 12
}

/* Currently only 32-bit integer BLAS routines are supported */
cublasStatus_t cublasXtSetCpuRoutine (cublasXtHandle_t handle, cublasXtBlasOp_t blasOp, cublasXtOpType_t type, void* blasFunctor);

/* Specified the percentage of work that should done by the CPU, default is 0 (no work) */
cublasStatus_t cublasXtSetCpuRatio (cublasXtHandle_t handle, cublasXtBlasOp_t blasOp, cublasXtOpType_t type, float ratio);

/* GEMM */
cublasStatus_t cublasXtSgemm (
    cublasXtHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    size_t m,
    size_t n,
    size_t k,
    const(float)* alpha,
    const(float)* A,
    size_t lda,
    const(float)* B,
    size_t ldb,
    const(float)* beta,
    float* C,
    size_t ldc);

cublasStatus_t cublasXtDgemm (
    cublasXtHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    size_t m,
    size_t n,
    size_t k,
    const(double)* alpha,
    const(double)* A,
    size_t lda,
    const(double)* B,
    size_t ldb,
    const(double)* beta,
    double* C,
    size_t ldc);

cublasStatus_t cublasXtCgemm (
    cublasXtHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    size_t m,
    size_t n,
    size_t k,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    size_t lda,
    const(cuComplex)* B,
    size_t ldb,
    const(cuComplex)* beta,
    cuComplex* C,
    size_t ldc);

cublasStatus_t cublasXtZgemm (
    cublasXtHandle_t handle,
    cublasOperation_t transa,
    cublasOperation_t transb,
    size_t m,
    size_t n,
    size_t k,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    size_t lda,
    const(cuDoubleComplex)* B,
    size_t ldb,
    const(cuDoubleComplex)* beta,
    cuDoubleComplex* C,
    size_t ldc);

/* ------------------------------------------------------- */
/* SYRK */
cublasStatus_t cublasXtSsyrk (
    cublasXtHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    size_t n,
    size_t k,
    const(float)* alpha,
    const(float)* A,
    size_t lda,
    const(float)* beta,
    float* C,
    size_t ldc);

cublasStatus_t cublasXtDsyrk (
    cublasXtHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    size_t n,
    size_t k,
    const(double)* alpha,
    const(double)* A,
    size_t lda,
    const(double)* beta,
    double* C,
    size_t ldc);

cublasStatus_t cublasXtCsyrk (
    cublasXtHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    size_t n,
    size_t k,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    size_t lda,
    const(cuComplex)* beta,
    cuComplex* C,
    size_t ldc);

cublasStatus_t cublasXtZsyrk (
    cublasXtHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    size_t n,
    size_t k,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    size_t lda,
    const(cuDoubleComplex)* beta,
    cuDoubleComplex* C,
    size_t ldc);

/* -------------------------------------------------------------------- */
/* HERK */
cublasStatus_t cublasXtCherk (
    cublasXtHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    size_t n,
    size_t k,
    const(float)* alpha,
    const(cuComplex)* A,
    size_t lda,
    const(float)* beta,
    cuComplex* C,
    size_t ldc);

cublasStatus_t cublasXtZherk (
    cublasXtHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    size_t n,
    size_t k,
    const(double)* alpha,
    const(cuDoubleComplex)* A,
    size_t lda,
    const(double)* beta,
    cuDoubleComplex* C,
    size_t ldc);

/* -------------------------------------------------------------------- */
/* SYR2K */
cublasStatus_t cublasXtSsyr2k (
    cublasXtHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    size_t n,
    size_t k,
    const(float)* alpha,
    const(float)* A,
    size_t lda,
    const(float)* B,
    size_t ldb,
    const(float)* beta,
    float* C,
    size_t ldc);

cublasStatus_t cublasXtDsyr2k (
    cublasXtHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    size_t n,
    size_t k,
    const(double)* alpha,
    const(double)* A,
    size_t lda,
    const(double)* B,
    size_t ldb,
    const(double)* beta,
    double* C,
    size_t ldc);

cublasStatus_t cublasXtCsyr2k (
    cublasXtHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    size_t n,
    size_t k,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    size_t lda,
    const(cuComplex)* B,
    size_t ldb,
    const(cuComplex)* beta,
    cuComplex* C,
    size_t ldc);

cublasStatus_t cublasXtZsyr2k (
    cublasXtHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    size_t n,
    size_t k,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    size_t lda,
    const(cuDoubleComplex)* B,
    size_t ldb,
    const(cuDoubleComplex)* beta,
    cuDoubleComplex* C,
    size_t ldc);

/* -------------------------------------------------------------------- */
/* HERKX : variant extension of HERK */
cublasStatus_t cublasXtCherkx (
    cublasXtHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    size_t n,
    size_t k,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    size_t lda,
    const(cuComplex)* B,
    size_t ldb,
    const(float)* beta,
    cuComplex* C,
    size_t ldc);

cublasStatus_t cublasXtZherkx (
    cublasXtHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    size_t n,
    size_t k,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    size_t lda,
    const(cuDoubleComplex)* B,
    size_t ldb,
    const(double)* beta,
    cuDoubleComplex* C,
    size_t ldc);

/* -------------------------------------------------------------------- */
/* TRSM */
cublasStatus_t cublasXtStrsm (
    cublasXtHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    size_t m,
    size_t n,
    const(float)* alpha,
    const(float)* A,
    size_t lda,
    float* B,
    size_t ldb);

cublasStatus_t cublasXtDtrsm (
    cublasXtHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    size_t m,
    size_t n,
    const(double)* alpha,
    const(double)* A,
    size_t lda,
    double* B,
    size_t ldb);

cublasStatus_t cublasXtCtrsm (
    cublasXtHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    size_t m,
    size_t n,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    size_t lda,
    cuComplex* B,
    size_t ldb);

cublasStatus_t cublasXtZtrsm (
    cublasXtHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    size_t m,
    size_t n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    size_t lda,
    cuDoubleComplex* B,
    size_t ldb);

/* -------------------------------------------------------------------- */
/* SYMM : Symmetric Multiply Matrix*/
cublasStatus_t cublasXtSsymm (
    cublasXtHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    size_t m,
    size_t n,
    const(float)* alpha,
    const(float)* A,
    size_t lda,
    const(float)* B,
    size_t ldb,
    const(float)* beta,
    float* C,
    size_t ldc);

cublasStatus_t cublasXtDsymm (
    cublasXtHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    size_t m,
    size_t n,
    const(double)* alpha,
    const(double)* A,
    size_t lda,
    const(double)* B,
    size_t ldb,
    const(double)* beta,
    double* C,
    size_t ldc);

cublasStatus_t cublasXtCsymm (
    cublasXtHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    size_t m,
    size_t n,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    size_t lda,
    const(cuComplex)* B,
    size_t ldb,
    const(cuComplex)* beta,
    cuComplex* C,
    size_t ldc);

cublasStatus_t cublasXtZsymm (
    cublasXtHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    size_t m,
    size_t n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    size_t lda,
    const(cuDoubleComplex)* B,
    size_t ldb,
    const(cuDoubleComplex)* beta,
    cuDoubleComplex* C,
    size_t ldc);

/* -------------------------------------------------------------------- */
/* HEMM : Hermitian Matrix Multiply */
cublasStatus_t cublasXtChemm (
    cublasXtHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    size_t m,
    size_t n,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    size_t lda,
    const(cuComplex)* B,
    size_t ldb,
    const(cuComplex)* beta,
    cuComplex* C,
    size_t ldc);

cublasStatus_t cublasXtZhemm (
    cublasXtHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    size_t m,
    size_t n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    size_t lda,
    const(cuDoubleComplex)* B,
    size_t ldb,
    const(cuDoubleComplex)* beta,
    cuDoubleComplex* C,
    size_t ldc);

/* -------------------------------------------------------------------- */
/* SYRKX : variant extension of SYRK  */
cublasStatus_t cublasXtSsyrkx (
    cublasXtHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    size_t n,
    size_t k,
    const(float)* alpha,
    const(float)* A,
    size_t lda,
    const(float)* B,
    size_t ldb,
    const(float)* beta,
    float* C,
    size_t ldc);

cublasStatus_t cublasXtDsyrkx (
    cublasXtHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    size_t n,
    size_t k,
    const(double)* alpha,
    const(double)* A,
    size_t lda,
    const(double)* B,
    size_t ldb,
    const(double)* beta,
    double* C,
    size_t ldc);

cublasStatus_t cublasXtCsyrkx (
    cublasXtHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    size_t n,
    size_t k,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    size_t lda,
    const(cuComplex)* B,
    size_t ldb,
    const(cuComplex)* beta,
    cuComplex* C,
    size_t ldc);

cublasStatus_t cublasXtZsyrkx (
    cublasXtHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    size_t n,
    size_t k,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    size_t lda,
    const(cuDoubleComplex)* B,
    size_t ldb,
    const(cuDoubleComplex)* beta,
    cuDoubleComplex* C,
    size_t ldc);

/* -------------------------------------------------------------------- */
/* HER2K : variant extension of HERK  */
cublasStatus_t cublasXtCher2k (
    cublasXtHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    size_t n,
    size_t k,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    size_t lda,
    const(cuComplex)* B,
    size_t ldb,
    const(float)* beta,
    cuComplex* C,
    size_t ldc);

cublasStatus_t cublasXtZher2k (
    cublasXtHandle_t handle,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    size_t n,
    size_t k,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    size_t lda,
    const(cuDoubleComplex)* B,
    size_t ldb,
    const(double)* beta,
    cuDoubleComplex* C,
    size_t ldc);

/* -------------------------------------------------------------------- */
/* SPMM : Symmetric Packed Multiply Matrix*/
cublasStatus_t cublasXtSspmm (
    cublasXtHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    size_t m,
    size_t n,
    const(float)* alpha,
    const(float)* AP,
    const(float)* B,
    size_t ldb,
    const(float)* beta,
    float* C,
    size_t ldc);

cublasStatus_t cublasXtDspmm (
    cublasXtHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    size_t m,
    size_t n,
    const(double)* alpha,
    const(double)* AP,
    const(double)* B,
    size_t ldb,
    const(double)* beta,
    double* C,
    size_t ldc);

cublasStatus_t cublasXtCspmm (
    cublasXtHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    size_t m,
    size_t n,
    const(cuComplex)* alpha,
    const(cuComplex)* AP,
    const(cuComplex)* B,
    size_t ldb,
    const(cuComplex)* beta,
    cuComplex* C,
    size_t ldc);

cublasStatus_t cublasXtZspmm (
    cublasXtHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    size_t m,
    size_t n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* AP,
    const(cuDoubleComplex)* B,
    size_t ldb,
    const(cuDoubleComplex)* beta,
    cuDoubleComplex* C,
    size_t ldc);

/* -------------------------------------------------------------------- */
/* TRMM */
cublasStatus_t cublasXtStrmm (
    cublasXtHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    size_t m,
    size_t n,
    const(float)* alpha,
    const(float)* A,
    size_t lda,
    const(float)* B,
    size_t ldb,
    float* C,
    size_t ldc);

cublasStatus_t cublasXtDtrmm (
    cublasXtHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    size_t m,
    size_t n,
    const(double)* alpha,
    const(double)* A,
    size_t lda,
    const(double)* B,
    size_t ldb,
    double* C,
    size_t ldc);

cublasStatus_t cublasXtCtrmm (
    cublasXtHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    size_t m,
    size_t n,
    const(cuComplex)* alpha,
    const(cuComplex)* A,
    size_t lda,
    const(cuComplex)* B,
    size_t ldb,
    cuComplex* C,
    size_t ldc);

cublasStatus_t cublasXtZtrmm (
    cublasXtHandle_t handle,
    cublasSideMode_t side,
    cublasFillMode_t uplo,
    cublasOperation_t trans,
    cublasDiagType_t diag,
    size_t m,
    size_t n,
    const(cuDoubleComplex)* alpha,
    const(cuDoubleComplex)* A,
    size_t lda,
    const(cuDoubleComplex)* B,
    size_t ldb,
    cuDoubleComplex* C,
    size_t ldc);
