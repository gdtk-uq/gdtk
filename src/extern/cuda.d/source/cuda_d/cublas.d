/*
  D bindings for CUDA.
  Authors:    Prasun Anand
  Copyright:  Copyright (c) 2017, Prasun Anand. All rights reserved.
  License:    BSD 3-Clause License
*/

module cuda_d.cublas;

import cuda_d.cublas_api;
import cuda_d.cuComplex;

extern (C):

/* CUBLAS data types */
alias cublasStatus_t cublasStatus;

cublasStatus_t cublasInit ();
cublasStatus_t cublasShutdown ();
cublasStatus_t cublasGetError ();

cublasStatus_t cublasGetVersion (int* version_);
cublasStatus_t cublasAlloc (int n, int elemSize, void** devicePtr);

cublasStatus_t cublasFree (void* devicePtr);

cublasStatus_t cublasSetKernelStream (cudaStream_t stream);

/* ---------------- CUBLAS BLAS1 functions ---------------- */
/* NRM2 */
float cublasSnrm2 (int n, const(float)* x, int incx);
double cublasDnrm2 (int n, const(double)* x, int incx);
float cublasScnrm2 (int n, const(cuComplex)* x, int incx);
double cublasDznrm2 (int n, const(cuDoubleComplex)* x, int incx);
/*------------------------------------------------------------------------*/
/* DOT */
float cublasSdot (int n, const(float)* x, int incx, const(float)* y, int incy);
double cublasDdot (
    int n,
    const(double)* x,
    int incx,
    const(double)* y,
    int incy);
cuComplex cublasCdotu (
    int n,
    const(cuComplex)* x,
    int incx,
    const(cuComplex)* y,
    int incy);
cuComplex cublasCdotc (
    int n,
    const(cuComplex)* x,
    int incx,
    const(cuComplex)* y,
    int incy);
cuDoubleComplex cublasZdotu (
    int n,
    const(cuDoubleComplex)* x,
    int incx,
    const(cuDoubleComplex)* y,
    int incy);
cuDoubleComplex cublasZdotc (
    int n,
    const(cuDoubleComplex)* x,
    int incx,
    const(cuDoubleComplex)* y,
    int incy);

/*------------------------------------------------------------------------*/
/* SCAL */
void cublasSscal (int n, float alpha, float* x, int incx);
void cublasDscal (int n, double alpha, double* x, int incx);
void cublasCscal (int n, cuComplex alpha, cuComplex* x, int incx);
void cublasZscal (int n, cuDoubleComplex alpha, cuDoubleComplex* x, int incx);

void cublasCsscal (int n, float alpha, cuComplex* x, int incx);
void cublasZdscal (int n, double alpha, cuDoubleComplex* x, int incx);
/*------------------------------------------------------------------------*/
/* AXPY */
void cublasSaxpy (
    int n,
    float alpha,
    const(float)* x,
    int incx,
    float* y,
    int incy);
void cublasDaxpy (
    int n,
    double alpha,
    const(double)* x,
    int incx,
    double* y,
    int incy);
void cublasCaxpy (
    int n,
    cuComplex alpha,
    const(cuComplex)* x,
    int incx,
    cuComplex* y,
    int incy);
void cublasZaxpy (
    int n,
    cuDoubleComplex alpha,
    const(cuDoubleComplex)* x,
    int incx,
    cuDoubleComplex* y,
    int incy);

/*------------------------------------------------------------------------*/
/* COPY */
void cublasScopy (int n, const(float)* x, int incx, float* y, int incy);
void cublasDcopy (int n, const(double)* x, int incx, double* y, int incy);
void cublasCcopy (int n, const(cuComplex)* x, int incx, cuComplex* y, int incy);
void cublasZcopy (
    int n,
    const(cuDoubleComplex)* x,
    int incx,
    cuDoubleComplex* y,
    int incy);

/*------------------------------------------------------------------------*/
/* SWAP */
void cublasSswap (int n, float* x, int incx, float* y, int incy);
void cublasDswap (int n, double* x, int incx, double* y, int incy);
void cublasCswap (int n, cuComplex* x, int incx, cuComplex* y, int incy);
void cublasZswap (int n, cuDoubleComplex* x, int incx, cuDoubleComplex* y, int incy);
/*------------------------------------------------------------------------*/
/* AMAX */
int cublasIsamax (int n, const(float)* x, int incx);
int cublasIdamax (int n, const(double)* x, int incx);
int cublasIcamax (int n, const(cuComplex)* x, int incx);
int cublasIzamax (int n, const(cuDoubleComplex)* x, int incx);
/*------------------------------------------------------------------------*/
/* AMIN */
int cublasIsamin (int n, const(float)* x, int incx);
int cublasIdamin (int n, const(double)* x, int incx);

int cublasIcamin (int n, const(cuComplex)* x, int incx);
int cublasIzamin (int n, const(cuDoubleComplex)* x, int incx);
/*------------------------------------------------------------------------*/
/* ASUM */
float cublasSasum (int n, const(float)* x, int incx);
double cublasDasum (int n, const(double)* x, int incx);
float cublasScasum (int n, const(cuComplex)* x, int incx);
double cublasDzasum (int n, const(cuDoubleComplex)* x, int incx);
/*------------------------------------------------------------------------*/
/* ROT */
void cublasSrot (
    int n,
    float* x,
    int incx,
    float* y,
    int incy,
    float sc,
    float ss);
void cublasDrot (
    int n,
    double* x,
    int incx,
    double* y,
    int incy,
    double sc,
    double ss);
void cublasCrot (
    int n,
    cuComplex* x,
    int incx,
    cuComplex* y,
    int incy,
    float c,
    cuComplex s);
void cublasZrot (
    int n,
    cuDoubleComplex* x,
    int incx,
    cuDoubleComplex* y,
    int incy,
    double sc,
    cuDoubleComplex cs);
void cublasCsrot (
    int n,
    cuComplex* x,
    int incx,
    cuComplex* y,
    int incy,
    float c,
    float s);
void cublasZdrot (
    int n,
    cuDoubleComplex* x,
    int incx,
    cuDoubleComplex* y,
    int incy,
    double c,
    double s);

/*------------------------------------------------------------------------*/
/* ROTG */
void cublasSrotg (float* sa, float* sb, float* sc, float* ss);
void cublasDrotg (double* sa, double* sb, double* sc, double* ss);
void cublasCrotg (cuComplex* ca, cuComplex cb, float* sc, cuComplex* cs);
void cublasZrotg (
    cuDoubleComplex* ca,
    cuDoubleComplex cb,
    double* sc,
    cuDoubleComplex* cs);

/*------------------------------------------------------------------------*/
/* ROTM */
void cublasSrotm (
    int n,
    float* x,
    int incx,
    float* y,
    int incy,
    const(float)* sparam);
void cublasDrotm (
    int n,
    double* x,
    int incx,
    double* y,
    int incy,
    const(double)* sparam);

/*------------------------------------------------------------------------*/
/* ROTMG */
void cublasSrotmg (
    float* sd1,
    float* sd2,
    float* sx1,
    const(float)* sy1,
    float* sparam);
void cublasDrotmg (
    double* sd1,
    double* sd2,
    double* sx1,
    const(double)* sy1,
    double* sparam);

/* --------------- CUBLAS BLAS2 functions  ---------------- */
/* GEMV */
void cublasSgemv (
    char trans,
    int m,
    int n,
    float alpha,
    const(float)* A,
    int lda,
    const(float)* x,
    int incx,
    float beta,
    float* y,
    int incy);
void cublasDgemv (
    char trans,
    int m,
    int n,
    double alpha,
    const(double)* A,
    int lda,
    const(double)* x,
    int incx,
    double beta,
    double* y,
    int incy);
void cublasCgemv (
    char trans,
    int m,
    int n,
    cuComplex alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* x,
    int incx,
    cuComplex beta,
    cuComplex* y,
    int incy);
void cublasZgemv (
    char trans,
    int m,
    int n,
    cuDoubleComplex alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* x,
    int incx,
    cuDoubleComplex beta,
    cuDoubleComplex* y,
    int incy);

/*------------------------------------------------------------------------*/
/* GBMV */
void cublasSgbmv (
    char trans,
    int m,
    int n,
    int kl,
    int ku,
    float alpha,
    const(float)* A,
    int lda,
    const(float)* x,
    int incx,
    float beta,
    float* y,
    int incy);
void cublasDgbmv (
    char trans,
    int m,
    int n,
    int kl,
    int ku,
    double alpha,
    const(double)* A,
    int lda,
    const(double)* x,
    int incx,
    double beta,
    double* y,
    int incy);
void cublasCgbmv (
    char trans,
    int m,
    int n,
    int kl,
    int ku,
    cuComplex alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* x,
    int incx,
    cuComplex beta,
    cuComplex* y,
    int incy);
void cublasZgbmv (
    char trans,
    int m,
    int n,
    int kl,
    int ku,
    cuDoubleComplex alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* x,
    int incx,
    cuDoubleComplex beta,
    cuDoubleComplex* y,
    int incy);

/*------------------------------------------------------------------------*/
/* TRMV */
void cublasStrmv (
    char uplo,
    char trans,
    char diag,
    int n,
    const(float)* A,
    int lda,
    float* x,
    int incx);
void cublasDtrmv (
    char uplo,
    char trans,
    char diag,
    int n,
    const(double)* A,
    int lda,
    double* x,
    int incx);
void cublasCtrmv (
    char uplo,
    char trans,
    char diag,
    int n,
    const(cuComplex)* A,
    int lda,
    cuComplex* x,
    int incx);
void cublasZtrmv (
    char uplo,
    char trans,
    char diag,
    int n,
    const(cuDoubleComplex)* A,
    int lda,
    cuDoubleComplex* x,
    int incx);

/*------------------------------------------------------------------------*/
/* TBMV */
void cublasStbmv (
    char uplo,
    char trans,
    char diag,
    int n,
    int k,
    const(float)* A,
    int lda,
    float* x,
    int incx);
void cublasDtbmv (
    char uplo,
    char trans,
    char diag,
    int n,
    int k,
    const(double)* A,
    int lda,
    double* x,
    int incx);
void cublasCtbmv (
    char uplo,
    char trans,
    char diag,
    int n,
    int k,
    const(cuComplex)* A,
    int lda,
    cuComplex* x,
    int incx);
void cublasZtbmv (
    char uplo,
    char trans,
    char diag,
    int n,
    int k,
    const(cuDoubleComplex)* A,
    int lda,
    cuDoubleComplex* x,
    int incx);

/*------------------------------------------------------------------------*/
/* TPMV */
void cublasStpmv (char uplo, char trans, char diag, int n, const(float)* AP, float* x, int incx);

void cublasDtpmv (char uplo, char trans, char diag, int n, const(double)* AP, double* x, int incx);

void cublasCtpmv (char uplo, char trans, char diag, int n, const(cuComplex)* AP, cuComplex* x, int incx);

void cublasZtpmv (char uplo, char trans, char diag, int n, const(cuDoubleComplex)* AP, cuDoubleComplex* x, int incx);
/*------------------------------------------------------------------------*/
/* TRSV */
void cublasStrsv (char uplo, char trans, char diag, int n, const(float)* A, int lda, float* x, int incx);

void cublasDtrsv (char uplo, char trans, char diag, int n, const(double)* A, int lda, double* x, int incx);

void cublasCtrsv (char uplo, char trans, char diag, int n, const(cuComplex)* A, int lda, cuComplex* x, int incx);

void cublasZtrsv (
    char uplo,
    char trans,
    char diag,
    int n,
    const(cuDoubleComplex)* A,
    int lda,
    cuDoubleComplex* x,
    int incx);

/*------------------------------------------------------------------------*/
/* TPSV */
void cublasStpsv (
    char uplo,
    char trans,
    char diag,
    int n,
    const(float)* AP,
    float* x,
    int incx);

void cublasDtpsv (char uplo, char trans, char diag, int n, const(double)* AP, double* x, int incx);

void cublasCtpsv (char uplo, char trans, char diag, int n, const(cuComplex)* AP, cuComplex* x, int incx);

void cublasZtpsv (
    char uplo,
    char trans,
    char diag,
    int n,
    const(cuDoubleComplex)* AP,
    cuDoubleComplex* x,
    int incx);

/*------------------------------------------------------------------------*/
/* TBSV */
void cublasStbsv (
    char uplo,
    char trans,
    char diag,
    int n,
    int k,
    const(float)* A,
    int lda,
    float* x,
    int incx);

void cublasDtbsv (
    char uplo,
    char trans,
    char diag,
    int n,
    int k,
    const(double)* A,
    int lda,
    double* x,
    int incx);
void cublasCtbsv (
    char uplo,
    char trans,
    char diag,
    int n,
    int k,
    const(cuComplex)* A,
    int lda,
    cuComplex* x,
    int incx);

void cublasZtbsv (
    char uplo,
    char trans,
    char diag,
    int n,
    int k,
    const(cuDoubleComplex)* A,
    int lda,
    cuDoubleComplex* x,
    int incx);

/*------------------------------------------------------------------------*/
/* SYMV/HEMV */
void cublasSsymv (
    char uplo,
    int n,
    float alpha,
    const(float)* A,
    int lda,
    const(float)* x,
    int incx,
    float beta,
    float* y,
    int incy);
void cublasDsymv (
    char uplo,
    int n,
    double alpha,
    const(double)* A,
    int lda,
    const(double)* x,
    int incx,
    double beta,
    double* y,
    int incy);
void cublasChemv (
    char uplo,
    int n,
    cuComplex alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* x,
    int incx,
    cuComplex beta,
    cuComplex* y,
    int incy);
void cublasZhemv (
    char uplo,
    int n,
    cuDoubleComplex alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* x,
    int incx,
    cuDoubleComplex beta,
    cuDoubleComplex* y,
    int incy);

/*------------------------------------------------------------------------*/
/* SBMV/HBMV */
void cublasSsbmv (
    char uplo,
    int n,
    int k,
    float alpha,
    const(float)* A,
    int lda,
    const(float)* x,
    int incx,
    float beta,
    float* y,
    int incy);
void cublasDsbmv (
    char uplo,
    int n,
    int k,
    double alpha,
    const(double)* A,
    int lda,
    const(double)* x,
    int incx,
    double beta,
    double* y,
    int incy);
void cublasChbmv (
    char uplo,
    int n,
    int k,
    cuComplex alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* x,
    int incx,
    cuComplex beta,
    cuComplex* y,
    int incy);
void cublasZhbmv (
    char uplo,
    int n,
    int k,
    cuDoubleComplex alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* x,
    int incx,
    cuDoubleComplex beta,
    cuDoubleComplex* y,
    int incy);

/*------------------------------------------------------------------------*/
/* SPMV/HPMV */
void cublasSspmv (
    char uplo,
    int n,
    float alpha,
    const(float)* AP,
    const(float)* x,
    int incx,
    float beta,
    float* y,
    int incy);
void cublasDspmv (
    char uplo,
    int n,
    double alpha,
    const(double)* AP,
    const(double)* x,
    int incx,
    double beta,
    double* y,
    int incy);
void cublasChpmv (
    char uplo,
    int n,
    cuComplex alpha,
    const(cuComplex)* AP,
    const(cuComplex)* x,
    int incx,
    cuComplex beta,
    cuComplex* y,
    int incy);
void cublasZhpmv (
    char uplo,
    int n,
    cuDoubleComplex alpha,
    const(cuDoubleComplex)* AP,
    const(cuDoubleComplex)* x,
    int incx,
    cuDoubleComplex beta,
    cuDoubleComplex* y,
    int incy);

/*------------------------------------------------------------------------*/
/* GER */
void cublasSger (
    int m,
    int n,
    float alpha,
    const(float)* x,
    int incx,
    const(float)* y,
    int incy,
    float* A,
    int lda);
void cublasDger (
    int m,
    int n,
    double alpha,
    const(double)* x,
    int incx,
    const(double)* y,
    int incy,
    double* A,
    int lda);

void cublasCgeru (
    int m,
    int n,
    cuComplex alpha,
    const(cuComplex)* x,
    int incx,
    const(cuComplex)* y,
    int incy,
    cuComplex* A,
    int lda);
void cublasCgerc (
    int m,
    int n,
    cuComplex alpha,
    const(cuComplex)* x,
    int incx,
    const(cuComplex)* y,
    int incy,
    cuComplex* A,
    int lda);
void cublasZgeru (
    int m,
    int n,
    cuDoubleComplex alpha,
    const(cuDoubleComplex)* x,
    int incx,
    const(cuDoubleComplex)* y,
    int incy,
    cuDoubleComplex* A,
    int lda);
void cublasZgerc (
    int m,
    int n,
    cuDoubleComplex alpha,
    const(cuDoubleComplex)* x,
    int incx,
    const(cuDoubleComplex)* y,
    int incy,
    cuDoubleComplex* A,
    int lda);

/*------------------------------------------------------------------------*/
/* SYR/HER */
void cublasSsyr (
    char uplo,
    int n,
    float alpha,
    const(float)* x,
    int incx,
    float* A,
    int lda);
void cublasDsyr (
    char uplo,
    int n,
    double alpha,
    const(double)* x,
    int incx,
    double* A,
    int lda);

void cublasCher (
    char uplo,
    int n,
    float alpha,
    const(cuComplex)* x,
    int incx,
    cuComplex* A,
    int lda);
void cublasZher (
    char uplo,
    int n,
    double alpha,
    const(cuDoubleComplex)* x,
    int incx,
    cuDoubleComplex* A,
    int lda);

/*------------------------------------------------------------------------*/
/* SPR/HPR */
void cublasSspr (
    char uplo,
    int n,
    float alpha,
    const(float)* x,
    int incx,
    float* AP);
void cublasDspr (
    char uplo,
    int n,
    double alpha,
    const(double)* x,
    int incx,
    double* AP);
void cublasChpr (
    char uplo,
    int n,
    float alpha,
    const(cuComplex)* x,
    int incx,
    cuComplex* AP);
void cublasZhpr (
    char uplo,
    int n,
    double alpha,
    const(cuDoubleComplex)* x,
    int incx,
    cuDoubleComplex* AP);

/*------------------------------------------------------------------------*/
/* SYR2/HER2 */
void cublasSsyr2 (
    char uplo,
    int n,
    float alpha,
    const(float)* x,
    int incx,
    const(float)* y,
    int incy,
    float* A,
    int lda);
void cublasDsyr2 (
    char uplo,
    int n,
    double alpha,
    const(double)* x,
    int incx,
    const(double)* y,
    int incy,
    double* A,
    int lda);
void cublasCher2 (
    char uplo,
    int n,
    cuComplex alpha,
    const(cuComplex)* x,
    int incx,
    const(cuComplex)* y,
    int incy,
    cuComplex* A,
    int lda);
void cublasZher2 (
    char uplo,
    int n,
    cuDoubleComplex alpha,
    const(cuDoubleComplex)* x,
    int incx,
    const(cuDoubleComplex)* y,
    int incy,
    cuDoubleComplex* A,
    int lda);

/*------------------------------------------------------------------------*/
/* SPR2/HPR2 */
void cublasSspr2 (
    char uplo,
    int n,
    float alpha,
    const(float)* x,
    int incx,
    const(float)* y,
    int incy,
    float* AP);
void cublasDspr2 (
    char uplo,
    int n,
    double alpha,
    const(double)* x,
    int incx,
    const(double)* y,
    int incy,
    double* AP);
void cublasChpr2 (
    char uplo,
    int n,
    cuComplex alpha,
    const(cuComplex)* x,
    int incx,
    const(cuComplex)* y,
    int incy,
    cuComplex* AP);
void cublasZhpr2 (
    char uplo,
    int n,
    cuDoubleComplex alpha,
    const(cuDoubleComplex)* x,
    int incx,
    const(cuDoubleComplex)* y,
    int incy,
    cuDoubleComplex* AP);

/* ------------------------BLAS3 Functions ------------------------------- */
/* GEMM */
void cublasSgemm (
    char transa,
    char transb,
    int m,
    int n,
    int k,
    float alpha,
    const(float)* A,
    int lda,
    const(float)* B,
    int ldb,
    float beta,
    float* C,
    int ldc);
void cublasDgemm (
    char transa,
    char transb,
    int m,
    int n,
    int k,
    double alpha,
    const(double)* A,
    int lda,
    const(double)* B,
    int ldb,
    double beta,
    double* C,
    int ldc);
void cublasCgemm (
    char transa,
    char transb,
    int m,
    int n,
    int k,
    cuComplex alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* B,
    int ldb,
    cuComplex beta,
    cuComplex* C,
    int ldc);
void cublasZgemm (
    char transa,
    char transb,
    int m,
    int n,
    int k,
    cuDoubleComplex alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* B,
    int ldb,
    cuDoubleComplex beta,
    cuDoubleComplex* C,
    int ldc);

/* -------------------------------------------------------*/
/* SYRK */
void cublasSsyrk (
    char uplo,
    char trans,
    int n,
    int k,
    float alpha,
    const(float)* A,
    int lda,
    float beta,
    float* C,
    int ldc);
void cublasDsyrk (
    char uplo,
    char trans,
    int n,
    int k,
    double alpha,
    const(double)* A,
    int lda,
    double beta,
    double* C,
    int ldc);

void cublasCsyrk (
    char uplo,
    char trans,
    int n,
    int k,
    cuComplex alpha,
    const(cuComplex)* A,
    int lda,
    cuComplex beta,
    cuComplex* C,
    int ldc);
void cublasZsyrk (
    char uplo,
    char trans,
    int n,
    int k,
    cuDoubleComplex alpha,
    const(cuDoubleComplex)* A,
    int lda,
    cuDoubleComplex beta,
    cuDoubleComplex* C,
    int ldc);

/* ------------------------------------------------------- */
/* HERK */
void cublasCherk (
    char uplo,
    char trans,
    int n,
    int k,
    float alpha,
    const(cuComplex)* A,
    int lda,
    float beta,
    cuComplex* C,
    int ldc);
void cublasZherk (
    char uplo,
    char trans,
    int n,
    int k,
    double alpha,
    const(cuDoubleComplex)* A,
    int lda,
    double beta,
    cuDoubleComplex* C,
    int ldc);

/* ------------------------------------------------------- */
/* SYR2K */
void cublasSsyr2k (
    char uplo,
    char trans,
    int n,
    int k,
    float alpha,
    const(float)* A,
    int lda,
    const(float)* B,
    int ldb,
    float beta,
    float* C,
    int ldc);

void cublasDsyr2k (
    char uplo,
    char trans,
    int n,
    int k,
    double alpha,
    const(double)* A,
    int lda,
    const(double)* B,
    int ldb,
    double beta,
    double* C,
    int ldc);
void cublasCsyr2k (
    char uplo,
    char trans,
    int n,
    int k,
    cuComplex alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* B,
    int ldb,
    cuComplex beta,
    cuComplex* C,
    int ldc);

void cublasZsyr2k (
    char uplo,
    char trans,
    int n,
    int k,
    cuDoubleComplex alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* B,
    int ldb,
    cuDoubleComplex beta,
    cuDoubleComplex* C,
    int ldc);

/* ------------------------------------------------------- */
/* HER2K */
void cublasCher2k (
    char uplo,
    char trans,
    int n,
    int k,
    cuComplex alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* B,
    int ldb,
    float beta,
    cuComplex* C,
    int ldc);

void cublasZher2k (
    char uplo,
    char trans,
    int n,
    int k,
    cuDoubleComplex alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* B,
    int ldb,
    double beta,
    cuDoubleComplex* C,
    int ldc);

/*------------------------------------------------------------------------*/
/* SYMM*/
void cublasSsymm (
    char side,
    char uplo,
    int m,
    int n,
    float alpha,
    const(float)* A,
    int lda,
    const(float)* B,
    int ldb,
    float beta,
    float* C,
    int ldc);
void cublasDsymm (
    char side,
    char uplo,
    int m,
    int n,
    double alpha,
    const(double)* A,
    int lda,
    const(double)* B,
    int ldb,
    double beta,
    double* C,
    int ldc);

void cublasCsymm (
    char side,
    char uplo,
    int m,
    int n,
    cuComplex alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* B,
    int ldb,
    cuComplex beta,
    cuComplex* C,
    int ldc);

void cublasZsymm (
    char side,
    char uplo,
    int m,
    int n,
    cuDoubleComplex alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* B,
    int ldb,
    cuDoubleComplex beta,
    cuDoubleComplex* C,
    int ldc);

/*------------------------------------------------------------------------*/
/* HEMM*/
void cublasChemm (
    char side,
    char uplo,
    int m,
    int n,
    cuComplex alpha,
    const(cuComplex)* A,
    int lda,
    const(cuComplex)* B,
    int ldb,
    cuComplex beta,
    cuComplex* C,
    int ldc);
void cublasZhemm (
    char side,
    char uplo,
    int m,
    int n,
    cuDoubleComplex alpha,
    const(cuDoubleComplex)* A,
    int lda,
    const(cuDoubleComplex)* B,
    int ldb,
    cuDoubleComplex beta,
    cuDoubleComplex* C,
    int ldc);

/*------------------------------------------------------------------------*/
/* TRSM*/
void cublasStrsm (
    char side,
    char uplo,
    char transa,
    char diag,
    int m,
    int n,
    float alpha,
    const(float)* A,
    int lda,
    float* B,
    int ldb);

void cublasDtrsm (
    char side,
    char uplo,
    char transa,
    char diag,
    int m,
    int n,
    double alpha,
    const(double)* A,
    int lda,
    double* B,
    int ldb);

void cublasCtrsm (
    char side,
    char uplo,
    char transa,
    char diag,
    int m,
    int n,
    cuComplex alpha,
    const(cuComplex)* A,
    int lda,
    cuComplex* B,
    int ldb);

void cublasZtrsm (
    char side,
    char uplo,
    char transa,
    char diag,
    int m,
    int n,
    cuDoubleComplex alpha,
    const(cuDoubleComplex)* A,
    int lda,
    cuDoubleComplex* B,
    int ldb);

/*------------------------------------------------------------------------*/
/* TRMM*/
void cublasStrmm (
    char side,
    char uplo,
    char transa,
    char diag,
    int m,
    int n,
    float alpha,
    const(float)* A,
    int lda,
    float* B,
    int ldb);
void cublasDtrmm (
    char side,
    char uplo,
    char transa,
    char diag,
    int m,
    int n,
    double alpha,
    const(double)* A,
    int lda,
    double* B,
    int ldb);
void cublasCtrmm (
    char side,
    char uplo,
    char transa,
    char diag,
    int m,
    int n,
    cuComplex alpha,
    const(cuComplex)* A,
    int lda,
    cuComplex* B,
    int ldb);
void cublasZtrmm (
    char side,
    char uplo,
    char transa,
    char diag,
    int m,
    int n,
    cuDoubleComplex alpha,
    const(cuDoubleComplex)* A,
    int lda,
    cuDoubleComplex* B,
    int ldb);
