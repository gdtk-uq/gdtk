/*
  D bindings for CUDA.
  Authors:    Prasun Anand
  Copyright:  Copyright (c) 2017, Prasun Anand. All rights reserved.
  License:    BSD 3-Clause License
*/

module cuda_d.cublas_v2;

import cuda_d.cublas_api;
import cuda_d.cuComplex;

alias cublasCreate = cublasCreate_v2;
alias cublasDestroy = cublasDestroy_v2;
alias cublasGetVersion = cublasGetVersion_v2;
alias cublasSetStream = cublasSetStream_v2;
alias cublasGetStream = cublasGetStream_v2;
alias cublasGetPointerMode = cublasGetPointerMode_v2;
alias cublasSetPointerMode = cublasSetPointerMode_v2;

/* Blas3 Routines   */

alias cublasSnrm2 = cublasSnrm2_v2;
alias cublasDnrm2 = cublasDnrm2_v2;
alias cublasScnrm2 = cublasScnrm2_v2;
alias cublasDznrm2 = cublasDznrm2_v2;

alias cublasSdot = cublasSdot_v2;
alias cublasDdot = cublasDdot_v2;
alias cublasCdotu = cublasCdotu_v2;
alias cublasCdotc = cublasCdotc_v2;
alias cublasZdotu = cublasZdotu_v2;
alias cublasZdotc = cublasZdotc_v2;

alias cublasSscal = cublasSscal_v2;
alias cublasDscal = cublasDscal_v2;
alias cublasCscal = cublasCscal_v2;
alias cublasCsscal = cublasCsscal_v2;
alias cublasZscal = cublasZscal_v2;
alias cublasZdscal = cublasZdscal_v2;

alias cublasSaxpy = cublasSaxpy_v2;
alias cublasDaxpy = cublasDaxpy_v2;
alias cublasCaxpy = cublasCaxpy_v2;
alias cublasZaxpy = cublasZaxpy_v2;

alias cublasScopy = cublasScopy_v2;
alias cublasDcopy = cublasDcopy_v2;
alias cublasCcopy = cublasCcopy_v2;
alias cublasZcopy = cublasZcopy_v2;

alias cublasSswap = cublasSswap_v2;
alias cublasDswap = cublasDswap_v2;
alias cublasCswap = cublasCswap_v2;
alias cublasZswap = cublasZswap_v2;

alias cublasIsamax = cublasIsamax_v2;
alias cublasIdamax = cublasIdamax_v2;
alias cublasIcamax = cublasIcamax_v2;
alias cublasIzamax = cublasIzamax_v2;

alias cublasIsamin = cublasIsamin_v2;
alias cublasIdamin = cublasIdamin_v2;
alias cublasIcamin = cublasIcamin_v2;
alias cublasIzamin = cublasIzamin_v2;

alias cublasSasum = cublasSasum_v2;
alias cublasDasum = cublasDasum_v2;
alias cublasScasum = cublasScasum_v2;
alias cublasDzasum = cublasDzasum_v2;

alias cublasSrot = cublasSrot_v2;
alias cublasDrot = cublasDrot_v2;
alias cublasCrot = cublasCrot_v2;
alias cublasCsrot = cublasCsrot_v2;
alias cublasZrot = cublasZrot_v2;
alias cublasZdrot = cublasZdrot_v2;

alias cublasSrotg = cublasSrotg_v2;
alias cublasDrotg = cublasDrotg_v2;
alias cublasCrotg = cublasCrotg_v2;
alias cublasZrotg = cublasZrotg_v2;

alias cublasSrotm = cublasSrotm_v2;
alias cublasDrotm = cublasDrotm_v2;

alias cublasSrotmg = cublasSrotmg_v2;
alias cublasDrotmg = cublasDrotmg_v2;

/* Blas2 Routines */

alias cublasSgemv = cublasSgemv_v2;
alias cublasDgemv = cublasDgemv_v2;
alias cublasCgemv = cublasCgemv_v2;
alias cublasZgemv = cublasZgemv_v2;

alias cublasSgbmv = cublasSgbmv_v2;
alias cublasDgbmv = cublasDgbmv_v2;
alias cublasCgbmv = cublasCgbmv_v2;
alias cublasZgbmv = cublasZgbmv_v2;

alias cublasStrmv = cublasStrmv_v2;
alias cublasDtrmv = cublasDtrmv_v2;
alias cublasCtrmv = cublasCtrmv_v2;
alias cublasZtrmv = cublasZtrmv_v2;

alias cublasStbmv = cublasStbmv_v2;
alias cublasDtbmv = cublasDtbmv_v2;
alias cublasCtbmv = cublasCtbmv_v2;
alias cublasZtbmv = cublasZtbmv_v2;

alias cublasStpmv = cublasStpmv_v2;
alias cublasDtpmv = cublasDtpmv_v2;
alias cublasCtpmv = cublasCtpmv_v2;
alias cublasZtpmv = cublasZtpmv_v2;

alias cublasStrsv = cublasStrsv_v2;
alias cublasDtrsv = cublasDtrsv_v2;
alias cublasCtrsv = cublasCtrsv_v2;
alias cublasZtrsv = cublasZtrsv_v2;

alias cublasStpsv = cublasStpsv_v2;
alias cublasDtpsv = cublasDtpsv_v2;
alias cublasCtpsv = cublasCtpsv_v2;
alias cublasZtpsv = cublasZtpsv_v2;

alias cublasStbsv = cublasStbsv_v2;
alias cublasDtbsv = cublasDtbsv_v2;
alias cublasCtbsv = cublasCtbsv_v2;
alias cublasZtbsv = cublasZtbsv_v2;

alias cublasSsymv = cublasSsymv_v2;
alias cublasDsymv = cublasDsymv_v2;
alias cublasCsymv = cublasCsymv_v2;
alias cublasZsymv = cublasZsymv_v2;
alias cublasChemv = cublasChemv_v2;
alias cublasZhemv = cublasZhemv_v2;

alias cublasSsbmv = cublasSsbmv_v2;
alias cublasDsbmv = cublasDsbmv_v2;
alias cublasChbmv = cublasChbmv_v2;
alias cublasZhbmv = cublasZhbmv_v2;

alias cublasSspmv = cublasSspmv_v2;
alias cublasDspmv = cublasDspmv_v2;
alias cublasChpmv = cublasChpmv_v2;
alias cublasZhpmv = cublasZhpmv_v2;

alias cublasSger = cublasSger_v2;
alias cublasDger = cublasDger_v2;
alias cublasCgeru = cublasCgeru_v2;
alias cublasCgerc = cublasCgerc_v2;
alias cublasZgeru = cublasZgeru_v2;
alias cublasZgerc = cublasZgerc_v2;

alias cublasSsyr = cublasSsyr_v2;
alias cublasDsyr = cublasDsyr_v2;
alias cublasCsyr = cublasCsyr_v2;
alias cublasZsyr = cublasZsyr_v2;
alias cublasCher = cublasCher_v2;
alias cublasZher = cublasZher_v2;

alias cublasSspr = cublasSspr_v2;
alias cublasDspr = cublasDspr_v2;
alias cublasChpr = cublasChpr_v2;
alias cublasZhpr = cublasZhpr_v2;

alias cublasSsyr2 = cublasSsyr2_v2;
alias cublasDsyr2 = cublasDsyr2_v2;
alias cublasCsyr2 = cublasCsyr2_v2;
alias cublasZsyr2 = cublasZsyr2_v2;
alias cublasCher2 = cublasCher2_v2;
alias cublasZher2 = cublasZher2_v2;

alias cublasSspr2 = cublasSspr2_v2;
alias cublasDspr2 = cublasDspr2_v2;
alias cublasChpr2 = cublasChpr2_v2;
alias cublasZhpr2 = cublasZhpr2_v2;

/* Blas3 Routines   */

alias cublasSgemm = cublasSgemm_v2;
alias cublasDgemm = cublasDgemm_v2;
alias cublasCgemm = cublasCgemm_v2;
alias cublasZgemm = cublasZgemm_v2;

alias cublasSsyrk = cublasSsyrk_v2;
alias cublasDsyrk = cublasDsyrk_v2;
alias cublasCsyrk = cublasCsyrk_v2;
alias cublasZsyrk = cublasZsyrk_v2;
alias cublasCherk = cublasCherk_v2;
alias cublasZherk = cublasZherk_v2;

alias cublasSsyr2k = cublasSsyr2k_v2;
alias cublasDsyr2k = cublasDsyr2k_v2;
alias cublasCsyr2k = cublasCsyr2k_v2;
alias cublasZsyr2k = cublasZsyr2k_v2;
alias cublasCher2k = cublasCher2k_v2;
alias cublasZher2k = cublasZher2k_v2;

alias cublasSsymm = cublasSsymm_v2;
alias cublasDsymm = cublasDsymm_v2;
alias cublasCsymm = cublasCsymm_v2;
alias cublasZsymm = cublasZsymm_v2;
alias cublasChemm = cublasChemm_v2;
alias cublasZhemm = cublasZhemm_v2;

alias cublasStrsm = cublasStrsm_v2;
alias cublasDtrsm = cublasDtrsm_v2;
alias cublasCtrsm = cublasCtrsm_v2;
alias cublasZtrsm = cublasZtrsm_v2;

alias cublasStrmm = cublasStrmm_v2;
alias cublasDtrmm = cublasDtrmm_v2;
alias cublasCtrmm = cublasCtrmm_v2;
alias cublasZtrmm = cublasZtrmm_v2;

/* !defined(CUBLAS_V2_H_) */
