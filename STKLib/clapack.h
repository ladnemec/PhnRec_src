#ifndef CLAPACK_H

#define CLAPACK_H
#include "cblas.h"

#ifndef ATL_INT
   #define ATL_INT int
#endif
#ifndef ATL_CINT
   #define ATL_CINT const ATL_INT
#endif
#ifndef ATLAS_ORDER
   #define ATLAS_ORDER CBLAS_ORDER
#endif
#ifndef ATLAS_UPLO
   #define ATLAS_UPLO CBLAS_UPLO
#endif
#ifndef ATLAS_DIAG
   #define ATLAS_DIAG CBLAS_DIAG
#endif
int clapack_sgesv(const  CBLAS_ORDER Order, const int N, const int NRHS,
                  float *A, const int lda, int *ipiv,
                  float *B, const int ldb);
int clapack_sgetrf(const  CBLAS_ORDER Order, const int M, const int N,
                   float *A, const int lda, int *ipiv);
int clapack_sgetrs
   (const  CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
    const int N, const int NRHS, const float *A, const int lda,
    const int *ipiv, float *B, const int ldb);
int clapack_sgetri(const  CBLAS_ORDER Order, const int N, float *A,
                   const int lda, const int *ipiv);
int clapack_sposv(const  ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                  const int N, const int NRHS, float *A, const int lda,
                  float *B, const int ldb);
int clapack_spotrf(const  ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                   const int N, float *A, const int lda);
int clapack_spotrs(const  CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, const int NRHS, const float *A, const int lda,
                   float *B, const int ldb);
int clapack_spotri(const  ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                   const int N, float *A, const int lda);
int clapack_slauum(const  ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                   const int N, float *A, const int lda);
int clapack_strtri(const  ATLAS_ORDER Order,const enum ATLAS_UPLO Uplo,
                   const  ATLAS_DIAG Diag, const int N, float *A,
                   const int lda);
int clapack_sgels(const  CBLAS_ORDER Order,
                  const  CBLAS_TRANSPOSE TA,
                  ATL_CINT M, ATL_CINT N, ATL_CINT NRHS, float *A,
                  ATL_CINT lda, float *B, const int ldb);
int clapack_sgelqf(const  CBLAS_ORDER Order, ATL_CINT M, ATL_CINT N,
                   float *A, ATL_CINT lda, float *TAU);
int clapack_sgeqlf(const  CBLAS_ORDER Order, ATL_CINT M, ATL_CINT N,
                   float *A, ATL_CINT lda, float *TAU);
int clapack_sgerqf(const  CBLAS_ORDER Order, ATL_CINT M, ATL_CINT N,
                   float *A, ATL_CINT lda, float *TAU);
int clapack_sgeqrf(const  CBLAS_ORDER Order, ATL_CINT M, ATL_CINT N,
                   float *A, ATL_CINT lda, float *TAU);

int clapack_dgesv(const  CBLAS_ORDER Order, const int N, const int NRHS,
                  double *A, const int lda, int *ipiv,
                  double *B, const int ldb);
int clapack_dgetrf(const  CBLAS_ORDER Order, const int M, const int N,
                   double *A, const int lda, int *ipiv);
int clapack_dgetrs
   (const  CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
    const int N, const int NRHS, const double *A, const int lda,
    const int *ipiv, double *B, const int ldb);
int clapack_dgetri(const  CBLAS_ORDER Order, const int N, double *A,
                   const int lda, const int *ipiv);
int clapack_dposv(const  ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                  const int N, const int NRHS, double *A, const int lda,
                  double *B, const int ldb);
int clapack_dpotrf(const  ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                   const int N, double *A, const int lda);
int clapack_dpotrs(const  CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, const int NRHS, const double *A, const int lda,
                   double *B, const int ldb);
int clapack_dpotri(const  ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                   const int N, double *A, const int lda);
int clapack_dlauum(const  ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                   const int N, double *A, const int lda);
int clapack_dtrtri(const  ATLAS_ORDER Order,const enum ATLAS_UPLO Uplo,
                   const  ATLAS_DIAG Diag, const int N, double *A,
                   const int lda);
int clapack_dgels(const  CBLAS_ORDER Order,
                  const  CBLAS_TRANSPOSE TA,
                  ATL_CINT M, ATL_CINT N, ATL_CINT NRHS, double *A,
                  ATL_CINT lda, double *B, const int ldb);
int clapack_dgelqf(const  CBLAS_ORDER Order, ATL_CINT M, ATL_CINT N,
                   double *A, ATL_CINT lda, double *TAU);
int clapack_dgeqlf(const  CBLAS_ORDER Order, ATL_CINT M, ATL_CINT N,
                   double *A, ATL_CINT lda, double *TAU);
int clapack_dgerqf(const  CBLAS_ORDER Order, ATL_CINT M, ATL_CINT N,
                   double *A, ATL_CINT lda, double *TAU);
int clapack_dgeqrf(const  CBLAS_ORDER Order, ATL_CINT M, ATL_CINT N,
                   double *A, ATL_CINT lda, double *TAU);

int clapack_cgesv(const  CBLAS_ORDER Order, const int N, const int NRHS,
                  void *A, const int lda, int *ipiv,
                  void *B, const int ldb);
int clapack_cgetrf(const  CBLAS_ORDER Order, const int M, const int N,
                   void *A, const int lda, int *ipiv);
int clapack_cgetrs
   (const  CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
    const int N, const int NRHS, const void *A, const int lda,
    const int *ipiv, void *B, const int ldb);
int clapack_cgetri(const  CBLAS_ORDER Order, const int N, void *A,
                   const int lda, const int *ipiv);
int clapack_cposv(const  ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                  const int N, const int NRHS, void *A, const int lda,
                  void *B, const int ldb);
int clapack_cpotrf(const  ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                   const int N, void *A, const int lda);
int clapack_cpotrs(const  CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, const int NRHS, const void *A, const int lda,
                   void *B, const int ldb);
int clapack_cpotri(const  ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                   const int N, void *A, const int lda);
int clapack_clauum(const  ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                   const int N, void *A, const int lda);
int clapack_ctrtri(const  ATLAS_ORDER Order,const enum ATLAS_UPLO Uplo,
                   const  ATLAS_DIAG Diag, const int N, void *A,
                   const int lda);
int clapack_cgels(const  CBLAS_ORDER Order,
                  const  CBLAS_TRANSPOSE TA,
                  ATL_CINT M, ATL_CINT N, ATL_CINT NRHS, void *A,
                  ATL_CINT lda, void *B, const int ldb);
int clapack_cgelqf(const  CBLAS_ORDER Order, ATL_CINT M, ATL_CINT N,
                   void *A, ATL_CINT lda, void *TAU);
int clapack_cgeqlf(const  CBLAS_ORDER Order, ATL_CINT M, ATL_CINT N,
                   void *A, ATL_CINT lda, void *TAU);
int clapack_cgerqf(const  CBLAS_ORDER Order, ATL_CINT M, ATL_CINT N,
                   void *A, ATL_CINT lda, void *TAU);
int clapack_cgeqrf(const  CBLAS_ORDER Order, ATL_CINT M, ATL_CINT N,
                   void *A, ATL_CINT lda, void *TAU);

int clapack_zgesv(const  CBLAS_ORDER Order, const int N, const int NRHS,
                  void *A, const int lda, int *ipiv,
                  void *B, const int ldb);
int clapack_zgetrf(const  CBLAS_ORDER Order, const int M, const int N,
                   void *A, const int lda, int *ipiv);
int clapack_zgetrs
   (const  CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE Trans,
    const int N, const int NRHS, const void *A, const int lda,
    const int *ipiv, void *B, const int ldb);
int clapack_zgetri(const  CBLAS_ORDER Order, const int N, void *A,
                   const int lda, const int *ipiv);
int clapack_zposv(const  ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                  const int N, const int NRHS, void *A, const int lda,
                  void *B, const int ldb);
int clapack_zpotrf(const  ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                   const int N, void *A, const int lda);
int clapack_zpotrs(const  CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                   const int N, const int NRHS, const void *A, const int lda,
                   void *B, const int ldb);
int clapack_zpotri(const  ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                   const int N, void *A, const int lda);
int clapack_zlauum(const  ATLAS_ORDER Order, const enum ATLAS_UPLO Uplo,
                   const int N, void *A, const int lda);
int clapack_ztrtri(const  ATLAS_ORDER Order,const enum ATLAS_UPLO Uplo,
                   const  ATLAS_DIAG Diag, const int N, void *A,
                   const int lda);
int clapack_zgels(const  CBLAS_ORDER Order,
                  const  CBLAS_TRANSPOSE TA,
                  ATL_CINT M, ATL_CINT N, ATL_CINT NRHS, void *A,
                  ATL_CINT lda, void *B, const int ldb);
int clapack_zgelqf(const  CBLAS_ORDER Order, ATL_CINT M, ATL_CINT N,
                   void *A, ATL_CINT lda, void *TAU);
int clapack_zgeqlf(const  CBLAS_ORDER Order, ATL_CINT M, ATL_CINT N,
                   void *A, ATL_CINT lda, void *TAU);
int clapack_zgerqf(const  CBLAS_ORDER Order, ATL_CINT M, ATL_CINT N,
                   void *A, ATL_CINT lda, void *TAU);
int clapack_zgeqrf(const  CBLAS_ORDER Order, ATL_CINT M, ATL_CINT N,
                   void *A, ATL_CINT lda, void *TAU);

#endif
