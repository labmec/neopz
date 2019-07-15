#include "pzreal.h"
#include "mkl.h"

#ifdef USING_TBB
#include "tbb/parallel_for.h"
#include "tbb/tick_count.h"
#endif

void MatrixMultiplication(bool trans, int *m, int *n, int *k, REAL *A, int *strideA, REAL *B, int *strideB, REAL *C, int *strideC, REAL alpha, int nmatrices) {

//    for (int imatrix = 0; imatrix < nmatrices; imatrix++) {
    
#ifdef USING_TBB
    tbb::parallel_for(size_t(0),size_t(nmatrices),size_t(1),[&](size_t imatrix)
#else
                      for (int imatrix = 0; imatrix < nmatrices; imatrix++)
#endif
                      {
                          int m_i = m[imatrix];
                          int n_i = n[imatrix];
                          int k_i = k[imatrix];
                          
                          int strideA_i = strideA[imatrix];
                          int strideB_i = strideB[imatrix];
                          int strideC_i = strideC[imatrix];
                          
                          //        int aux1, aux2;
                          int lda_i, ldb_i, ldc_i;
                          CBLAS_TRANSPOSE transpose;
                          if (trans == false) {
                              //            aux1 = 1;
                              //            aux2 = k_i;
                              
                              lda_i = k_i;
                              ldb_i = n_i;
                              ldc_i = n_i;
                              
                              transpose = CblasNoTrans;
                          } else {
                              //            aux1 = m_i;
                              //            aux2 = 1;
                              
                              lda_i = m_i;
                              ldb_i = n_i;
                              ldc_i = n_i;
                              
                              transpose = CblasTrans;
                          }
                          
                          //        for (<#initialization#>; <#condition#>; <#increment#>) {
                          //            <#statements#>
                          //        }
                          
                          cblas_dgemm(CblasRowMajor, transpose, CblasNoTrans, m_i, n_i, k_i, alpha, &A[strideA_i], lda_i,  &B[strideB_i], ldb_i, 0., &C[strideC_i], ldc_i);
                          
                          //        //            //ROW MAJOR
                          //        for (int i = 0; i < m_i; i++) {
                          //            for (int j = 0; j < n_i; j++) {
                          //                C[j + i * n_i + strideC_i] = 0;
                          //                for (int l = 0; l < k_i; l++) {
                          //                    C[j + i * n_i + strideC_i] += alpha * A[l * aux1 + i * aux2 + strideA_i] * B[j + l * n_i + strideB_i];
                          //                }
                          //            }
                          //        }
                          
                          //     //COL MAJOR
                          // for (int i = 0; i < m_i; i++) {
                          //     for (int j = 0; j < n_i; j++) {
                          //         C[j * m_i + i] = 0;
                          //         for (int l = 0; l < k_i; l++) {
                          //             C[j * m_i + i + strideC_i] += alpha * A[l * aux1 + i * aux2 + strideA_i] * B[j * k_i + l + strideB_i];
                          //         }
                          //     }
                          // }
                      }
#ifdef USING_TBB
                      );
#endif
}

void SpMV(bool trans, int m, int k, REAL alpha, REAL *csrVal, int *csrRowPtr, int *csrColInd, REAL *B, REAL *C) {

    char transpose;
    if(trans == false) {
        transpose = 'N';
    } else {
        transpose = 'T';

    }
    char matdescra[] = {'G',' ',' ','C'};
    REAL beta = 0.;
    mkl_dcsrmv(&transpose, &m, &k, &alpha, matdescra , csrVal, &csrColInd[0], &csrRowPtr[0], &csrRowPtr[1], B, &beta, C);
}
