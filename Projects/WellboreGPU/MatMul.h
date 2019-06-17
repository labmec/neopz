#include "pzreal.h"

void MatrixMultiplication(bool trans, int *m, int *n, int *k, REAL *A, int *strideA, REAL *B, int *strideB, REAL *C, int *strideC, REAL alpha, int nmatrices) {

	for (int imatrix = 0; imatrix < nmatrices; imatrix++) {
		int m_i = m[imatrix];
		int n_i = n[imatrix];
		int k_i = k[imatrix];
		int strideA_i = strideA[imatrix]; 
		int strideB_i = strideB[imatrix]; 
		int strideC_i = strideC[imatrix]; 

		int aux1, aux2;

		if (trans == false) {
			aux1 = 1;
			aux2 = k_i;

		} else {
			aux1 = m_i;
			aux2 = 1;
		}

			//ROW MAJOR
		for (int i = 0; i < m_i; i++) {
			for (int j = 0; j < n_i; j++) {
				C[j + i * n_i + strideC_i] = 0;
				for (int l = 0; l < k_i; l++) {
					C[j + i * n_i + strideC_i] += alpha * A[l * aux1 + i * aux2 + strideA_i] * B[j + l * n_i + strideB_i];
				}
			}
		}

		// 	//COL MAJOR
		// for (int i = 0; i < m_i; i++) {
		// 	for (int j = 0; j < n_i; j++) {
		// 		C[j * m_i + i] = 0;
		// 		for (int l = 0; l < k_i; l++) {
		// 			C[j * m_i + i + strideC_i] += alpha * A[l * aux1 + i * aux2 + strideA_i] * B[j * k_i + l + strideB_i];
		// 		}
		// 	}
		// }
	}
}