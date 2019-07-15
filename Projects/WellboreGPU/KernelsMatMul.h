#include "pzreal.h"

__device__ void MultAddDevice (bool trans, int m, int n, int k, REAL *A, REAL *B, REAL *C, REAL alpha, REAL beta) {

	int aux1, aux2;

	if (trans == false) {
		aux1 = 1;
		aux2 = k;

	} else {
		aux1 = m;
		aux2 = 1;
	}

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			REAL c = C[j + i * n];
			REAL sum = 0.;
			for (int l = 0; l < k; l++) {
				sum += A[l * aux1 + i * aux2] * B[j + l * n];
			}
			C[j + i * n] = alpha * sum + beta * c; 
		}
	}
}

__global__ void MatrixMultiplicationKernel (bool trans, int *m, int *n, int *k, REAL *A, int *strideA, REAL *B, int *strideB, REAL *C, int *strideC, REAL alpha, int nmatrices) {

	int imatrix = blockIdx.x * blockDim.x + threadIdx.x;

	if (imatrix < nmatrices) {
		int m_i = m[imatrix];
		int n_i = n[imatrix];
		int k_i = k[imatrix];
		int strideA_i = strideA[imatrix]; 
		int strideB_i = strideB[imatrix]; 
		int strideC_i = strideC[imatrix]; 

		int beta = 0;

		MultAddDevice (trans, m_i, n_i, k_i, &A[strideA_i], &B[strideB_i], &C[strideC_i], alpha, beta);

		// 	//ROW MAJOR
		// for (int i = 0; i < m_i; i++) {
		// 	for (int j = 0; j < n_i; j++) {
		// 		C[j + i * n_i + strideC_i] = 0;
		// 		for (int l = 0; l < k_i; l++) {
		// 			C[j + i * n_i + strideC_i] += alpha * A[l * aux1 + i * aux2 + strideA_i] * B[j + l * n_i + strideB_i];
		// 		}
		// 	}
		// }

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