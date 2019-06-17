#include "pzreal.h"

//#define TILE_WIDTH 16

//extern "C" {
//__global__ void matrixMultiply(REAL * A, REAL * B, REAL * C,
//  		       int numARows, int numAColumns,
//			       int numBRows, int numBColumns,
//			       int numCRows, int numCColumns) {
//    //@@ Insert code to implement matrix multiplication here
//    __shared__ REAL ds_M[TILE_WIDTH][TILE_WIDTH];
//    __shared__ REAL ds_N[TILE_WIDTH][TILE_WIDTH];
//    int bx = blockIdx.x, by = blockIdx.y,
//       tx = threadIdx.x, ty = threadIdx.y,
//       Row = by * TILE_WIDTH + ty,
//       Col = bx * TILE_WIDTH + tx;
//    REAL Pvalue = 0;
//
//    for (int m = 0; m < (numAColumns-1)/TILE_WIDTH+1; ++m) {
//       if (Row < numARows && m*TILE_WIDTH+tx < numAColumns)
//          ds_M[ty][tx] = A[Row + m*TILE_WIDTH+tx*numARows];
//       else
//          ds_M[ty][tx] = 0;
//       if (Col < numBColumns && m*TILE_WIDTH+ty < numBRows)
//          ds_N[ty][tx] = B[(m*TILE_WIDTH+ty)*numBColumns+Col];
//       else
//          ds_N[ty][tx] = 0;
//
//       __syncthreads();
//       for (int k = 0; k < TILE_WIDTH; ++k)
//          Pvalue += ds_M[ty][k] * ds_N[k][tx];
//       __syncthreads();
//    }
//    if (Row < numCRows && Col < numCColumns)
//       C[Row*numCColumns+Col] = Pvalue;
//}
//}

__global__ void MatrixMultiplicationKernel (bool trans, int *m, int *n, int *k, REAL *A, int *strideA, REAL *B, int *strideB, REAL *C, int *strideC, REAL alpha, int nmatrices) {

	int imatrix = blockIdx.x;

	if (imatrix < nmatrices) {
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

// extern "C" {
// __global__ void MatMulcuBLASKernel(cublasOperation_t trans, int64_t nelem,
// 		REAL *A, int *rowsizes, int *colsizes, int *matrixpos,
// 		int *rowfirstindex, int* colfirstindex, int npts, int nphis, REAL *B,
// 		REAL *C) {

// 	int iel = blockIdx.x * blockDim.x + threadIdx.x;

// 	REAL alpha;
// 	REAL beta;

// 	int lda, ldb, ldc;
// 	int Bpos, Cpos;
// 	int Boffset, Coffset;
// 	int m, n, k;
// 	int Apos;

// 	if (iel < nelem) {
// 		cublasHandle_t cnpHandle; //each thread must have its own handle
// 		cublasCreate(&cnpHandle);

// 		Apos = matrixpos[iel];

// 		if (trans == CUBLAS_OP_N) {
// 			m = rowsizes[iel];
// 			n = 1;
// 			k = colsizes[iel];

// 			alpha = 1.;
// 			beta = 0;

// 			lda = m;
// 			ldb = k;
// 			ldc = m;

// 			Bpos = colfirstindex[iel];
// 			Boffset = nphis;
// 			Cpos = rowfirstindex[iel];
// 			Coffset = npts;

// 		} else if (trans == CUBLAS_OP_T) {
// 			m = colsizes[iel];
// 			n = 1;
// 			k = rowsizes[iel];

// 			alpha = -1.;
// 			beta = 0;

// 			lda = k;
// 			ldb = k;
// 			ldc = m;

// 			Bpos = rowfirstindex[iel];
// 			Boffset = npts;
// 			Cpos = colfirstindex[iel];
// 			Coffset = nphis;
// 		}
// 		cublasDgemm(cnpHandle, trans, CUBLAS_OP_N, m, n, k, &alpha, &A[Apos],
// 				lda, &B[Bpos], ldb, &beta, &C[Cpos], ldc);

//// 		cublasDgemm(cnpHandle, trans, CUBLAS_OP_N, m, n, k, &alpha, &A[Apos],
// //				lda, &B[Bpos + Boffset], ldb, &beta, &C[Cpos + Coffset], ldc);

// 		__syncthreads();
// 		cublasDestroy(cnpHandle);
//		}
//}
//}



