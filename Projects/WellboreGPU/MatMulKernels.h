#include "pzreal.h"

extern "C" {
__global__ void MatMulcuBLASKernel(cublasOperation_t trans, int64_t nelem,
		REAL *A, int *rowsizes, int *colsizes, int *matrixpos,
		int *rowfirstindex, int* colfirstindex, int npts, int nphis, REAL *B,
		REAL *C) {

	int iel = blockIdx.x * blockDim.x + threadIdx.x;

	REAL alpha;
	REAL beta;

	int lda, ldb, ldc;
	int Bpos, Cpos;
	int Boffset, Coffset;
	int m, n, k;
	int Apos;

	if (iel < nelem) {
		cublasHandle_t cnpHandle; //each thread must have its own handle
		cublasCreate(&cnpHandle);

		Apos = matrixpos[iel];

		if (trans == CUBLAS_OP_N) {
			m = rowsizes[iel];
			n = 1;
			k = colsizes[iel];

			alpha = 1.;
			beta = 0;

			lda = m;
			ldb = k;
			ldc = m;

			Bpos = colfirstindex[iel];
			Boffset = nphis;

			Cpos = rowfirstindex[iel];
			Coffset = npts;

		} else if (trans == CUBLAS_OP_T) {
			m = colsizes[iel];
			n = 1;
			k = rowsizes[iel];

			alpha = -1.;
			beta = 0;

			lda = k;
			ldb = k;
			ldc = m;

			Bpos = rowfirstindex[iel];
			Boffset = npts;

			Cpos = colfirstindex[iel];
			Coffset = nphis;
		}
		cublasDgemm(cnpHandle, trans, CUBLAS_OP_N, m, n, k, &alpha, &A[Apos],
				lda, &B[Bpos], ldb, &beta, &C[Cpos], ldc);

		cublasDgemm(cnpHandle, trans, CUBLAS_OP_N, m, n, k, &alpha, &A[Apos],
				lda, &B[Bpos + Boffset], ldb, &beta, &C[Cpos + Coffset], ldc);

		__syncthreads();
		cublasDestroy(cnpHandle);

	}
}
}

__global__ void MatMulKernel(bool trans, int64_t nelem, REAL *A, int *rowsizes,
		int *colsizes, int *matrixpos, int *rowfirstindex, int* colfirstindex,
		int npts, int nphis, REAL *B, REAL *C) {
	int iel = blockIdx.x * blockDim.x + threadIdx.x;

	__shared__ REAL alpha;

	int Bpos, Cpos;
	int Boffset, Coffset;
	int m, k;
	int Apos;
	int aux1;
	int aux2;

	if (iel < nelem) {
		Apos = matrixpos[iel];

		if (trans == false) {
			m = rowsizes[iel];
			k = colsizes[iel];

			aux1 = rowsizes[iel];
			aux2 = 1;

			alpha = 1.;

			Bpos = colfirstindex[iel];
			Boffset = nphis;

			Cpos = rowfirstindex[iel];
			Coffset = npts;

		} else if (trans == true) {
			m = colsizes[iel];
			k = rowsizes[iel];

			aux1 = 1;
			aux2 = rowsizes[iel];

			alpha = -1.;

			Bpos = rowfirstindex[iel];
			Boffset = npts;

			Cpos = colfirstindex[iel];
			Coffset = nphis;
		}
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < k; j++) {
				C[i + Cpos] += alpha * A[j * aux1 + i * aux2 + Apos] * B[j + Bpos];
				C[i + Cpos + Coffset] += alpha * A[j * aux1 + i * aux2 + Apos] * B[j + Bpos + Boffset];
			}
		}
	}
}
