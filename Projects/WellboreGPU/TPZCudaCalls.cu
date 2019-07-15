#include "TPZCudaCalls.h"
#include "pzreal.h"
#include "pzvec.h"

// #include "MatMulKernels.h"
#include "KernelsComputeSigma.h"
#include "KernelsMatMul.h"
#include "KernelsMatrixAssemble.h"

#define NT  64

TPZCudaCalls::TPZCudaCalls() {
	cusparse_h = false;
	cublas_h = false;
}

TPZCudaCalls::~TPZCudaCalls() {
	if(cublas_h == true) {
		cublasDestroy(handle_cublas);
	}
	if(cusparse_h == true) {
		cusparseDestroy(handle_cusparse);			
	}
}

TPZCudaCalls &TPZCudaCalls::operator=(const TPZCudaCalls &copy) {
	if(&copy == this){
		return *this;
	}
	handle_cusparse = copy.handle_cusparse;
	cusparse_h = copy.cusparse_h;
	handle_cublas = copy.handle_cublas;
	cublas_h = copy.cublas_h;

	return *this;
}

void TPZCudaCalls::Multiply(bool trans, int *m, int *n, int *k, REAL *A, int *strideA, 
	REAL *B, int *strideB,  REAL *C, int *strideC, REAL alpha, int nmatrices) {
	int numBlocks = (nmatrices + NT - 1) / NT;
	MatrixMultiplicationKernel<<<numBlocks,NT>>> (trans, m, n, k, A, strideA, B, strideB, C, strideC, alpha, nmatrices);
	cudaDeviceSynchronize();
	cudaError_t error = cudaGetLastError();
	if (error != cudaSuccess) {
		std::string error_string = cudaGetErrorString(error);
		std::string error_message = "failed to perform MatrixMultiplicationKernel: " + error_string;
		throw std::runtime_error(error_message);      
	}

}

void TPZCudaCalls::GatherOperation(int n, REAL *x, REAL *y, int *id) {
	if(cusparse_h == false) {
		cusparse_h = true;
		cusparseStatus_t result = cusparseCreate(&handle_cusparse);
		if (result != CUSPARSE_STATUS_SUCCESS) {
			throw std::runtime_error("failed to initialize cuSparse");      
		}			
	}
	cusparseStatus_t result = cusparseDgthr(handle_cusparse, n, x, y, id, CUSPARSE_INDEX_BASE_ZERO);
	if (result != CUSPARSE_STATUS_SUCCESS) {
		throw std::runtime_error("failed to perform cusparseDgthr");      
	}	
}

void TPZCudaCalls::ScatterOperation(int n, REAL *x, REAL *y, int *id) {
	if(cusparse_h == false) {
		cusparse_h = true;
		cusparseStatus_t result = cusparseCreate(&handle_cusparse);
		if (result != CUSPARSE_STATUS_SUCCESS) {
			throw std::runtime_error("failed to initialize cuSparse");      
		}			
	}
	cusparseStatus_t result = cusparseDsctr(handle_cusparse, n, x, id, y, CUSPARSE_INDEX_BASE_ZERO);
	if (result != CUSPARSE_STATUS_SUCCESS) {
		throw std::runtime_error("failed to perform cusparseDsctr");      
	}	
}

void TPZCudaCalls::DaxpyOperation(int n, double alpha, double *x, double *y) {
	if(cublas_h == false) {
		cublas_h = true;
		cublasStatus_t result = cublasCreate(&handle_cublas);
		if (result != CUBLAS_STATUS_SUCCESS) {
			throw std::runtime_error("failed to initialize cuBLAS");      
		}			
	}
	cublasStatus_t result = cublasDaxpy(handle_cublas, n, &alpha, x, 1., y, 1.);
	if (result != CUBLAS_STATUS_SUCCESS) {
		throw std::runtime_error("failed to perform cublasDaxpy");      
	}	
}

void TPZCudaCalls::SpMV(int opt, int sym, int m, int k, int nnz, REAL alpha, REAL *csrVal, int *csrRowPtr, int *csrColInd, REAL *B, REAL *C) {
	if(cusparse_h == false) {
		cusparse_h = true;
		cusparseStatus_t result = cusparseCreate(&handle_cusparse);
		if (result != CUSPARSE_STATUS_SUCCESS) {
			throw std::runtime_error("failed to initialize cuSparse");      
		}			
	}
	cusparseMatDescr_t descr;
	cusparseCreateMatDescr(&descr);
	cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
    if(sym == 0) {
	   cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);        
    } 
    else {
        cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_SYMMETRIC); 
    }
    cusparseOperation_t op;
    if(opt == 0) { 
        op = CUSPARSE_OPERATION_NON_TRANSPOSE;
    } else {
        op = CUSPARSE_OPERATION_TRANSPOSE;
    }

	REAL beta = 0.;
	cusparseStatus_t result = cusparseDcsrmv(handle_cusparse, op, m, k, nnz, &alpha, descr, csrVal, csrRowPtr, csrColInd, B, &beta, C);
	if (result != CUSPARSE_STATUS_SUCCESS) {
		throw std::runtime_error("failed to perform cusparseDcsrmv");      
	}	
}

void TPZCudaCalls::SpMSpM(int opt, int sym, int m, int n, int k, int nnzA, REAL *csrValA, int *csrRowPtrA, int *csrColIndA, 
	int nnzB, REAL *csrValB, int *csrRowPtrB, int *csrColIndB, 
	int nnzC, REAL *csrValC, int *csrRowPtrC) {
	if(cusparse_h == false) {
		cusparse_h = true;
		cusparseStatus_t result = cusparseCreate(&handle_cusparse);
		if (result != CUSPARSE_STATUS_SUCCESS) {
			throw std::runtime_error("failed to initialize cuSparse");      
		}			
	}

    cusparseMatDescr_t descr;
    cusparseCreateMatDescr(&descr);
    cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
    if(sym == 0) {
       cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);        
    } 
    else {
        cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_SYMMETRIC); 
    }

    cusparseOperation_t op;
    if(opt == 0) { 
        op = CUSPARSE_OPERATION_NON_TRANSPOSE;
    } else {
        op = CUSPARSE_OPERATION_TRANSPOSE;
    }

	int *csrColIndC;
	cudaMalloc((void**)&csrColIndC, sizeof(int)*nnzC);

	cusparseStatus_t result = cusparseDcsrgemm(handle_cusparse, op, CUSPARSE_OPERATION_NON_TRANSPOSE, m, n, k, 
		descr, nnzA, csrValA, csrRowPtrA, csrColIndA, 
		descr, nnzB, csrValB, csrRowPtrB, csrColIndB,
		descr, csrValC, csrRowPtrC, csrColIndC);
	if (result != CUSPARSE_STATUS_SUCCESS) {
		throw std::runtime_error("failed to perform cusparseDcsrgemm");      
	}	
}

void TPZCudaCalls::ComputeSigma(bool update_mem, int npts, REAL *glob_delta_strain, REAL *glob_sigma, REAL lambda, REAL mu, REAL mc_phi, REAL mc_psi, REAL mc_cohesion, REAL *dPlasticStrain,  
	REAL *dMType, REAL *dAlpha, REAL *dSigma, REAL *dStrain, REAL *weight) {
	
	int numBlocks = (npts + NT - 1) / NT;
	ComputeSigmaKernel<<<numBlocks,NT>>> (update_mem, npts, glob_delta_strain, glob_sigma, lambda, mu, mc_phi, mc_psi, mc_cohesion, dPlasticStrain, dMType, dAlpha, dSigma, dStrain, weight);
	cudaDeviceSynchronize();
	cudaError_t error = cudaGetLastError();
	if (error != cudaSuccess) {
		std::string error_string = cudaGetErrorString(error);
		std::string error_message = "failed to perform ComputeSigmaKernel: " + error_string;
		throw std::runtime_error(error_message);      
	}
}

void TPZCudaCalls::MatrixAssembleGS(REAL *Kc, int first_el, int last_el, int64_t *el_color_index, REAL *weight, int *dof_indexes,
	REAL *storage, int *rowsizes, int *colsizes, int *rowfirstindex, int *colfirstindex, int *matrixposition) {
	int nel = last_el - first_el;
	int numBlocks = (nel + NT_sm - 1) / NT_sm;

	MatrixAssembleKernelGS<<<numBlocks,NT_sm>>>(nel, Kc, el_color_index, weight, dof_indexes, 
	storage, rowsizes, colsizes, rowfirstindex, colfirstindex, matrixposition);
	cudaDeviceSynchronize();
	cudaError_t error = cudaGetLastError();
	if (error != cudaSuccess) {
		std::string error_string = cudaGetErrorString(error);
		std::string error_message = "failed to perform MatrixAssembleKernel: " + error_string;
		throw std::runtime_error(error_message);      
	}
}

void TPZCudaCalls::MatrixAssemble(REAL *Kg, int first_el, int last_el, int64_t *el_color_index,  REAL *weight, int *dof_indexes,
	REAL *storage, int *rowsizes, int *colsizes, int *rowfirstindex, int *colfirstindex, int *matrixposition, int *ia_to_sequence, int *ja_to_sequence) {
	int nel = last_el - first_el;
	int numBlocks = (nel + NT_sm - 1) / NT_sm;
	MatrixAssembleKernel<<<numBlocks,NT_sm>>> (nel, Kg, el_color_index, weight, dof_indexes, storage, rowsizes, colsizes, rowfirstindex, colfirstindex,
	matrixposition, ia_to_sequence, ja_to_sequence); 
	cudaDeviceSynchronize();
	cudaError_t error = cudaGetLastError();
	if (error != cudaSuccess) {
		std::string error_string = cudaGetErrorString(error);
		std::string error_message = "failed to perform MatrixAssembleKernel: " + error_string;
		throw std::runtime_error(error_message);      
	}
}

void TPZCudaCalls::DeToDevice(REAL lambda, REAL mu) {
		REAL De_host[] = {lambda + 2.0*mu, 0, lambda, 0, mu, 0, lambda, 0, lambda + 2.0*mu};
		cudaMemcpyToSymbol(De, &De_host, 9 * sizeof(REAL));
	}


void TPZCudaCalls::SolveCG(int n, int nnzA, REAL *csrValA, int *csrRowPtrA, int *csrColIndA, REAL *r, REAL *x) {
    if(cusparse_h == false) {
        cusparse_h = true;
        cusparseStatus_t result = cusparseCreate(&handle_cusparse);
        if (result != CUSPARSE_STATUS_SUCCESS) {
            throw std::runtime_error("failed to initialize cuSparse");      
        }           
    }

    if(cublas_h == false) {
        cublas_h = true;
        cublasStatus_t result = cublasCreate(&handle_cublas);
        if (result != CUBLAS_STATUS_SUCCESS) {
            throw std::runtime_error("failed to initialize cuBLAS");      
        }           
    }

    cusparseMatDescr_t descr;
    cusparseCreateMatDescr(&descr);
    cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_SYMMETRIC);
    cusparseSetMatFillMode(descr, CUSPARSE_FILL_MODE_UPPER);
    cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);

    REAL alpha = 1.0;
    REAL alpham1 = -1.0;
    REAL beta = 0.0;
    REAL r0 = 0.;
    REAL b;
    REAL r1;
    REAL dot;
    REAL a;
    REAL na;

    REAL *d_Ax;
    REAL *d_p;
    cudaMalloc((void **)&d_Ax, n*sizeof(REAL));
    cudaMalloc((void **)&d_p, n*sizeof(REAL));

    cusparseDcsrmv(handle_cusparse,CUSPARSE_OPERATION_NON_TRANSPOSE, n, n, nnzA, &alpha, descr, csrValA, csrRowPtrA, csrColIndA, x, &beta, d_Ax);
    cublasDaxpy(handle_cublas, n, &alpham1, d_Ax, 1, r, 1);


    cublasDdot(handle_cublas, n, r, 1, r, 1, &r1);

    const REAL tol = 1.e-5;
    const int max_iter = 10000;
    int k;

    k = 1;

    while (r1 > tol*tol && k <= max_iter)
    {
        if (k > 1)
        {
            b = r1 / r0;
            cublasDscal(handle_cublas, n, &b, d_p, 1);
            cublasDaxpy(handle_cublas, n, &alpha, r, 1, d_p, 1);
        }
        else
        {
            cublasDcopy(handle_cublas, n, r, 1, d_p, 1);
        }
        cusparseDcsrmv(handle_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE, n, n, nnzA, &alpha, descr, csrValA, csrRowPtrA, csrColIndA, d_p, &beta, d_Ax);
        cublasDdot(handle_cublas, n, d_p, 1, d_Ax, 1, &dot);
        a = r1 / dot;

        cublasDaxpy(handle_cublas, n, &a, d_p, 1, x, 1);
        na = -a;
        cublasDaxpy(handle_cublas, n, &na, d_Ax, 1, r, 1);

        r0 = r1;
        cublasDdot(handle_cublas, n, r, 1, r, 1, &r1);
        cudaThreadSynchronize();
        k++;
    }
    cudaFree(d_p);
    cudaFree(d_Ax);
}