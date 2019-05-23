#include "TPZIntPointsFEM.h"
#include "TPZTensor.h"
#include "pzmatrix.h"
#include <stdlib.h>
#include "TPZTensor.h"
#include "TPZVTKGeoMesh.h"
#include "pzintel.h"
#include "pzskylstrmatrix.h"
#include "Timer.h"

#include "SpectralDecompKernels.h"
#include "MatMulKernels.h"
#include "StressStrainKernels.h"
#include "SigmaProjectionKernels.h"

REAL timeGatherSolution = 0;
REAL timeDeltaStrain = 0;
REAL timeElasticStrain = 0;
REAL timePlasticStrain = 0;
REAL timeComputeStress = 0;
REAL timeComputeStrain = 0;
REAL timeSpectralDecomposition = 0;
REAL timeProjectSigma = 0;
REAL timeStressCompleteTensor = 0;
REAL timeNodalForces = 0;
REAL timeColoredAssemble = 0;

#define NT 128
#define NT_MULT 128

TPZIntPointsFEM::TPZIntPointsFEM() :
		fDim(-1), fBoundaryElements(), fCmesh(0), fNpts(-1), fNphis(-1), fElemColor(
				0), fMaterial(0), fRhs(0, 0), fRhsBoundary(0, 0), fSolution(0,
				0), fPlasticStrain(0, 0), fStorage(0), fRowSizes(0), fColSizes(
				0), fMatrixPosition(0), fRowFirstIndex(0), fColFirstIndex(0), fIndexes(
				0), fIndexesColor(0), fWeight(0), fRowPtr(0), fColInd(0),
				fTimer() {

	dRhs = new REAL[0];
	dRhsBoundary = new REAL[0];
	dSolution = new REAL[0];
	dPlasticStrain = new REAL[0];
	dStorage = new REAL[0];
	dRowSizes = new int[0];
	dColSizes = new int[0];
	dMatrixPosition = new int[0];
	dRowFirstIndex = new int[0];
	dColFirstIndex = new int[0];
	dIndexes = new int[0];
	dIndexesColor = new int[0];
	dWeight = new REAL[0];

	dRowPtr = new int[0];
	dColInd = new int[0];


}

TPZIntPointsFEM::TPZIntPointsFEM(TPZCompMesh *cmesh, int materialid) :
		fDim(-1), fBoundaryElements(), fCmesh(0), fNpts(-1), fNphis(-1), fElemColor(
				0), fMaterial(0), fRhs(0, 0), fRhsBoundary(0, 0), fSolution(0,
				0), fPlasticStrain(0, 0), fStorage(0), fRowSizes(0), fColSizes(
				0), fMatrixPosition(0), fRowFirstIndex(0), fColFirstIndex(0), fIndexes(
				0), fIndexesColor(0), fWeight(0), fRowPtr(0), fColInd(0),
				fTimer () {
	SetCompMesh(cmesh);
	SetMaterialId(materialid);

	dRhs = new REAL[0];
	dRhsBoundary = new REAL[0];
	dSolution = new REAL[0];
	dPlasticStrain = new REAL[0];
	dStorage = new REAL[0];
	dRowSizes = new int[0];
	dColSizes = new int[0];
	dMatrixPosition = new int[0];
	dRowFirstIndex = new int[0];
	dColFirstIndex = new int[0];
	dIndexes = new int[0];
	dIndexesColor = new int[0];
	dWeight = new REAL[0];

	dRowPtr = new int[0];
	dColInd = new int[0];
}

TPZIntPointsFEM::~TPZIntPointsFEM() {
	cudaFree(dRhs);
	cudaFree(dRhsBoundary);
	cudaFree(dSolution);
	cudaFree(dPlasticStrain);
	cudaFree(dStorage);
	cudaFree(dRowSizes);
	cudaFree(dColSizes);
	cudaFree(dMatrixPosition);
	cudaFree(dRowFirstIndex);
	cudaFree(dColFirstIndex);
	cudaFree(dIndexes);
	cudaFree(dIndexesColor);
	cudaFree(dWeight);

	cudaFree(dRowPtr);
	cudaFree(dColInd);

	cublasDestroy(handle_cublas);
	cusparseDestroy(handle_cusparse);
}

TPZIntPointsFEM::TPZIntPointsFEM(const TPZIntPointsFEM &copy) {
	fDim = copy.fDim;
	fBoundaryElements = copy.fBoundaryElements;
	fCmesh = copy.fCmesh;
	fNpts = copy.fNpts;
	fNphis = copy.fNphis;
	fElemColor = copy.fElemColor;
	fMaterial = copy.fMaterial;

	fTimer = copy.fTimer;

	fRhs = copy.fRhs;
	fRhsBoundary = copy.fRhsBoundary;
	fSolution = copy.fSolution;
	fPlasticStrain = copy.fPlasticStrain;
	fStorage = copy.fStorage;
	fColSizes = copy.fColSizes;
	fRowSizes = copy.fRowSizes;
	fMatrixPosition = copy.fMatrixPosition;
	fRowFirstIndex = copy.fRowFirstIndex;
	fColFirstIndex = copy.fColFirstIndex;
	fIndexes = copy.fIndexes;
	fIndexesColor = copy.fIndexesColor;
	fWeight = copy.fWeight;

	fRowPtr = copy.fRowPtr;
	fColInd = copy.fColInd;

	handle_cusparse = copy.handle_cusparse;
	handle_cublas = copy.handle_cublas;

	dRhs = copy.dRhs;
	dRhsBoundary = copy.dRhsBoundary;
	dSolution = copy.dSolution;
	dPlasticStrain = copy.dPlasticStrain;
	dStorage = copy.dStorage;
	dRowSizes = copy.dRowSizes;
	dColSizes = copy.dColSizes;
	dMatrixPosition = copy.dMatrixPosition;
	dRowFirstIndex = copy.dRowFirstIndex;
	dColFirstIndex = copy.dColFirstIndex;
	dIndexes = copy.dIndexes;
	dIndexesColor = copy.dIndexesColor;
	dWeight = copy.dWeight;

	dRowPtr = copy.dRowPtr;
	dColInd = copy.dColInd;
}

TPZIntPointsFEM &TPZIntPointsFEM::operator=(const TPZIntPointsFEM &copy) {
	if (&copy == this) {
		return *this;
	}

	fDim = copy.fDim;
	fBoundaryElements = copy.fBoundaryElements;
	fCmesh = copy.fCmesh;
	fNpts = copy.fNpts;
	fNphis = copy.fNphis;
	fElemColor = copy.fElemColor;
	fMaterial = copy.fMaterial;

	fTimer = copy.fTimer;

	fRhs = copy.fRhs;
	fRhsBoundary = copy.fRhsBoundary;
	fSolution = copy.fSolution;
	fPlasticStrain = copy.fPlasticStrain;
	fStorage = copy.fStorage;
	fColSizes = copy.fColSizes;
	fRowSizes = copy.fRowSizes;
	fMatrixPosition = copy.fMatrixPosition;
	fRowFirstIndex = copy.fRowFirstIndex;
	fColFirstIndex = copy.fColFirstIndex;
	fIndexes = copy.fIndexes;
	fIndexesColor = copy.fIndexesColor;
	fWeight = copy.fWeight;

	fRowPtr = copy.fRowPtr;
	fColInd = copy.fColInd;

	handle_cusparse = copy.handle_cusparse;
	handle_cublas = copy.handle_cublas;

	dRhs = copy.dRhs;
	dRhsBoundary = copy.dRhsBoundary;
	dSolution = copy.dSolution;
	dPlasticStrain = copy.dPlasticStrain;
	dStorage = copy.dStorage;
	dRowSizes = copy.dRowSizes;
	dColSizes = copy.dColSizes;
	dMatrixPosition = copy.dMatrixPosition;
	dRowFirstIndex = copy.dRowFirstIndex;
	dColFirstIndex = copy.dColFirstIndex;
	dIndexes = copy.dIndexes;
	dIndexesColor = copy.dIndexesColor;
	dWeight = copy.dWeight;

	dRowPtr = copy.dRowPtr;
	dColInd = copy.dColInd;

	return *this;
}

void TPZIntPointsFEM::SetTimerConfig(Timer::WhichUnit unit) {
	fTimer.TimerConfig(unit);
}

void TPZIntPointsFEM::TransferDataStructure() {

	std::cout << "Initializing libraries contexts ... " << std::endl;
	cublasCreate(&handle_cublas);
	cusparseCreate(&handle_cusparse);

	int64_t neq = fCmesh->NEquations();
	int64_t nelem = fColSizes.size();
    int64_t nnz = fStorage.size();


	std::cout << "Allocating and transfering data to GPU ... " << std::endl;
	fTimer.Start();
	cudaMalloc((void**) &dRhs, neq * sizeof(REAL));

	cudaMalloc((void**) &dRhsBoundary, neq * sizeof(REAL));
	cudaMemcpy(dRhsBoundary, &fRhsBoundary(0, 0), neq * sizeof(REAL), cudaMemcpyHostToDevice);

	cudaMalloc((void**) &dSolution, neq * sizeof(REAL));

	cudaMalloc((void**) &dPlasticStrain, fDim * fNpts * sizeof(REAL));

	cudaMalloc((void**) &dStorage, fStorage.size() * sizeof(REAL));
	cudaMemcpy(dStorage, &fStorage[0], fStorage.size() * sizeof(REAL), cudaMemcpyHostToDevice);

	cudaMalloc((void**) &dIndexes, fIndexes.size() * sizeof(int));
	cudaMemcpy(dIndexes, &fIndexes[0], fIndexes.size() * sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc((void**) &dIndexesColor, fIndexesColor.size() * sizeof(int));
	cudaMemcpy(dIndexesColor, &fIndexesColor[0],
			fIndexesColor.size() * sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc((void**) &dWeight, fWeight.size() * sizeof(REAL));
	cudaMemcpy(dWeight, &fWeight[0], fWeight.size() * sizeof(REAL), cudaMemcpyHostToDevice);

#ifdef USING_CSRMV_MULT
	cudaMalloc((void**) &dRowPtr, (fNpts + 1) * sizeof(int));
	cudaMemcpy(dRowPtr, &fRowPtr[0], (fNpts + 1) * sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc((void**) &dColInd, nnz * sizeof(int));
	cudaMemcpy(dColInd, &fColInd[0], nnz * sizeof(int), cudaMemcpyHostToDevice);
#else
	cudaMalloc((void**) &dRowSizes, nelem * sizeof(int));
	cudaMemcpy(dRowSizes, &fRowSizes[0], nelem * sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc((void**) &dColSizes, nelem * sizeof(int));
	cudaMemcpy(dColSizes, &fColSizes[0], nelem * sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc((void**) &dMatrixPosition, nelem * sizeof(int));
	cudaMemcpy(dMatrixPosition, &fMatrixPosition[0], nelem * sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc((void**) &dRowFirstIndex, nelem * sizeof(int));
	cudaMemcpy(dRowFirstIndex, &fRowFirstIndex[0], nelem * sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc((void**) &dColFirstIndex, nelem * sizeof(int));
	cudaMemcpy(dColFirstIndex, &fColFirstIndex[0], nelem * sizeof(int), cudaMemcpyHostToDevice);

#endif
}

void TPZIntPointsFEM::GatherSolution(REAL *gather_solution) {
	fTimer.Start();
	cusparseDgthr(handle_cusparse, fDim * fNphis, dSolution, gather_solution, dIndexes, CUSPARSE_INDEX_BASE_ZERO);
    fTimer.Stop();
    timeGatherSolution+= fTimer.ElapsedTime();
}

void TPZIntPointsFEM::DeltaStrain(REAL *gather_solution, REAL *delta_strain) {
	int64_t nelem = fRowSizes.size();
	int numBlocks = (nelem + NT - 1) / NT;

#ifdef USING_CUBLAS_MULT //Using cuBLAS matrix-multiplication (each multiplication is done in one thread through cuBLAS library)
	cublasOperation_t trans = CUBLAS_OP_N;

	fTimer.Start();
	MatMulcuBLASKernel<<<numBlocks, NT>>>(trans, nelem, dStorage, dRowSizes, dColSizes, dMatrixPosition, dRowFirstIndex, dColFirstIndex, fNpts, fNphis, gather_solution, delta_strain);
	cudaDeviceSynchronize();
	fTimer.Stop();
	timeDeltaStrain+= fTimer.ElapsedTime();

#elif USING_CSRMV_MULT //Using cuSPARSE Spmv
	int nnz = fStorage.size();
	REAL alpha = 1.;
	REAL beta = 0.;
	cusparseMatDescr_t descr;
	cusparseCreateMatDescr(&descr);
	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);
	int m = fNpts;
	int n = 1;
	int k = fNphis;

	fTimer.Start();
	cusparseDcsrmm(handle_cusparse,CUSPARSE_OPERATION_NON_TRANSPOSE, m, n, k, nnz, &alpha, descr, dStorage, dRowPtr, dColInd, &gather_solution[0], k, &beta, &delta_strain[0], m);
	cusparseDcsrmm(handle_cusparse,CUSPARSE_OPERATION_NON_TRANSPOSE, m, n, k, nnz, &alpha, descr, dStorage, dRowPtr, dColInd, &gather_solution[fNphis], k, &beta, &delta_strain[fNpts], m);
//	cusparseDcsrmv(handle_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE, m, k, nnz, &alpha, descr, dStorage, dRowPtr, dColInd, &gather_solution[0], &beta, &delta_strain[0]);
//	cusparseDcsrmv(handle_cusparse, CUSPARSE_OPERATION_NON_TRANSPOSE, m, k, nnz, &alpha, descr, dStorage, dRowPtr, dColInd, &gather_solution[fNphis], &beta, &delta_strain[fNpts]);
	fTimer.Stop();
	timeDeltaStrain+= fTimer.ElapsedTime();
#elif USING_CUBLASBATCHED_MULT //Using cuBlas matrix-multiplication using gemmBatched
    int64_t cols = fColSizes[0];
    int64_t rows = fRowSizes[0];
	REAL alpha = 1.;
	REAL beta = 0.;
	cublasOperation_t transA = CUBLAS_OP_N;
	cublasOperation_t transB = CUBLAS_OP_N;

	fTimer.Start();
    cublasDgemmStridedBatched(handle_cublas, transA, transB, rows, 1, cols, &alpha, dStorage, rows, rows*cols, &gather_solution[0], cols, cols*1, &beta, &delta_strain[0], rows, rows*1, nelem);
    cublasDgemmStridedBatched(handle_cublas, transA, transB, rows, 1, cols, &alpha, dStorage, rows, rows*cols, &gather_solution[fNphis], cols, cols*1, &beta,  &delta_strain[fNpts], rows, rows*1, nelem);
	fTimer.Stop();
	timeDeltaStrain+= fTimer.ElapsedTime();

#else //Using a loop over each line of the matrices
	bool trans = false;

	fTimer.Start();
	MatMulKernel<<<numBlocks, NT>>>(trans, nelem, dStorage, dRowSizes, dColSizes, dMatrixPosition, dRowFirstIndex, dColFirstIndex, fNpts, fNphis, gather_solution, delta_strain);
	cudaDeviceSynchronize();
	fTimer.Stop();
	timeDeltaStrain+= fTimer.ElapsedTime();
#endif

}

void TPZIntPointsFEM::ElasticStrain(REAL *delta_strain, REAL *plastic_strain, REAL *elastic_strain) {
	cudaMemcpy(elastic_strain, &delta_strain[0], fDim * fNpts * sizeof(REAL), cudaMemcpyDeviceToDevice);
	cudaMemset(plastic_strain, 0, fDim * fNpts * sizeof(REAL));

    fTimer.Start();
	REAL a = -1.;
	cublasDaxpy(handle_cublas, fDim * fNpts, &a, &plastic_strain[0], 1, &elastic_strain[0], 1);
	fTimer.Stop();
	timeElasticStrain+= fTimer.ElapsedTime();
}

void TPZIntPointsFEM::PlasticStrain(REAL *delta_strain, REAL *elastic_strain, REAL *plastic_strain) {
	cudaMemcpy(plastic_strain, &delta_strain[0], fDim * fNpts * sizeof(REAL), cudaMemcpyDeviceToDevice);

    fTimer.Start();
	REAL a = -1.;
	cublasDaxpy(handle_cublas, fDim * fNpts, &a, &elastic_strain[0], 1, &plastic_strain[0], 1);
	fTimer.Stop();
	timePlasticStrain+= fTimer.ElapsedTime();
}

void TPZIntPointsFEM::ComputeStress(REAL *elastic_strain, REAL *sigma) {
	REAL lambda = fMaterial->GetPlasticModel().fER.Lambda();
	REAL mu = fMaterial->GetPlasticModel().fER.Mu();
	int numBlocks = (fNpts / fDim + NT - 1) / NT;

	fTimer.Start();
	ComputeStressKernel<<<numBlocks, NT>>>(fNpts, fDim, elastic_strain, sigma, mu, lambda);
	cudaDeviceSynchronize();
	fTimer.Stop();
	timeComputeStress+= fTimer.ElapsedTime();
}

void TPZIntPointsFEM::ComputeStrain(REAL *sigma, REAL *elastic_strain) {
	REAL E = fMaterial->GetPlasticModel().fER.E();
	REAL nu = fMaterial->GetPlasticModel().fER.Poisson();
	int numBlocks = (fNpts / fDim + NT - 1) / NT;

	fTimer.Start();
	ComputeStrainKernel<<<numBlocks, NT>>>(fNpts, fDim, sigma, elastic_strain, nu, E, dWeight);
	cudaDeviceSynchronize();
	fTimer.Stop();
	timeComputeStrain+= fTimer.ElapsedTime();
}

void TPZIntPointsFEM::SpectralDecomposition(REAL *sigma_trial, REAL *eigenvalues, REAL *eigenvectors) {
	int numBlocks = (fNpts / fDim + NT - 1) / NT;

	fTimer.Start();
	SpectralDecompositionKernel<<<numBlocks, NT>>>(fNpts, fDim, sigma_trial, eigenvalues, eigenvectors);
	cudaDeviceSynchronize();
	fTimer.Stop();
	timeSpectralDecomposition+= fTimer.ElapsedTime();
}

void TPZIntPointsFEM::ProjectSigma(REAL *eigenvalues, REAL *sigma_projected) {

	REAL mc_phi = fMaterial->GetPlasticModel().fYC.Phi();
	REAL mc_psi = fMaterial->GetPlasticModel().fYC.Psi();
	REAL mc_cohesion = fMaterial->GetPlasticModel().fYC.Cohesion();
	REAL K = fMaterial->GetPlasticModel().fER.K();
	REAL G = fMaterial->GetPlasticModel().fER.G();

	REAL *m_type;
	cudaMalloc((void**) &m_type, fNpts / fDim * sizeof(REAL));

	REAL *alpha;
	cudaMalloc((void**) &alpha, fNpts / fDim * sizeof(REAL));

	int numBlocks = (fNpts / fDim + NT - 1) / NT;

	fTimer.Start();
	ProjectSigmaKernel<<<numBlocks, NT>>>(fNpts, fDim, mc_phi, mc_psi, mc_cohesion, K, G, eigenvalues, sigma_projected, m_type, alpha);
	cudaDeviceSynchronize();
	fTimer.Stop();
	timeProjectSigma+= fTimer.ElapsedTime();
}

void TPZIntPointsFEM::StressCompleteTensor(REAL *sigma_projected, REAL *eigenvectors, REAL *sigma) {
	int numBlocks = (fNpts / fDim + NT - 1) / NT;

	fTimer.Start();
	StressCompleteTensorKernel<<<numBlocks, NT>>>(fNpts, fDim, sigma_projected, eigenvectors, sigma, dWeight);
	cudaDeviceSynchronize();
	fTimer.Stop();
	timeStressCompleteTensor+= fTimer.ElapsedTime();
}

void TPZIntPointsFEM::NodalForces(REAL *sigma, REAL *nodal_forces) {
	int64_t nelem = fRowSizes.size();
	int numBlocks = (nelem + NT - 1) / NT;


#ifdef USING_CUBLAS_MULT //Using cuBLAS matrix-multiplication (each multiplication is done in one thread through cuBLAS library)
	cublasOperation_t transA = CUBLAS_OP_T;

    fTimer.Start();
	MatMulcuBLASKernel<<<numBlocks, NT>>>(transA, nelem, dStorage, dRowSizes, dColSizes, dMatrixPosition, dRowFirstIndex, dColFirstIndex, fNpts, fNphis, sigma, nodal_forces);
	cudaDeviceSynchronize();
	fTimer.Stop();
	timeNodalForces+= fTimer.ElapsedTime();

#elif USING_CSRMV_MULT //Using cuSPARSE Spmv
	int nnz = fStorage.size();
	REAL alpha = -1.;
	REAL beta = 0.;
	int m = fNpts;
	int n = 1;
	int k = fNphis;

	cusparseMatDescr_t descr;
	cusparseCreateMatDescr(&descr);
	cusparseSetMatType(descr, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descr, CUSPARSE_INDEX_BASE_ZERO);


	fTimer.Start();
	cusparseDcsrmm(handle_cusparse,CUSPARSE_OPERATION_TRANSPOSE, m, n, k, nnz, &alpha, descr, dStorage, dRowPtr, dColInd, &sigma[0], m, &beta, &nodal_forces[0], k);
	cusparseDcsrmm(handle_cusparse,CUSPARSE_OPERATION_TRANSPOSE, m, n, k, nnz, &alpha, descr, dStorage, dRowPtr, dColInd, &sigma[fNpts], m, &beta, &nodal_forces[fNphis], k);
//	cusparseDcsrmv(handle_cusparse, CUSPARSE_OPERATION_TRANSPOSE, m, k, nnz, &alpha, descr, dStorage, dRowPtr, dColInd, &sigma[0], &beta, &nodal_forces[0]);
//	cusparseDcsrmv(handle_cusparse, CUSPARSE_OPERATION_TRANSPOSE, m, k, nnz, &alpha, descr, dStorage, dRowPtr, dColInd, &sigma[fNpts], &beta, &nodal_forces[fNphis]);
	fTimer.Stop();
	timeNodalForces+= fTimer.ElapsedTime();

#elif USING_CUBLASBATCHED_MULT //Using cuBlas matrix-multiplication using gemmBatched
    int64_t cols = fColSizes[0];
    int64_t rows = fRowSizes[0];
	REAL alpha = -1.;
	REAL beta = 0.;

	cublasOperation_t transA = CUBLAS_OP_T;
	cublasOperation_t transB = CUBLAS_OP_N;

	fTimer.Start();
	cublasDgemmStridedBatched(handle_cublas, CUBLAS_OP_T, CUBLAS_OP_N, cols, 1, rows, &alpha, dStorage, rows, rows*cols, &sigma[0], rows, rows*1, &beta, &nodal_forces[0], cols, cols*1, nelem);
    cublasDgemmStridedBatched(handle_cublas, CUBLAS_OP_T, CUBLAS_OP_N, cols, 1, rows, &alpha, dStorage, rows, rows*cols, &sigma[fNpts], rows, rows*1, &beta,  &nodal_forces[fNphis], cols, cols*1, nelem);
	fTimer.Stop();
	timeNodalForces+= fTimer.ElapsedTime();
#else //Using a loop over each line of the matrices
	bool trans = true;

	fTimer.Start();
	MatMulKernel<<<numBlocks, NT>>>(trans, nelem, dStorage, dRowSizes, dColSizes, dMatrixPosition, dRowFirstIndex, dColFirstIndex, fNpts, fNphis, sigma, nodal_forces);
	cudaDeviceSynchronize();
	fTimer.Stop();
	timeNodalForces+= fTimer.ElapsedTime();
#endif

}

void TPZIntPointsFEM::ColoredAssemble(REAL *nodal_forces, REAL *residual) {
	int64_t ncolor = *std::max_element(fElemColor.begin(), fElemColor.end())
			+ 1;
	int64_t sz = fIndexes.size();
	int64_t neq = fCmesh->NEquations();

    fTimer.Start();
	cusparseDsctr(handle_cusparse, sz, &nodal_forces[0], &dIndexesColor[0], &residual[0], CUSPARSE_INDEX_BASE_ZERO);

	int64_t colorassemb = ncolor / 2.;
	REAL alpha = 1.;
	while (colorassemb > 0) {

		int64_t firsteq = (ncolor - colorassemb) * neq;
		cublasDaxpy(handle_cublas, colorassemb * neq, &alpha, &residual[firsteq], 1., &residual[0], 1.);

		ncolor -= colorassemb;
		colorassemb = ncolor / 2;
	}
	fTimer.Stop();
	timeColoredAssemble+= fTimer.ElapsedTime();
}

void TPZIntPointsFEM::AssembleResidual() {
	int64_t ncolor = *std::max_element(fElemColor.begin(), fElemColor.end()) + 1;
	int64_t neq = fCmesh->NEquations();

	REAL *gather_solution;
	cudaMalloc((void**) &gather_solution, fNpts * sizeof(REAL));
	cudaMemset(gather_solution, 0, fNpts * sizeof(REAL));

	REAL *delta_strain;
	cudaMalloc((void**) &delta_strain, fDim * fNpts * sizeof(REAL));
	cudaMemset(delta_strain, 0, fDim * fNpts * sizeof(REAL));

	REAL *elastic_strain;
	cudaMalloc((void**) &elastic_strain, fDim * fNpts * sizeof(REAL));

	REAL *sigma_trial;
	cudaMalloc((void**) &sigma_trial, fDim * fNpts * sizeof(REAL));

	REAL *eigenvalues;
	cudaMalloc((void**) &eigenvalues, 3 * fNpts / fDim * sizeof(REAL));

	REAL *eigenvectors;
	cudaMalloc((void**) &eigenvectors, 9 * fNpts / fDim * sizeof(REAL));

	REAL *sigma_projected;
	cudaMalloc((void**) &sigma_projected, 3 * fNpts / fDim * sizeof(REAL));

	REAL *sigma;
	cudaMalloc((void**) &sigma, fDim * fNpts * sizeof(REAL));

	REAL *nodal_forces;
	cudaMalloc((void**) &nodal_forces, fDim * fNpts * sizeof(REAL));
	cudaMemset(nodal_forces, 0, fDim * fNpts * sizeof(REAL));

	REAL *residual;
	cudaMalloc((void**) &residual, neq * ncolor * sizeof(REAL));
	cudaMemset(residual, 0, neq * ncolor * sizeof(REAL));

	cudaMemcpy(dSolution, &fSolution(0, 0), neq * sizeof(REAL), cudaMemcpyHostToDevice);
	GatherSolution(gather_solution);
	DeltaStrain(gather_solution, delta_strain);
	ElasticStrain(delta_strain, dPlasticStrain, elastic_strain);
	ComputeStress(elastic_strain, sigma_trial);
	SpectralDecomposition(sigma_trial, eigenvalues, eigenvectors);
	ProjectSigma(eigenvalues, sigma_projected);
	StressCompleteTensor(sigma_projected, eigenvectors, sigma);
	NodalForces(sigma, nodal_forces);
	ColoredAssemble(nodal_forces, residual);

//update strain
	ComputeStrain(sigma, elastic_strain);
	PlasticStrain(delta_strain, elastic_strain, dPlasticStrain);

    //add boundary contribution
	REAL a = 1.;
	cublasDaxpy(handle_cublas, neq, &a, &dRhsBoundary[0], 1, &residual[0], 1);

	fRhs.Resize(neq, 1);
	cudaMemcpy(&fRhs(0,0), residual, neq * sizeof(REAL), cudaMemcpyDeviceToHost);

    //cumulated time elapsed to each method
    ofstream file("timing-gpu.txt");
    file << "GatherSolution	    " << timeGatherSolution << fTimer.Unit() << std::endl;
    file << "DeltaStrain	    " << timeDeltaStrain << fTimer.Unit() << std::endl;
    file << "ElasticStrain	    " << timeElasticStrain <<  fTimer.Unit() << std::endl;
    file << "ComputeStress	    " << timeComputeStress <<  fTimer.Unit() << std::endl;
    file << "SpectralDecomp	    " << timeSpectralDecomposition <<  fTimer.Unit() << std::endl;
    file << "ProjectSigma	    " << timeProjectSigma <<  fTimer.Unit() << std::endl;
    file << "StressTensor	    " << timeStressCompleteTensor <<  fTimer.Unit() << std::endl;
    file << "NodalForces	    " << timeNodalForces <<  fTimer.Unit() << std::endl;
    file << "ColoredAssemble	" << timeColoredAssemble << fTimer.Unit() << std::endl;

    file << "ComputeStrain	    " << timeComputeStrain <<  fTimer.Unit() << std::endl;
    file << "PlasticStrain	    " << timeElasticStrain <<  fTimer.Unit() << std::endl;

	cudaFree(gather_solution);
	cudaFree(delta_strain);
	cudaFree(elastic_strain);
	cudaFree(sigma_trial);
	cudaFree(eigenvalues);
	cudaFree(eigenvectors);
	cudaFree(sigma_projected);
	cudaFree(sigma);
	cudaFree(nodal_forces);
	cudaFree(residual);
}

void TPZIntPointsFEM::SetDataStructure() {

	int dim_mesh = (fCmesh->Reference())->Dimension(); // Mesh dimension
	this->SetMeshDimension(dim_mesh);
	int64_t nelem_c = fCmesh->NElements(); // Number of computational elements
	std::vector<int64_t> cel_indexes;

// Number of domain geometric elements
	for (int64_t i = 0; i < nelem_c; i++) {
		TPZCompEl *cel = fCmesh->Element(i);
		if (!cel)
			continue;
		TPZGeoEl *gel = fCmesh->Element(i)->Reference();
		if (!gel)
			continue;
		if (gel->Dimension() == dim_mesh)
			cel_indexes.push_back(cel->Index());
		if (gel->Dimension() < dim_mesh)
			fBoundaryElements.Push(cel->Index());
	}

	if (cel_indexes.size() == 0) {
		DebugStop();
	}

// RowSizes and ColSizes vectors
	int64_t nelem = cel_indexes.size();
	TPZVec < MKL_INT > rowsizes(nelem);
	TPZVec < MKL_INT > colsizes(nelem);

	int64_t npts_tot = 0;
	int64_t nf_tot = 0;
	int it = 0;
	for (auto iel : cel_indexes) {
		//Verification
		TPZCompEl *cel = fCmesh->Element(iel);

		//Integration rule
		TPZInterpolatedElement *cel_inter =
				dynamic_cast<TPZInterpolatedElement *>(cel);
		if (!cel_inter)
			DebugStop();
		TPZIntPoints *int_rule = &(cel_inter->GetIntegrationRule());

		int64_t npts = int_rule->NPoints(); // number of integration points of the element
		int64_t dim = cel_inter->Dimension(); //dimension of the element
		int64_t nf = cel_inter->NShapeF(); // number of shape functions of the element

		rowsizes[it] = dim * npts;
		colsizes[it] = nf;

		it++;

		npts_tot += npts;
		nf_tot += nf;
	}
	this->SetNumberofIntPoints(dim_mesh * npts_tot);
	this->SetNumberofPhis(nf_tot);
	this->SetRowandColSizes(rowsizes, colsizes);

// Dphi matrix, weight and indexes vectors
	TPZFMatrix < REAL > elmatrix;
	TPZVec < REAL > weight(npts_tot);
	TPZManVector < MKL_INT > indexes(dim_mesh * nf_tot);

	int64_t cont1 = 0;
	int64_t cont2 = 0;
	it = 0;
	int64_t contw = 0;
	for (auto iel : cel_indexes) {
		//Verification
		TPZCompEl *cel = fCmesh->Element(iel);

		//Integration rule
		TPZInterpolatedElement *cel_inter =
				dynamic_cast<TPZInterpolatedElement *>(cel);
		if (!cel_inter)
			DebugStop();
		TPZIntPoints *int_rule = &(cel_inter->GetIntegrationRule());

		int64_t npts = int_rule->NPoints(); // number of integration points of the element
		int64_t dim = cel_inter->Dimension(); //dimension of the element
		int64_t nf = cel_inter->NShapeF(); // number of shape functions of the element

		TPZMaterialData data;
		cel_inter->InitMaterialData(data);

		elmatrix.Resize(dim * npts, nf);
		for (int64_t inpts = 0; inpts < npts; inpts++) {
			TPZManVector < REAL > qsi(dim, 1);
			REAL w;
			int_rule->Point(inpts, qsi, w);
			cel_inter->ComputeRequiredData(data, qsi);
//			weight.Push(w * std::abs(data.detjac)); //weight = w * detjac
			weight[contw] = w * std::abs(data.detjac);
			contw++;

			TPZFMatrix < REAL > axes = data.axes;
			TPZFMatrix < REAL > dphix = data.dphix;
			TPZFMatrix < REAL > dphiXY;
			axes.Transpose();
			axes.Multiply(dphix, dphiXY);

			for (int inf = 0; inf < nf; inf++) {
				for (int idim = 0; idim < dim; idim++)
					elmatrix(inpts * dim + idim, inf) = dphiXY(idim, inf);
			}
		}

#ifdef USING_CSRMV_MULT
		elmatrix.Transpose();
#endif
		this->SetElementMatrix(it, elmatrix);
		it++;

		//Indexes vector
		int64_t ncon = cel->NConnects();
		for (int64_t icon = 0; icon < ncon; icon++) {
			int64_t id = cel->ConnectIndex(icon);
			TPZConnect &df = fCmesh->ConnectVec()[id];
			int64_t conid = df.SequenceNumber();
			if (df.NElConnected() == 0 || conid < 0
					|| fCmesh->Block().Size(conid) == 0)
				continue;
			else {
				int64_t pos = fCmesh->Block().Position(conid);
				int64_t nsize = fCmesh->Block().Size(conid);
				for (int64_t isize = 0; isize < nsize; isize++) {
					if (isize % 2 == 0) {
						indexes[cont1] = pos + isize;
						cont1++;
					} else {
						indexes[cont2 + nf_tot] = pos + isize;
						cont2++;
					}
				}
			}
		}
	}
	this->SetIndexes(indexes);
	this->SetWeightVector(weight);
	this->ColoringElements();
	this->AssembleRhsBoundary();
#ifdef USING_CSRMV_MULT
	this->CSRInfo();
#endif
	this->TransferDataStructure();
}

void TPZIntPointsFEM::ColoringElements() const {
	int64_t nelem_c = fCmesh->NElements();
	int64_t nconnects = fCmesh->NConnects();
	TPZVec<int64_t> connects_vec(nconnects, 0);

	int64_t contcolor = 0;
	bool needstocontinue = true;

	while (needstocontinue) {
		int it = 0;
		needstocontinue = false;
		for (int64_t iel = 0; iel < nelem_c; iel++) {
			TPZCompEl *cel = fCmesh->Element(iel);
			if (!cel || cel->Dimension() != fCmesh->Dimension())
				continue;

			it++;
			if (fElemColor[it - 1] != -1)
				continue;

			TPZStack<int64_t> connectlist;
			fCmesh->Element(iel)->BuildConnectList(connectlist);
			int64_t ncon = connectlist.size();

			int64_t icon;
			for (icon = 0; icon < ncon; icon++) {
				if (connects_vec[connectlist[icon]] != 0)
					break;
			}
			if (icon != ncon) {
				needstocontinue = true;
				continue;
			}
			fElemColor[it - 1] = contcolor;
//            cel->Reference()->SetMaterialId(contcolor);

			for (icon = 0; icon < ncon; icon++) {
				connects_vec[connectlist[icon]] = 1;
			}
		}
		contcolor++;
		connects_vec.Fill(0);
	}
//    ofstream file("colored.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(fCmesh->Reference(),file);

	int64_t nelem = fRowSizes.size();
	int64_t neq = fCmesh->NEquations();
	for (int64_t iel = 0; iel < nelem; iel++) {
		int64_t cols = fColSizes[iel];
		int64_t cont_cols = fColFirstIndex[iel];

		for (int64_t icols = 0; icols < cols; icols++) {
			fIndexesColor[cont_cols + icols] = fIndexes[cont_cols + icols]
					+ fElemColor[iel] * neq;
			fIndexesColor[cont_cols + fNphis + icols] = fIndexes[cont_cols
					+ fNphis + icols] + fElemColor[iel] * neq;
		}
	}
}

void TPZIntPointsFEM::AssembleRhsBoundary() {
	int64_t neq = fCmesh->NEquations();
	fRhsBoundary.Resize(neq, 1);
	fRhsBoundary.Zero();

	for (auto iel : fBoundaryElements) {
		TPZCompEl *cel = fCmesh->Element(iel);
		if (!cel)
			continue;
		TPZElementMatrix ef(fCmesh, TPZElementMatrix::EF);
		cel->CalcResidual(ef);
		ef.ComputeDestinationIndices();
		fRhsBoundary.AddFel(ef.fMat, ef.fSourceIndex, ef.fDestinationIndex);
	}
}
